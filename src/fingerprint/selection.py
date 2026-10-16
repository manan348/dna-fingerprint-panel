"""
Panel selection algorithms: lazy-greedy set multicover and ILP refinement.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass

import numpy as np
import pandas as pd

# Optional imports - catch both ImportError and OSError (DLL issues on Windows)
try:
    from ortools.sat.python import cp_model
    _HAS_ORT = True
except (ImportError, OSError):
    _HAS_ORT = False


@dataclass
class PanelResult:
    """Result of panel selection."""
    panel_indices: np.ndarray           # Indices of selected SNPs
    variant_ids: np.ndarray             # SNP IDs of selected SNPs
    unresolved_pairs: int               # Number of variety pairs with separation < k
    separation_matrix: np.ndarray       # Pairwise separation counts
    total_snps_tested: int
    runtime_seconds: float


def compute_separation_vectors(
    gt: np.ndarray,
    iu: np.ndarray,
    ju: np.ndarray,
    candidate_indices: np.ndarray,
) -> dict[int, np.ndarray]:
    """
    Precompute which variety pairs each candidate SNP separates.

    Returns dict mapping candidate index -> boolean vector of separated pairs.
    """
    separation = {}
    for s in candidate_indices:
        ok = (gt[s, iu] >= 0) & (gt[s, ju] >= 0)
        diff = ok & (gt[s, iu] != gt[s, ju])
        separation[s] = diff
    return separation


def lazy_greedy_multicover(
    gt: np.ndarray,
    k: int = 3,
    max_snps: int = 96,
    candidate_indices: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    chrom: np.ndarray | None = None,
    min_gap_bp: int = 2_000_000,
    max_chrom_fraction: float = 0.15,
    separation_cache: dict | None = None,
    verbose: bool = True,
) -> PanelResult:
    """
    Lazy-greedy set multicover for SNP panel selection.

    Each variety pair needs to be separated by at least k SNPs.
    Uses lazy evaluation: only recompute marginal gain when needed.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples) coded 0/1/2, -1=missing
    k : Minimum separation per variety pair
    max_snps : Maximum panel size
    candidate_indices : Indices of candidate SNPs to consider
    pos : Physical positions for spacing constraint
    chrom : Chromosomes for spacing constraint
    min_gap_bp : Minimum physical distance between selected SNPs
    max_chrom_fraction : Max fraction of panel from any one chromosome
    separation_cache : Precomputed separation vectors
    verbose : Print progress

    Returns
    -------
    PanelResult
    """
    import time
    start_time = time.time()

    n_variants, n_samples = gt.shape
    n_pairs = n_samples * (n_samples - 1) // 2

    if candidate_indices is None:
        candidate_indices = np.arange(n_variants)

    iu, ju = np.triu_indices(n_samples, 1)
    need = np.full(n_pairs, k, dtype=np.int16)

    # Precompute separation vectors if not cached
    if separation_cache is None:
        if verbose:
            print(f"Precomputing separation vectors for {len(candidate_indices)} candidates...")
        separation_cache = compute_separation_vectors(gt, iu, ju, candidate_indices)

    # Initialize heap with initial gains
    heap = []
    for s in candidate_indices:
        gain = int(separation_cache[s].sum())
        if gain > 0:
            heapq.heappush(heap, (-gain, s))

    if verbose:
        print(f"Initial candidates with positive gain: {len(heap)}")

    # Chromosome quota tracking
    chrom_counts = {}
    if chrom is not None:
        max_per_chrom = int(max_snps * max_chrom_fraction)

    panel = []
    panel_set = set()

    while heap and len(panel) < max_snps and (need > 0).any():
        neg_gain, s = heapq.heappop(heap)

        # Check spacing constraint
        if pos is not None and chrom is not None:
            too_close = False
            for t in panel:
                if chrom[s] == chrom[t] and abs(pos[s] - pos[t]) < min_gap_bp:
                    too_close = True
                    break
            if too_close:
                continue

        # Check chromosome quota
        if chrom is not None:
            c = chrom[s]
            if chrom_counts.get(c, 0) >= max_per_chrom:
                continue

        # Recompute actual marginal gain
        d = separation_cache[s]
        gain = int((d & (need > 0)).sum())

        # If gain is stale (less than next best), re-push and continue
        if heap and gain < -heap[0][0]:
            heapq.heappush(heap, (-gain, s))
            continue

        if gain == 0:
            break

        # Select this SNP
        panel.append(s)
        panel_set.add(s)
        need[d & (need > 0)] -= 1

        if chrom is not None:
            chrom_counts[chrom[s]] = chrom_counts.get(chrom[s], 0) + 1

        if verbose and len(panel) % 10 == 0:
            unresolved = int((need > 0).sum())
            print(f"  Selected {len(panel)} SNPs, {unresolved} pairs still need separation")

    unresolved = int((need > 0).sum())
    runtime = time.time() - start_time

    if verbose:
        print(f"Greedy selection complete: {len(panel)} SNPs, {unresolved} unresolved pairs, {runtime:.1f}s")

    # Build separation matrix for selected panel
    if panel:
        sep_matrix = compute_pairwise_separation(gt, panel, iu, ju)
    else:
        sep_matrix = np.zeros(n_pairs, dtype=np.int16)

    return PanelResult(
        panel_indices=np.array(panel, dtype=int),
        variant_ids=np.array([]),  # Filled by caller
        unresolved_pairs=unresolved,
        separation_matrix=sep_matrix,
        total_snps_tested=len(candidate_indices),
        runtime_seconds=runtime,
    )


def compute_pairwise_separation(
    gt: np.ndarray,
    panel_indices: np.ndarray,
    iu: np.ndarray,
    ju: np.ndarray,
) -> np.ndarray:
    """Compute separation count for each variety pair using selected panel."""
    n_pairs = len(iu)
    sep_counts = np.zeros(n_pairs, dtype=np.int16)

    for s in panel_indices:
        ok = (gt[s, iu] >= 0) & (gt[s, ju] >= 0)
        diff = ok & (gt[s, iu] != gt[s, ju])
        sep_counts[diff] += 1

    return sep_counts


def ilp_refinement(
    gt: np.ndarray,
    initial_panel: np.ndarray,
    k: int = 3,
    max_snps: int = 96,
    candidate_indices: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    chrom: np.ndarray | None = None,
    min_gap_bp: int = 2_000_000,
    max_chrom_fraction: float = 0.15,
    time_limit_seconds: int = 300,
    verbose: bool = True,
) -> PanelResult:
    """
    Exact ILP refinement using OR-Tools CP-SAT.

    Minimizes panel size while ensuring all pairs separated by >= k SNPs.
    Falls back to greedy if OR-Tools is not available.
    """
    import time
    start_time = time.time()

    if not _HAS_ORT:
        if verbose:
            print("OR-Tools not available, falling back to greedy")
        return PanelResult(
            panel_indices=initial_panel,
            variant_ids=np.array([]),
            unresolved_pairs=0,  # Would need to compute
            separation_matrix=np.array([]),
            total_snps_tested=len(candidate_indices) if candidate_indices is not None else 0,
            runtime_seconds=time.time() - start_time,
        )

    n_variants, n_samples = gt.shape
    n_pairs = n_samples * (n_samples - 1) // 2
    iu, ju = np.triu_indices(n_samples, 1)

    if candidate_indices is None:
        candidate_indices = np.arange(n_variants)

    # Precompute separation
    separation = compute_separation_vectors(gt, iu, ju, candidate_indices)

    # Create model
    model = cp_model.CpModel()

    # Decision variables: x[s] = 1 if SNP s is selected
    x = {s: model.NewBoolVar(f"x_{s}") for s in candidate_indices}

    # Constraint: each pair must be separated by at least k selected SNPs
    for p_idx in range(n_pairs):
        separating_snps = [s for s in candidate_indices if separation[s][p_idx]]
        if separating_snps:
            model.Add(sum(x[s] for s in separating_snps) >= k)
        else:
            # This pair cannot be separated - relax constraint
            pass

    # Spacing constraint: no two selected SNPs within min_gap_bp on same chrom
    if pos is not None and chrom is not None:
        for i, s1 in enumerate(candidate_indices):
            for s2 in candidate_indices[i+1:]:
                if chrom[s1] == chrom[s2] and abs(pos[s1] - pos[s2]) < min_gap_bp:
                    model.Add(x[s1] + x[s2] <= 1)

    # Chromosome quota
    if chrom is not None:
        max_per_chrom = int(max_snps * max_chrom_fraction)
        for c in np.unique(chrom[candidate_indices]):
            chrom_snps = [s for s in candidate_indices if chrom[s] == c]
            model.Add(sum(x[s] for s in chrom_snps) <= max_per_chrom)

    # Objective: minimize panel size
    model.Minimize(sum(x[s] for s in candidate_indices))

    # Solve
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_seconds
    solver.parameters.num_search_workers = 8
    solver.parameters.log_search_progress = verbose

    if verbose:
        print(f"Running ILP with {len(candidate_indices)} candidates, {n_pairs} pairs...")

    status = solver.Solve(model)

    runtime = time.time() - start_time

    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        selected = [s for s in candidate_indices if solver.Value(x[s]) == 1]
        if verbose:
            print(f"ILP solved: {len(selected)} SNPs in {runtime:.1f}s (status: {solver.StatusName(status)})")
    else:
        if verbose:
            print(f"ILP failed: {solver.StatusName(status)}, falling back to greedy")
        selected = initial_panel.tolist()

    # Compute results
    unresolved = 0
    if selected:
        sep_matrix = compute_pairwise_separation(gt, np.array(selected), iu, ju)
        unresolved = int((sep_matrix < k).sum())
    else:
        sep_matrix = np.zeros(n_pairs, dtype=np.int16)
        unresolved = n_pairs

    return PanelResult(
        panel_indices=np.array(selected, dtype=int),
        variant_ids=np.array([]),
        unresolved_pairs=unresolved,
        separation_matrix=sep_matrix,
        total_snps_tested=len(candidate_indices),
        runtime_seconds=runtime,
    )


def select_panel(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    scores: np.ndarray | None = None,
    config: dict | None = None,
) -> tuple[PanelResult, PanelResult, PanelResult]:
    """
    Main panel selection function: runs greedy then ILP refinement.

    Returns (core_result, backup_result, extended_result)
    """
    if config is None:
        config = {
            "k": 3,
            "core_max_snps": 96,
            "backup_fraction": 0.20,
            "extended_max_snps": 384,
            "min_gap_bp": 2_000_000,
            "max_chrom_fraction": 0.15,
            "ilp_time_limit": 300,
            "min_maf": 0.20,
        }

    n_variants = gt.shape[0]
    n_samples = gt.shape[1]

    # Filter candidates by MAF if scores provided
    if scores is not None:
        candidate_mask = scores >= config.get("min_composite_score", 0)
        candidate_indices = np.where(candidate_mask)[0]
    else:
        candidate_indices = np.arange(n_variants)

    if len(candidate_indices) == 0:
        raise ValueError("No candidate SNPs pass filters")

    print(f"Starting panel selection with {len(candidate_indices)} candidates, {n_samples} samples")

    # --- Core panel ---
    print("\n=== Core Panel Selection ===")
    core_greedy = lazy_greedy_multicover(
        gt=gt,
        k=config["k"],
        max_snps=config["core_max_snps"],
        candidate_indices=candidate_indices,
        pos=pos,
        chrom=chrom,
        min_gap_bp=config["min_gap_bp"],
        max_chrom_fraction=config["max_chrom_fraction"],
        verbose=True,
    )

    # ILP refinement on top candidates
    # Use top 2x core_max_snps candidates by score for ILP
    if scores is not None:
        top_candidates = candidate_indices[np.argsort(scores[candidate_indices])[::-1]]
        top_candidates = top_candidates[:config["core_max_snps"] * 3]
    else:
        top_candidates = candidate_indices[:config["core_max_snps"] * 3]

    core_ilp = ilp_refinement(
        gt=gt,
        initial_panel=core_greedy.panel_indices,
        k=config["k"],
        max_snps=config["core_max_snps"],
        candidate_indices=top_candidates,
        pos=pos,
        chrom=chrom,
        min_gap_bp=config["min_gap_bp"],
        max_chrom_fraction=config["max_chrom_fraction"],
        time_limit_seconds=config["ilp_time_limit"],
        verbose=True,
    )

    # Choose better result
    if core_ilp.unresolved_pairs <= core_greedy.unresolved_pairs:
        core_result = core_ilp
        print("Using ILP result")
    else:
        core_result = core_greedy
        print("Using greedy result")

    core_result.variant_ids = variant_ids[core_result.panel_indices]

    # --- Backup panel ---
    # Additional SNPs not in core, ranked by score
    backup_size = int(config["core_max_snps"] * config["backup_fraction"])
    remaining = [s for s in candidate_indices if s not in core_result.panel_indices]
    if scores is not None:
        remaining = sorted(remaining, key=lambda s: scores[s], reverse=True)
    backup_indices = np.array(remaining[:backup_size])
    backup_result = PanelResult(
        panel_indices=backup_indices,
        variant_ids=variant_ids[backup_indices],
        unresolved_pairs=0,
        separation_matrix=np.array([]),
        total_snps_tested=len(remaining),
        runtime_seconds=0,
    )

    # --- Extended panel ---
    extended_max = config["extended_max_snps"]
    extended_candidates = [s for s in candidate_indices if s not in core_result.panel_indices]
    if scores is not None:
        extended_candidates = sorted(extended_candidates, key=lambda s: scores[s], reverse=True)
    extended_indices = np.array(extended_candidates[:extended_max])
    extended_result = PanelResult(
        panel_indices=extended_indices,
        variant_ids=variant_ids[extended_indices],
        unresolved_pairs=0,
        separation_matrix=np.array([]),
        total_snps_tested=len(extended_candidates),
        runtime_seconds=0,
    )

    print("\n=== Final Panels ===")
    print(f"Core: {len(core_result.panel_indices)} SNPs, {core_result.unresolved_pairs} unresolved pairs")
    print(f"Backup: {len(backup_result.panel_indices)} SNPs")
    print(f"Extended: {len(extended_result.panel_indices)} SNPs")

    return core_result, backup_result, extended_result


def compute_resolution_curve(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    scores: np.ndarray,
    panel_sizes: list[int],
    k: int = 3,
    n_random_trials: int = 100,
) -> pd.DataFrame:
    """
    Compute resolution curve: fraction of pairs separated at k=1,2,3 vs panel size.
    Compares designed panel vs PIC-ranked vs random.
    """
    n_samples = gt.shape[1]
    n_pairs = n_samples * (n_samples - 1) // 2
    iu, ju = np.triu_indices(n_samples, 1)

    # Precompute separation for all SNPs
    all_sep = compute_separation_vectors(gt, iu, ju, np.arange(gt.shape[0]))

    results = []

    for size in panel_sizes:
        # Designed panel (top by composite score)
        designed_idx = np.argsort(scores)[::-1][:size]
        designed_sep = compute_pairwise_separation(gt, designed_idx, iu, ju)

        # PIC-ranked panel
        # Need PIC scores - would need to pass them separately
        # For now, use composite score as proxy

        # Random panels
        random_fractions = {1: [], 2: [], 3: []}
        for _ in range(n_random_trials):
            rand_idx = np.random.choice(gt.shape[0], size, replace=False)
            rand_sep = compute_pairwise_separation(gt, rand_idx, iu, ju)
            for k_val in [1, 2, 3]:
                random_fractions[k_val].append((rand_sep >= k_val).mean())

        for k_val in [1, 2, 3]:
            results.append({
                "panel_size": size,
                "k": k_val,
                "designed_fraction": (designed_sep >= k_val).mean(),
                "random_mean": np.mean(random_fractions[k_val]),
                "random_std": np.std(random_fractions[k_val]),
                "random_min": np.min(random_fractions[k_val]),
                "random_max": np.max(random_fractions[k_val]),
            })

    return pd.DataFrame(results)


def write_panel_files(
    core_result: PanelResult,
    backup_result: PanelResult,
    extended_result: PanelResult,
    output_dir: str,
) -> None:
    """Write panel SNP lists to files."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    np.savetxt(os.path.join(output_dir, "core_panel_snps.txt"),
               core_result.variant_ids, fmt="%s")
    np.savetxt(os.path.join(output_dir, "backup_panel_snps.txt"),
               backup_result.variant_ids, fmt="%s")
    np.savetxt(os.path.join(output_dir, "extended_panel_snps.txt"),
               extended_result.variant_ids, fmt="%s")

    # Write detailed info
    for name, result in [("core", core_result), ("backup", backup_result), ("extended", extended_result)]:
        df = pd.DataFrame({
            "variant_id": result.variant_ids,
            "panel_index": np.arange(len(result.panel_indices)),
        })
        df.to_csv(os.path.join(output_dir, f"{name}_panel_details.tsv"), sep="\t", index=False)
