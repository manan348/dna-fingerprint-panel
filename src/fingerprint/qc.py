"""
Quality control and filtering for genotype data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .io import GenotypeData


@dataclass
class QCMetrics:
    """Container for QC metrics."""
    n_samples_initial: int
    n_variants_initial: int
    n_samples_after_sample_qc: int
    n_variants_after_snp_qc: int
    sample_missing_rates: np.ndarray
    variant_missing_rates: np.ndarray
    variant_maf: np.ndarray
    variant_het: np.ndarray
    clone_groups: pd.DataFrame
    removed_samples: list[str]
    removed_variants: list[str]


def sample_qc(
    gt: np.ndarray,
    sample_ids: np.ndarray,
    max_missing: float = 0.10,
    max_het_rate: float = 0.10,
    het_warning_only: bool = True,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Perform sample-level QC.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples) coded 0/1/2, -1=missing
    sample_ids : Sample identifiers
    max_missing : Maximum missing rate per sample
    max_het_rate : Maximum heterozygosity rate for inbreds
    het_warning_only : If True, flag high-het samples but don't remove

    Returns
    -------
    (keep_mask, removed_samples, metrics_dict)
    """
    n_samples = gt.shape[1]
    missing_per_sample = (gt == -1).mean(axis=0)
    het_per_sample = (gt == 1).sum(axis=0) / np.maximum((gt >= 0).sum(axis=0), 1)

    keep = missing_per_sample <= max_missing
    removed = sample_ids[~keep].tolist()

    # Heterozygosity check
    high_het = het_per_sample > max_het_rate
    if het_warning_only:
        pass  # Just flag
    else:
        keep = keep & ~high_het
        removed.extend(sample_ids[high_het & (missing_per_sample <= max_missing)].tolist())

    metrics = {
        "n_initial": n_samples,
        "n_kept": int(keep.sum()),
        "n_removed_missing": int((missing_per_sample > max_missing).sum()),
        "n_flagged_high_het": int(high_het.sum()),
        "missing_rates": missing_per_sample,
        "het_rates": het_per_sample,
    }

    return keep, np.array(removed), metrics


def variant_qc(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    max_missing: float = 0.05,
    min_maf: float = 0.20,
    max_het: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Perform variant-level QC.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples) coded 0/1/2, -1=missing
    variant_ids : Variant identifiers
    max_missing : Maximum missing rate per variant
    min_maf : Minimum minor allele frequency
    max_het : Maximum heterozygosity rate (for inbreds)

    Returns
    -------
    (keep_mask, removed_variants, metrics_dict)
    """
    n_variants = gt.shape[0]
    n_samples = gt.shape[1]

    called = gt >= 0
    n_called = called.sum(axis=1)
    missing_rate = 1 - n_called / n_samples

    # Allele frequency
    alt_count = np.where(called, gt, 0).sum(axis=1)
    alt_freq = alt_count / (2 * np.maximum(n_called, 1))
    maf = np.minimum(alt_freq, 1 - alt_freq)

    # Heterozygosity
    het_rate = np.where(called, gt == 1, 0).sum(axis=1) / np.maximum(n_called, 1)

    keep = (missing_rate <= max_missing) & (maf >= min_maf) & (het_rate <= max_het)
    removed = variant_ids[~keep].tolist()

    metrics = {
        "n_initial": n_variants,
        "n_kept": int(keep.sum()),
        "n_removed_missing": int((missing_rate > max_missing).sum()),
        "n_removed_low_maf": int((maf < min_maf).sum()),
        "n_removed_high_het": int((het_rate > max_het).sum()),
        "missing_rates": missing_rate,
        "maf": maf,
        "het_rates": het_rate,
    }

    return keep, np.array(removed), metrics


def filter_mappability(
    chrom: np.ndarray,
    pos: np.ndarray,
    repeat_mask_bed: str,
) -> np.ndarray:
    """
    Filter variants in repeat-masked or low-mappability regions.

    Parameters
    ----------
    chrom : Chromosome array
    pos : Position array
    repeat_mask_bed : Path to BED file with repeat regions

    Returns
    -------
    Boolean mask (True = keep)
    """
    # Use bedtools or pybedtools for intersection
    # Placeholder implementation
    try:
        import pybedtools
        variants_bed = pybedtools.BedTool()
        # Create temporary BED from coordinates
        # This is simplified - real implementation needs proper BED creation
        return np.ones(len(chrom), dtype=bool)
    except ImportError:
        print("Warning: pybedtools not available, skipping mappability filter")
        return np.ones(len(chrom), dtype=bool)


def detect_clones(
    gt: np.ndarray,
    sample_ids: np.ndarray,
    min_shared_snps: int = 10000,
    identity_threshold: float = 0.99,
) -> pd.DataFrame:
    """
    Detect clonal/duplicate samples by pairwise IBS.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples) coded 0/1/2, -1=missing
    sample_ids : Sample identifiers
    min_shared_snps : Minimum shared SNPs for comparison
    identity_threshold : IBS threshold for calling clones

    Returns
    -------
    DataFrame with clone groups
    """
    n_samples = gt.shape[1]
    clone_groups = []

    # Sample a subset of SNPs for efficiency if many variants
    n_variants = gt.shape[0]
    if n_variants > min_shared_snps:
        subset = np.random.choice(n_variants, min_shared_snps, replace=False)
        gt_sub = gt[subset]
    else:
        gt_sub = gt

    # Compute pairwise IBS (identity by state)
    # For efficiency, use matrix operations
    # IBS = proportion of shared alleles among called sites
    called = gt_sub >= 0

    # Compare all pairs
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            both_called = called[:, i] & called[:, j]
            n_shared = both_called.sum()
            if n_shared < min_shared_snps:
                continue
            matches = (gt_sub[both_called, i] == gt_sub[both_called, j]).sum()
            ibs = matches / n_shared
            if ibs >= identity_threshold:
                clone_groups.append({
                    "sample1": sample_ids[i],
                    "sample2": sample_ids[j],
                    "ibs": ibs,
                    "n_shared_snps": n_shared,
                })

    if clone_groups:
        df = pd.DataFrame(clone_groups)
        # Group into clusters (simple transitive closure)
        groups = assign_clone_groups(df)
        return groups
    else:
        return pd.DataFrame(columns=["group_id", "sample", "representative"])


def assign_clone_groups(pairs_df: pd.DataFrame) -> pd.DataFrame:
    """Assign clone group IDs from pairwise clone relationships."""
    # Build adjacency
    from collections import defaultdict

    adj = defaultdict(set)
    for _, row in pairs_df.iterrows():
        adj[row["sample1"]].add(row["sample2"])
        adj[row["sample2"]].add(row["sample1"])

    visited = set()
    groups = []
    group_id = 0

    for sample in adj:
        if sample in visited:
            continue
        # BFS to find connected component
        component = set()
        stack = [sample]
        while stack:
            s = stack.pop()
            if s in visited:
                continue
            visited.add(s)
            component.add(s)
            stack.extend(adj[s] - visited)

        # Choose representative (first alphabetically)
        rep = sorted(component)[0]
        for s in component:
            groups.append({
                "group_id": group_id,
                "sample": s,
                "representative": rep,
            })
        group_id += 1

    return pd.DataFrame(groups)


def collapse_clones(
    gt: np.ndarray,
    sample_ids: np.ndarray,
    clone_groups: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Collapse clone groups to a single representative sample.

    Uses majority vote for genotype calls within each group.
    """
    if clone_groups.empty:
        return gt, sample_ids

    # Get unique representatives
    reps = clone_groups.drop_duplicates("representative")["representative"].values
    rep_to_idx = {r: i for i, r in enumerate(sample_ids) if r in reps}

    # For each group, compute consensus
    n_variants = gt.shape[0]
    new_gt = np.zeros((n_variants, len(reps)), dtype=np.int8)
    new_ids = []

    for i, rep in enumerate(reps):
        group_samples = clone_groups[clone_groups["representative"] == rep]["sample"].values
        group_indices = [np.where(sample_ids == s)[0][0] for s in group_samples]

        # Majority vote per variant
        for v in range(n_variants):
            calls = gt[v, group_indices]
            called = calls[calls >= 0]
            if len(called) == 0:
                new_gt[v, i] = -1
            else:
                # Most common genotype
                vals, counts = np.unique(called, return_counts=True)
                new_gt[v, i] = vals[np.argmax(counts)]

        new_ids.append(rep)

    return new_gt, np.array(new_ids)


def run_full_qc(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    sample_ids: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    config: dict,
) -> tuple[GenotypeData, QCMetrics]:
    """
    Run complete QC pipeline.

    Parameters
    ----------
    gt, variant_ids, sample_ids, chrom, pos : Genotype data arrays
    config : Dictionary with QC parameters

    Returns
    -------
    (filtered_GenotypeData, QCMetrics)
    """
    # Sample QC
    sample_keep, removed_samples, sample_metrics = sample_qc(
        gt, sample_ids,
        max_missing=config.get("sample_missing_threshold", 0.10),
        max_het_rate=config.get("het_threshold_inbreds", 0.10),
        het_warning_only=True,
    )

    gt = gt[:, sample_keep]
    sample_ids = sample_ids[sample_keep]

    # Variant QC
    variant_keep, removed_variants, variant_metrics = variant_qc(
        gt, variant_ids,
        max_missing=config.get("snp_missing_threshold", 0.05),
        min_maf=config.get("maf_threshold_fingerprint", 0.20),
        max_het=config.get("het_threshold_inbreds", 0.05),
    )

    gt = gt[variant_keep]
    variant_ids = variant_ids[variant_keep]
    chrom = chrom[variant_keep]
    pos = pos[variant_keep]

    # Mappability filter
    if "repeat_mask_bed" in config:
        mappability_keep = filter_mappability(chrom, pos, config["repeat_mask_bed"])
        gt = gt[mappability_keep]
        variant_ids = variant_ids[mappability_keep]
        chrom = chrom[mappability_keep]
        pos = pos[mappability_keep]
    else:
        mappability_keep = np.ones(len(variant_ids), dtype=bool)

    # Clone detection
    clone_groups = detect_clones(
        gt, sample_ids,
        min_shared_snps=config.get("clone_min_shared_snps", 10000),
        identity_threshold=config.get("clone_identity_threshold", 0.99),
    )

    # Collapse clones
    gt_collapsed, sample_ids_collapsed = collapse_clones(gt, sample_ids, clone_groups)

    # Build metrics object
    metrics = QCMetrics(
        n_samples_initial=len(sample_ids),
        n_variants_initial=len(variant_ids),
        n_samples_after_sample_qc=int(sample_keep.sum()),
        n_variants_after_snp_qc=int(variant_keep.sum() & mappability_keep.sum()),
        sample_missing_rates=sample_metrics["missing_rates"],
        variant_missing_rates=variant_metrics["missing_rates"],
        variant_maf=variant_metrics["maf"],
        variant_het=variant_metrics["het_rates"],
        clone_groups=clone_groups,
        removed_samples=removed_samples,
        removed_variants=removed_variants.tolist(),
    )

    # Create filtered GenotypeData
    filtered_data = GenotypeData(
        gt=gt_collapsed,
        variant_ids=variant_ids[mappability_keep] if "repeat_mask_bed" in config else variant_ids,
        sample_ids=sample_ids_collapsed,
        chrom=chrom[mappability_keep] if "repeat_mask_bed" in config else chrom,
        pos=pos[mappability_keep] if "repeat_mask_bed" in config else pos,
        ref=np.array([]),  # Would need to pass through
        alt=np.array([]),
        metadata={"qc_config": config},
    )

    return filtered_data, metrics
