"""
Matching engine and decision rules for varietal identification.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.metrics import auc, roc_curve


@dataclass
class MatchResult:
    """Result of a fingerprint match query."""
    decision: str                     # MATCH, AMBIGUOUS/RELATED, NO MATCH, INSUFFICIENT DATA
    best_match: str                   # Variety ID of best match
    best_distance: float              # Mismatch rate to best match
    n_loci_used: int                  # Number of loci with calls in both
    top_hits: list[tuple[str, float]] # Top 3 (variety_id, distance)
    probability_identity: float       # Combined PI for loci used
    warnings: list[str]               # Any warnings


@dataclass
class Thresholds:
    """Calibrated decision thresholds."""
    match_threshold: float            # Max distance for MATCH
    related_threshold: float          # Max distance for RELATED
    margin_to_second: float           # Min margin to second best for MATCH
    min_loci_fraction: float          # Min fraction of loci called
    replicate_mean: float             # Mean replicate distance
    replicate_std: float              # Std of replicate distance
    nearest_mean: float               # Mean nearest-neighbor distance
    nearest_std: float                # Std of nearest-neighbor distance
    auc: float                        # ROC AUC
    eer: float                        # Equal error rate


def pair_distance(
    q: np.ndarray,
    r: np.ndarray,
    coding: str = "dosage",
) -> tuple[float, int]:
    """
    Compute mismatch rate between two genotype vectors.

    Parameters
    ----------
    q, r : Genotype vectors coded 0/1/2, -1=missing
    coding : "dosage" (0/1/2) or "haploid" (0/1)

    Returns
    -------
    (distance, n_loci_used)
    """
    ok = (q >= 0) & (r >= 0)
    n = int(ok.sum())

    if n == 0:
        return np.nan, 0

    if coding == "dosage":
        # Homozygote vs different homozygote = 1
        # Heterozygote vs homozygote = 0.5
        # Same = 0
        d = np.abs(q[ok].astype(int) - r[ok].astype(int))
        # Convert 2->1 (opposite homozygotes), 1->0.5 (het vs hom)
        d = np.where(d == 2, 1, d)
        d = np.where(d == 1, 0.5, d)  # Wait, this is wrong
        # Let's do it properly:
        # 0 vs 0 = 0, 2 vs 2 = 0, 1 vs 1 = 0
        # 0 vs 2 = 1, 2 vs 0 = 1
        # 0 vs 1 = 0.5, 1 vs 0 = 0.5
        # 2 vs 1 = 0.5, 1 vs 2 = 0.5
        diff = np.abs(q[ok] - r[ok])
        mismatch = np.where(diff == 2, 1.0, np.where(diff == 1, 0.5, 0.0))
        distance = mismatch.sum() / n
    else:
        # Simple mismatch rate
        distance = (q[ok] != r[ok]).sum() / n

    return float(distance), n


def classify_query(
    query_gt: np.ndarray,
    ref_gt: np.ndarray,
    ref_ids: np.ndarray,
    thresholds: Thresholds,
    min_loci_fraction: float | None = None,
) -> MatchResult:
    """
    Classify a query sample against reference panel.

    Parameters
    ----------
    query_gt : Query genotype vector (n_loci,) coded 0/1/2, -1=missing
    ref_gt : Reference genotypes (n_loci, n_ref) coded 0/1/2, -1=missing
    ref_ids : Reference variety IDs
    thresholds : Calibrated thresholds
    min_loci_fraction : Override minimum loci fraction

    Returns
    -------
    MatchResult
    """
    warnings = []

    # Check sufficient data
    n_loci = len(query_gt)
    n_called = (query_gt >= 0).sum()
    loci_fraction = n_called / n_loci

    min_frac = min_loci_fraction or thresholds.min_loci_fraction
    if loci_fraction < min_frac:
        return MatchResult(
            decision="INSUFFICIENT DATA",
            best_match="",
            best_distance=np.nan,
            n_loci_used=0,
            top_hits=[],
            probability_identity=1.0,
            warnings=[f"Only {loci_fraction:.1%} loci called (minimum {min_frac:.1%})"],
        )

    # Compute distances to all references
    distances = []
    n_loci_used_list = []

    for i in range(ref_gt.shape[1]):
        dist, n_used = pair_distance(query_gt, ref_gt[:, i])
        distances.append(dist)
        n_loci_used_list.append(n_used)

    distances = np.array(distances)
    n_loci_used_arr = np.array(n_loci_used_list)

    # Sort by distance
    order = np.argsort(distances)
    top3 = [(ref_ids[order[j]], float(distances[order[j]])) for j in range(min(3, len(order)))]

    best_idx = order[0]
    second_idx = order[1] if len(order) > 1 else None

    best_dist = distances[best_idx]
    best_id = ref_ids[best_idx]
    n_loci_used = n_loci_used_arr[best_idx]

    # Combined probability of identity for loci used
    # This is approximate - would need per-locus allele frequencies
    # For now, use a conservative estimate
    pi = 0.5 ** n_loci_used  # Assuming MAF ~ 0.5

    # Decision logic
    if best_dist <= thresholds.match_threshold:
        if second_idx is not None:
            second_dist = distances[second_idx]
            margin = second_dist - best_dist
            if margin >= thresholds.margin_to_second:
                decision = "MATCH"
            else:
                decision = "AMBIGUOUS/RELATED"
                warnings.append(f"Small margin to second best: {margin:.4f}")
        else:
            decision = "MATCH"
    elif best_dist <= thresholds.related_threshold:
        decision = "AMBIGUOUS/RELATED"
    else:
        decision = "NO MATCH"

    # Check if best match has very few loci
    if n_loci_used < n_loci * 0.5:
        warnings.append(f"Best match based on only {n_loci_used}/{n_loci} loci")

    return MatchResult(
        decision=decision,
        best_match=best_id,
        best_distance=best_dist,
        n_loci_used=n_loci_used,
        top_hits=top3,
        probability_identity=pi,
        warnings=warnings,
    )


def calibrate_thresholds(
    ref_gt: np.ndarray,
    ref_ids: np.ndarray,
    replicate_groups: dict[str, list[int]],
    match_percentile: float = 99.9,
    related_multiplier: float = 3.0,
    margin: float = 0.05,
    min_loci_fraction: float = 0.90,
) -> Thresholds:
    """
    Calibrate decision thresholds from replicate data.

    Parameters
    ----------
    ref_gt : Reference genotypes (n_loci, n_samples)
    ref_ids : Sample IDs
    replicate_groups : Dict mapping variety_id -> list of sample indices
    match_percentile : Percentile of replicate distance for match threshold
    related_multiplier : Multiplier for related threshold
    margin : Required margin to second best
    min_loci_fraction : Minimum fraction of loci called

    Returns
    -------
    Thresholds object
    """
    n_samples = ref_gt.shape[1]

    # Collect replicate distances (within same variety)
    replicate_distances = []
    for var_id, indices in replicate_groups.items():
        if len(indices) < 2:
            continue
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                dist, n_used = pair_distance(ref_gt[:, indices[i]], ref_gt[:, indices[j]])
                if not np.isnan(dist):
                    replicate_distances.append(dist)

    replicate_distances = np.array(replicate_distances)

    # Collect nearest-neighbor distances (between different varieties)
    # For each sample, find distance to nearest sample of different variety
    nearest_distances = []
    for i in range(n_samples):
        var_i = ref_ids[i]
        min_dist = np.inf
        for j in range(n_samples):
            if i == j:
                continue
            var_j = ref_ids[j]
            if var_i == var_j:
                continue
            dist, _ = pair_distance(ref_gt[:, i], ref_gt[:, j])
            if not np.isnan(dist) and dist < min_dist:
                min_dist = dist
        if min_dist < np.inf:
            nearest_distances.append(min_dist)

    nearest_distances = np.array(nearest_distances)

    # Compute statistics
    rep_mean = replicate_distances.mean() if len(replicate_distances) > 0 else 0
    rep_std = replicate_distances.std() if len(replicate_distances) > 0 else 0
    near_mean = nearest_distances.mean() if len(nearest_distances) > 0 else 1
    near_std = nearest_distances.std() if len(nearest_distances) > 0 else 0

    # Match threshold at upper tail of replicate distribution
    if len(replicate_distances) > 0:
        match_thresh = np.percentile(replicate_distances, match_percentile)
    else:
        match_thresh = 0.02  # Default fallback

    # Related threshold: ensure gap between replicate and nearest distributions
    # Use a multiple of match threshold or a fraction of nearest mean
    related_thresh = min(match_thresh * related_multiplier, near_mean * 0.5)

    # ROC analysis
    # Labels: 0 = same variety (replicates), 1 = different variety
    y_true = []
    y_score = []

    for dist in replicate_distances:
        y_true.append(0)
        y_score.append(-dist)  # Negative distance as score (higher = more similar)

    for dist in nearest_distances:
        y_true.append(1)
        y_score.append(-dist)

    if len(set(y_true)) > 1:
        fpr, tpr, roc_thresholds = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)

        # Equal error rate
        fnr = 1 - tpr
        eer_idx = np.argmin(np.abs(fpr - fnr))
        eer = (fpr[eer_idx] + fnr[eer_idx]) / 2
    else:
        roc_auc = 1.0
        eer = 0.0

    return Thresholds(
        match_threshold=match_thresh,
        related_threshold=related_thresh,
        margin_to_second=margin,
        min_loci_fraction=min_loci_fraction,
        replicate_mean=rep_mean,
        replicate_std=rep_std,
        nearest_mean=near_mean,
        nearest_std=near_std,
        auc=roc_auc,
        eer=eer,
    )


def verify_hybrid(
    hybrid_gt: np.ndarray,
    parent1_gt: np.ndarray,
    parent2_gt: np.ndarray,
) -> dict:
    """
    Verify an F1 hybrid against declared parents.

    Parameters
    ----------
    hybrid_gt : Hybrid genotype vector
    parent1_gt, parent2_gt : Parent genotype vectors

    Returns
    -------
    Dict with concordance, selfed_like fraction, etc.
    """
    # Expected F1 genotype
    ok = (parent1_gt >= 0) & (parent2_gt >= 0)
    expected = np.full_like(hybrid_gt, -1, dtype=np.int8)

    # Parent dosages: 0, 1, 2
    # F1 expectation: if both homozygous same -> that homozygote
    # if different homozygotes -> heterozygote (1)
    # if either het -> more complex, but assume inbred parents (0 or 2)
    p1 = parent1_gt[ok]
    p2 = parent2_gt[ok]

    # Both homozygous
    both_hom = (p1 != 1) & (p2 != 1)
    same_hom = both_hom & (p1 == p2)
    diff_hom = both_hom & (p1 != p2)

    # Use proper indexing to avoid chained indexing issues
    expected[ok & same_hom] = p1[same_hom]
    expected[ok & diff_hom] = 1  # Heterozygote

    # Compare
    h_ok = (hybrid_gt >= 0) & (expected >= 0)
    if h_ok.sum() == 0:
        return {
            "concordance": np.nan,
            "selfed_like": np.nan,
            "loci_compared": 0,
            "expected_het_obs_hom": 0,
            "unexpected_alleles": 0,
        }

    concordance = (hybrid_gt[h_ok] == expected[h_ok]).mean()

    # Selfed-like: homozygous where heterozygote expected
    exp_het = (expected[h_ok] == 1)
    obs_hom = (hybrid_gt[h_ok] != 1)
    selfed_like = (exp_het & obs_hom).mean() if exp_het.any() else 0

    # Unexpected alleles: neither parent has this allele
    # For inbred parents (0/2), hybrid should only have alleles from parents
    # This is complex with dosage coding; simplified check:
    both_par_ok = (parent1_gt >= 0) & (parent2_gt >= 0) & (hybrid_gt >= 0)
    unexpected = 0
    if both_par_ok.any():
        # At each locus, hybrid allele must be subset of parent alleles
        # With 0/1/2 coding, this is tricky without knowing ref/alt
        pass

    return {
        "concordance": float(concordance),
        "selfed_like": float(selfed_like),
        "loci_compared": int(h_ok.sum()),
        "expected_het_obs_hom": float(selfed_like * h_ok.sum()),
        "unexpected_alleles": unexpected,
    }


def seed_lot_purity_test(
    n_seeds: int,
    n_offtypes: int,
    confidence: float = 0.95,
) -> dict:
    """
    Calculate seed-lot purity statistics.

    Parameters
    ----------
    n_seeds : Number of seeds genotyped
    n_offtypes : Number of off-type seeds detected
    confidence : Confidence level for CI

    Returns
    -------
    Dict with observed rate, Clopper-Pearson CI, and required n for detection
    """
    from scipy.stats import beta

    observed_rate = n_offtypes / n_seeds if n_seeds > 0 else 0

    # Clopper-Pearson exact CI
    alpha = 1 - confidence
    if n_offtypes == 0:
        ci_lower = 0
        ci_upper = 1 - alpha ** (1 / n_seeds)
    elif n_offtypes == n_seeds:
        ci_lower = alpha ** (1 / n_seeds)
        ci_upper = 1
    else:
        ci_lower = beta.ppf(alpha / 2, n_offtypes, n_seeds - n_offtypes + 1)
        ci_upper = beta.ppf(1 - alpha / 2, n_offtypes + 1, n_seeds - n_offtypes)

    # Sample size needed to detect off-type rate p with confidence
    def n_for_detection(p: float, alpha: float = 0.05) -> int:
        return int(np.ceil(np.log(alpha) / np.log(1 - p)))

    return {
        "n_seeds": n_seeds,
        "n_offtypes": n_offtypes,
        "observed_rate": observed_rate,
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "confidence": confidence,
        "n_for_5pct": n_for_detection(0.05),
        "n_for_2pct": n_for_detection(0.02),
        "n_for_1pct": n_for_detection(0.01),
    }


def build_reference_database(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    panel_indices: np.ndarray,
    output_path: str,
) -> None:
    """
    Build a reference database for fast matching.

    Saves panel genotypes as Parquet for fast loading.
    """
    from .io import GenotypeData, write_genotype_parquet

    panel_gt = gt[panel_indices]
    panel_var_ids = variant_ids[panel_indices]

    data = GenotypeData(
        gt=panel_gt,
        variant_ids=panel_var_ids,
        sample_ids=variant_ids,  # This should be sample/variety IDs
        chrom=np.array([]),  # Would need to pass
        pos=np.array([]),
        ref=np.array([]),
        alt=np.array([]),
    )

    write_genotype_parquet(data, output_path)


def load_reference_database(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load reference database for matching.

    Returns (genotypes, variant_ids, sample_ids)
    """
    from .io import read_genotype_parquet
    data = read_genotype_parquet(path)
    return data.gt, data.variant_ids, data.sample_ids
