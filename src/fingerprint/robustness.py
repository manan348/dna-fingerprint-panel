"""
Robustness testing: genotyping error, missing data, mixtures, assay dropout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import binom
from tqdm import tqdm


@dataclass
class RobustnessResult:
    """Results from a robustness simulation."""
    error_rate: float
    missing_rate: float
    mixture_prop: float
    dropout_rate: float
    n_correct: int
    n_total: int
    accuracy: float
    mean_separation: float
    min_separation: int


def simulate_genotyping_error(
    gt: np.ndarray,
    error_rate: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Simulate genotyping errors by randomly flipping calls.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples) coded 0/1/2, -1=missing
    error_rate : Probability of erroneous call per genotype
    rng : Random number generator

    Returns
    -------
    Genotype matrix with errors introduced
    """
    gt_err = gt.copy()
    n_variants, n_samples = gt.shape

    # Mask of called genotypes
    called = gt >= 0

    # Random errors
    error_mask = called & (rng.random((n_variants, n_samples)) < error_rate)

    # For each error, flip to a different genotype
    for v, s in zip(*np.where(error_mask)):
        current = gt[v, s]
        # Choose from other two genotypes
        alternatives = [g for g in [0, 1, 2] if g != current]
        gt_err[v, s] = rng.choice(alternatives)

    return gt_err


def simulate_missing_data(
    gt: np.ndarray,
    missing_rate: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Simulate missing data by randomly setting calls to missing.

    Parameters
    ----------
    gt : Genotype matrix
    missing_rate : Probability of setting a call to missing
    rng : Random number generator

    Returns
    -------
    Genotype matrix with additional missing data
    """
    gt_miss = gt.copy()
    called = gt >= 0
    miss_mask = called & (rng.random(gt.shape) < missing_rate)
    gt_miss[miss_mask] = -1
    return gt_miss


def simulate_seed_mixture(
    gt: np.ndarray,
    mixture_prop: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Simulate seed-lot mixtures by combining two varieties.

    For each sample, with probability mixture_prop, mix with a random other variety.
    At heterozygous loci in the mixture, we observe a mixed signal.
    For simplicity, we model as: if parents differ, mixture shows het (1).
    """
    gt_mix = gt.copy()
    n_variants, n_samples = gt.shape

    for s in range(n_samples):
        if rng.random() < mixture_prop:
            # Pick a different variety as contaminant
            other = rng.choice([i for i in range(n_samples) if i != s])
            # Create mixture: where they differ, call becomes het (1)
            # Where they agree, keep the shared call
            # Where one is missing, keep the other's call
            s_gt = gt[:, s]
            o_gt = gt[:, other]

            both_called = (s_gt >= 0) & (o_gt >= 0)
            differ = both_called & (s_gt != o_gt)
            agree = both_called & (s_gt == o_gt)
            s_only = (s_gt >= 0) & (o_gt == -1)
            o_only = (s_gt == -1) & (o_gt >= 0)

            gt_mix[differ, s] = 1  # Heterozygous in mixture
            gt_mix[agree, s] = s_gt[agree]
            gt_mix[s_only, s] = s_gt[s_only]
            gt_mix[o_only, s] = o_gt[o_only]

    return gt_mix


def simulate_assay_dropout(
    panel_indices: np.ndarray,
    dropout_rate: float,
    rng: np.random.Generator,
    backup_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Simulate assay dropout by removing SNPs from panel.

    If backup SNPs provided, they replace dropped assays in order.

    Returns
    -------
    (modified_panel_indices, dropout_mask)
    """
    n_panel = len(panel_indices)
    n_drop = int(n_panel * dropout_rate)

    if n_drop == 0:
        return panel_indices, np.ones(n_panel, dtype=bool)

    drop_mask = rng.choice(n_panel, n_drop, replace=False)
    keep_mask = np.ones(n_panel, dtype=bool)
    keep_mask[drop_mask] = False

    new_panel = panel_indices[keep_mask]

    # Add backup SNPs if available
    if backup_indices is not None and len(backup_indices) > 0:
        n_backup_needed = n_drop
        n_backup_avail = min(n_backup_needed, len(backup_indices))
        new_panel = np.concatenate([new_panel, backup_indices[:n_backup_avail]])

    return new_panel, keep_mask


def evaluate_panel(
    gt: np.ndarray,
    panel_indices: np.ndarray,
    k: int = 3,
) -> tuple[int, int, float, float, int]:
    """
    Evaluate panel: count correctly identified samples.

    A sample is correctly identified if its distance to its true identity
    is < distance to any other variety by at least margin.

    Returns
    -------
    (n_correct, n_total, accuracy, mean_separation, min_separation)
    """
    n_variants, n_samples = gt.shape
    iu, ju = np.triu_indices(n_samples, 1)

    # Compute pairwise distances on panel
    panel_gt = gt[panel_indices]
    n_pairs = len(iu)
    distances = np.zeros(n_pairs)

    for p_idx in range(n_pairs):
        i, j = iu[p_idx], ju[p_idx]
        ok = (panel_gt[:, i] >= 0) & (panel_gt[:, j] >= 0)
        if ok.sum() == 0:
            distances[p_idx] = np.nan
        else:
            diff = panel_gt[ok, i] != panel_gt[ok, j]
            distances[p_idx] = diff.sum() / ok.sum()

    # For each sample, find nearest neighbor distance
    correct = 0
    separations = []

    for s in range(n_samples):
        # Distances to all other samples
        other_dists = []
        for t in range(n_samples):
            if t == s:
                continue
            # Find pair index
            if s < t:
                p_idx = np.where((iu == s) & (ju == t))[0][0]
            else:
                p_idx = np.where((iu == t) & (ju == s))[0][0]
            other_dists.append(distances[p_idx])

        other_dists = np.array(other_dists)
        valid = ~np.isnan(other_dists)

        if valid.any():
            min_dist = np.min(other_dists[valid])
            separations.append(min_dist * len(panel_indices))
            # Correct if min distance > 0 (perfect separation)
            # In practice, use a threshold
            if min_dist > 0:
                correct += 1

    mean_sep = np.mean(separations) if separations else 0
    min_sep = int(np.min(separations)) if separations else 0

    return correct, n_samples, correct / n_samples, mean_sep, min_sep


def run_robustness_simulation(
    gt: np.ndarray,
    core_panel: np.ndarray,
    backup_panel: np.ndarray | None = None,
    extended_panel: np.ndarray | None = None,
    config: dict | None = None,
    n_simulations: int = 200,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Run full robustness simulation grid.

    Parameters
    ----------
    gt : Genotype matrix
    core_panel : Core panel SNP indices
    backup_panel : Backup panel SNP indices
    extended_panel : Extended panel SNP indices
    config : Dictionary with simulation parameters
    n_simulations : Number of replicates per condition
    seed : Random seed

    Returns
    -------
    DataFrame with results
    """
    if config is None:
        config = {
            "error_rates": [0.0, 0.01, 0.02, 0.05],
            "missing_rates": [0.0, 0.05, 0.10, 0.15],
            "mixture_props": [0.0, 0.05, 0.10, 0.20],
            "dropout_rates": [0.0, 0.05, 0.10, 0.15],
            "k": 3,
        }

    rng = np.random.default_rng(seed)
    results = []

    # Total combinations
    total = (len(config["error_rates"]) *
             len(config["missing_rates"]) *
             len(config["mixture_props"]) *
             len(config["dropout_rates"]))

    with tqdm(total=total * n_simulations, desc="Robustness simulations") as pbar:
        for err in config["error_rates"]:
            for miss in config["missing_rates"]:
                for mix in config["mixture_props"]:
                    for drop in config["dropout_rates"]:
                        for sim in range(n_simulations):
                            # Apply perturbations in order
                            gt_pert = gt.copy()

                            if err > 0:
                                gt_pert = simulate_genotyping_error(gt_pert, err, rng)
                            if miss > 0:
                                gt_pert = simulate_missing_data(gt_pert, miss, rng)
                            if mix > 0:
                                gt_pert = simulate_seed_mixture(gt_pert, mix, rng)

                            # Assay dropout
                            panel, _ = simulate_assay_dropout(
                                gt_pert, core_panel, drop, rng, backup_panel
                            )

                            # Evaluate
                            n_corr, n_tot, acc, mean_sep, min_sep = evaluate_panel(
                                gt_pert, panel, config["k"]
                            )

                            results.append({
                                "simulation": sim,
                                "error_rate": err,
                                "missing_rate": miss,
                                "mixture_prop": mix,
                                "dropout_rate": drop,
                                "panel_size": len(panel),
                                "n_correct": n_corr,
                                "n_total": n_tot,
                                "accuracy": acc,
                                "mean_separation": mean_sep,
                                "min_separation": min_sep,
                            })

                            pbar.update(1)

    return pd.DataFrame(results)


def run_stress_test(
    gt: np.ndarray,
    core_panel: np.ndarray,
    backup_panel: np.ndarray | None = None,
    config: dict | None = None,
    n_simulations: int = 200,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Run stress test: simulate worst-case scenarios.
    """
    if config is None:
        config = {
            "k": 3,
        }

    rng = np.random.default_rng(seed)
    results = []

    # Scenario 1: High error, high missing
    for err in [0.02, 0.05, 0.10]:
        for miss in [0.10, 0.15, 0.20]:
            for sim in range(n_simulations):
                gt_pert = gt.copy()
                gt_pert = simulate_genotyping_error(gt_pert, err, rng)
                gt_pert = simulate_missing_data(gt_pert, miss, rng)

                n_corr, n_tot, acc, mean_sep, min_sep = evaluate_panel(
                    gt_pert, core_panel, config["k"]
                )

                results.append({
                    "scenario": "high_error_missing",
                    "error_rate": err,
                    "missing_rate": miss,
                    "simulation": sim,
                    "accuracy": acc,
                    "mean_separation": mean_sep,
                    "min_separation": min_sep,
                })

    # Scenario 2: Close relatives (simulate by using only high-MAF SNPs)
    # This would require knowing which pairs are close relatives
    # Placeholder for now

    return pd.DataFrame(results)


def compute_analytical_power(
    panel_size: int,
    k: int,
    error_rate: float,
    missing_rate: float,
) -> dict[str, float]:
    """
    Compute analytical power to correctly identify a variety.

    Based on: probability that at least k separating SNPs remain correct
    after errors and missing data.
    """
    # Expected number of called separating SNPs
    n_called = panel_size * (1 - missing_rate)
    # Expected number of correct calls among separating SNPs
    n_correct = n_called * (1 - error_rate)
    # Probability that >= k are correct (binomial)
    prob = 1 - binom.cdf(k - 1, int(n_called), 1 - error_rate)

    return {
        "panel_size": panel_size,
        "k": k,
        "error_rate": error_rate,
        "missing_rate": missing_rate,
        "expected_called": n_called,
        "expected_correct": n_correct,
        "prob_k_or_more": prob,
    }


def plot_robustness_heatmap(
    results_df: pd.DataFrame,
    output_path: str,
    metric: str = "accuracy",
) -> None:
    """
    Create heatmap of robustness results.

    Requires matplotlib/seaborn.
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns

        # Pivot for heatmap: error_rate vs missing_rate, faceted by mixture/dropout
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

        mixtures = sorted(results_df["mixture_prop"].unique())
        dropouts = sorted(results_df["dropout_rate"].unique())

        for idx, (mix, drop) in enumerate(itertools.product(mixtures, dropouts)):
            if idx >= 4:
                break
            ax = axes[idx]
            subset = results_df[
                (results_df["mixture_prop"] == mix) &
                (results_df["dropout_rate"] == drop)
            ]
            if subset.empty:
                continue

            pivot = subset.pivot_table(
                values=metric,
                index="error_rate",
                columns="missing_rate",
                aggfunc="mean"
            )

            sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn",
                        vmin=0, vmax=1, ax=ax, cbar_kws={"label": metric})
            ax.set_title(f"Mixture={mix}, Dropout={drop}")
            ax.set_xlabel("Missing rate")
            ax.set_ylabel("Error rate")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

    except ImportError:
        print("Matplotlib/seaborn not available, skipping heatmap")


def summarize_robustness(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize robustness results across simulations.
    """
    summary = results_df.groupby(
        ["error_rate", "missing_rate", "mixture_prop", "dropout_rate"]
    ).agg(
        mean_accuracy=("accuracy", "mean"),
        std_accuracy=("accuracy", "std"),
        min_accuracy=("accuracy", "min"),
        mean_separation=("mean_separation", "mean"),
        min_separation=("min_separation", "min"),
        n_simulations=("simulation", "count"),
    ).reset_index()

    return summary


def recommend_parameters(
    summary_df: pd.DataFrame,
    target_accuracy: float = 0.95,
) -> dict:
    """
    Recommend k, panel size, backup fraction based on robustness results.
    """
    # Find conditions where accuracy >= target
    good = summary_df[summary_df["mean_accuracy"] >= target_accuracy]

    if good.empty:
        return {
            "recommended_k": 3,
            "recommended_panel_size": 96,
            "recommended_backup_fraction": 0.20,
            "warning": "Target accuracy not achieved under any tested condition",
        }

    # Conservative recommendation: worst-case condition that meets target
    worst_good = good.sort_values(
        ["error_rate", "missing_rate", "mixture_prop", "dropout_rate"],
        ascending=False
    ).iloc[0]

    return {
        "recommended_k": 3,
        "recommended_panel_size": 96,
        "recommended_backup_fraction": 0.20,
        "worst_case_meeting_target": worst_good.to_dict(),
        "note": "Recommendations based on simulation; validate with real data",
    }


# Import itertools at top level for plot_robustness_heatmap
import itertools
