"""
SNP scoring for informativeness and assay-friendliness.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

# Optional imports
try:
    import pysam
    _HAS_PYSAM = True
except ImportError:
    _HAS_PYSAM = False


@dataclass
class SNPScores:
    """Container for per-SNP scores."""
    variant_ids: np.ndarray
    chrom: np.ndarray
    pos: np.ndarray
    ref: np.ndarray
    alt: np.ndarray

    # Informativeness scores
    maf: np.ndarray
    pic: np.ndarray
    he: np.ndarray
    pi_locus: np.ndarray

    # Assay-friendliness scores
    flank_uniqueness: np.ndarray
    flank_gc_content: np.ndarray
    nearby_variants_count: np.ndarray
    repeat_overlap: np.ndarray

    # Composite score
    composite_score: np.ndarray

    # Metadata
    n_samples: int
    call_rates: np.ndarray
    het_rates: np.ndarray


def calculate_informativeness(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    min_called: int = 10,
) -> dict[str, np.ndarray]:
    """
    Calculate per-SNP informativeness metrics.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples) coded 0/1/2, -1=missing
    variant_ids : Variant identifiers
    min_called : Minimum called samples for valid calculation

    Returns
    -------
    Dictionary with maf, pic, he, pi_locus, call_rate, het_rate
    """
    n_variants, n_samples = gt.shape
    called = gt >= 0
    n_called = called.sum(axis=1)

    # Allele frequency (dosage of alt allele)
    alt_count = np.where(called, gt, 0).sum(axis=1)
    alt_freq = alt_count / (2 * np.maximum(n_called, 1))
    maf = np.minimum(alt_freq, 1 - alt_freq)

    # Heterozygosity rate
    het_rate = np.where(called, gt == 1, 0).sum(axis=1) / np.maximum(n_called, 1)

    # Call rate
    call_rate = n_called / n_samples

    # Expected heterozygosity / discrimination power (for inbreds)
    he = 2 * alt_freq * (1 - alt_freq)

    # Polymorphism Information Content
    # PIC = 2pq(1 - pq) where p = alt_freq, q = 1 - alt_freq
    pq = alt_freq * (1 - alt_freq)
    pic = 2 * pq * (1 - pq)

    # Probability of identity at this locus
    pi_locus = alt_freq**2 + (1 - alt_freq)**2

    # Mask low-call variants
    valid = n_called >= min_called
    maf[~valid] = 0
    pic[~valid] = 0
    he[~valid] = 0
    pi_locus[~valid] = 1.0

    return {
        "maf": maf,
        "pic": pic,
        "he": he,
        "pi_locus": pi_locus,
        "call_rate": call_rate,
        "het_rate": het_rate,
        "alt_freq": alt_freq,
        "n_called": n_called,
    }


def extract_flank_sequence(
    fasta_path: str,
    chrom: str,
    pos: int,
    flank_length: int = 60,
) -> str | None:
    """
    Extract flanking sequence around a SNP using samtools faidx.

    Parameters
    ----------
    fasta_path : Path to indexed FASTA file
    chrom : Chromosome name
    pos : 1-based position
    flank_length : Number of bases on each side

    Returns
    -------
    Flank sequence (uppercase) or None if error
    """
    import subprocess

    start = max(1, pos - flank_length)
    end = pos + flank_length
    region = f"{chrom}:{start}-{end}"

    try:
        result = subprocess.run(
            ["samtools", "faidx", fasta_path, region],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        )
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 2:
            seq = "".join(lines[1:]).upper()
            # Verify SNP position
            snp_idx = flank_length  # 0-based in extracted sequence
            if len(seq) == 2 * flank_length + 1:
                return seq
    except (subprocess.SubprocessError, FileNotFoundError, IndexError):
        pass
    return None


def check_flank_uniqueness(
    fasta_path: str,
    chrom: np.ndarray,
    pos: np.ndarray,
    flank_length: int = 60,
    min_identity: float = 0.95,
) -> np.ndarray:
    """
    Check flank uniqueness by BLASTing against reference genome.

    Returns a score 0-1 where 1 = unique, 0 = repetitive.
    """
    # This is a simplified placeholder
    # Real implementation would use BLAST or minimap2
    # For now, return ones (assume unique)
    return np.ones(len(chrom))


def calculate_flank_gc_content(
    fasta_path: str,
    chrom: np.ndarray,
    pos: np.ndarray,
    flank_length: int = 60,
) -> np.ndarray:
    """
    Calculate GC content of flanking sequences.
    """
    gc_scores = np.zeros(len(chrom))

    for i, (c, p) in enumerate(zip(chrom, pos)):
        seq = extract_flank_sequence(fasta_path, c, p, flank_length)
        if seq:
            gc = (seq.count("G") + seq.count("C")) / len(seq)
            gc_scores[i] = gc
        else:
            gc_scores[i] = 0.5  # Neutral if can't extract

    return gc_scores


def count_nearby_variants(
    vcf_path: str,
    chrom: np.ndarray,
    pos: np.ndarray,
    window_bp: int = 50,
) -> np.ndarray:
    """
    Count other variants within a window of each SNP.
    """
    counts = np.zeros(len(chrom), dtype=int)

    # Use tabix to query region
    if not _HAS_PYSAM:
        return counts

    try:
        tbx = pysam.TabixFile(vcf_path)
        for i, (c, p) in enumerate(zip(chrom, pos)):
            region = f"{c}:{max(1, p - window_bp)}-{p + window_bp}"
            try:
                variants = list(tbx.fetch(region=region))
                counts[i] = len(variants) - 1  # Exclude self
            except ValueError:
                counts[i] = 0
        tbx.close()
    except (OSError, ValueError):
        pass

    return counts


def check_repeat_overlap(
    chrom: np.ndarray,
    pos: np.ndarray,
    repeat_bed: str,
) -> np.ndarray:
    """
    Check if SNPs overlap repeat regions.
    Returns 1 if no overlap, 0 if in repeat.
    """
    # Placeholder - would use pybedtools
    return np.ones(len(chrom))


def score_snps(
    gt: np.ndarray,
    variant_ids: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    ref: np.ndarray,
    alt: np.ndarray,
    fasta_path: str,
    vcf_path: str | None = None,
    repeat_bed: str | None = None,
    config: dict | None = None,
) -> SNPScores:
    """
    Score all candidate SNPs for informativeness and assay-friendliness.

    Parameters
    ----------
    gt : Genotype matrix (n_variants, n_samples)
    variant_ids : SNP IDs
    chrom, pos, ref, alt : Variant coordinates and alleles
    fasta_path : Path to reference genome FASTA (indexed)
    vcf_path : Path to VCF for nearby variant counting
    repeat_bed : Path to repeat mask BED file
    config : Configuration dictionary with weights

    Returns
    -------
    SNPScores object
    """
    if config is None:
        config = {
            "pic_weight": 0.4,
            "he_weight": 0.3,
            "flank_uniqueness_weight": 0.2,
            "spacing_weight": 0.1,
            "min_maf": 0.20,
            "max_het": 0.05,
            "min_call_rate": 0.95,
            "flank_length": 60,
        }

    n_variants = len(variant_ids)

    # 1. Informativeness scores
    info = calculate_informativeness(gt, variant_ids)

    # 2. Assay-friendliness scores
    flank_unique = check_flank_uniqueness(
        fasta_path, chrom, pos, config["flank_length"]
    )
    flank_gc = calculate_flank_gc_content(
        fasta_path, chrom, pos, config["flank_length"]
    )

    nearby_vars = np.zeros(n_variants)
    if vcf_path:
        nearby_vars = count_nearby_variants(vcf_path, chrom, pos)

    repeat_ovlp = np.ones(n_variants)
    if repeat_bed:
        repeat_ovlp = check_repeat_overlap(chrom, pos, repeat_bed)

    # 3. Normalize scores to 0-1 range for composite
    # PIC and HE are already 0-1
    # Flank uniqueness: already 0-1
    # GC content: penalize extreme values (target 40-60%)
    gc_score = 1 - 2 * np.abs(flank_gc - 0.5)  # 1 at 0.5, 0 at 0 or 1
    gc_score = np.clip(gc_score, 0, 1)

    # Nearby variants: fewer is better
    nearby_score = 1 / (1 + nearby_vars / 5)  # Decay with count

    # Repeat overlap: 1 if no overlap, 0 if overlap
    repeat_score = repeat_ovlp

    # 4. Composite score (weighted sum of normalized components)
    composite = (
        config["pic_weight"] * info["pic"] +
        config["he_weight"] * info["he"] +
        config["flank_uniqueness_weight"] * flank_unique +
        config["spacing_weight"] * (gc_score + nearby_score + repeat_score) / 3
    )

    # Zero out composite for SNPs failing hard filters
    hard_filter = (
        (info["maf"] >= config["min_maf"]) &
        (info["het_rate"] <= config["max_het"]) &
        (info["call_rate"] >= config["min_call_rate"])
    )
    composite[~hard_filter] = 0

    return SNPScores(
        variant_ids=variant_ids,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        maf=info["maf"],
        pic=info["pic"],
        he=info["he"],
        pi_locus=info["pi_locus"],
        flank_uniqueness=flank_unique,
        flank_gc_content=flank_gc,
        nearby_variants_count=nearby_vars,
        repeat_overlap=repeat_ovlp,
        composite_score=composite,
        n_samples=gt.shape[1],
        call_rates=info["call_rate"],
        het_rates=info["het_rate"],
    )


def write_scores_table(scores: SNPScores, output_path: str) -> None:
    """Write SNP scores to TSV."""
    df = pd.DataFrame({
        "variant_id": scores.variant_ids,
        "chrom": scores.chrom,
        "pos": scores.pos,
        "ref": scores.ref,
        "alt": scores.alt,
        "maf": scores.maf,
        "pic": scores.pic,
        "he": scores.he,
        "pi_locus": scores.pi_locus,
        "flank_uniqueness": scores.flank_uniqueness,
        "flank_gc_content": scores.flank_gc_content,
        "nearby_variants_count": scores.nearby_variants_count,
        "repeat_overlap": scores.repeat_overlap,
        "composite_score": scores.composite_score,
        "call_rate": scores.call_rates,
        "het_rate": scores.het_rates,
    })
    df.to_csv(output_path, sep="\t", index=False)


def read_scores_table(input_path: str) -> SNPScores:
    """Read SNP scores from TSV."""
    df = pd.read_csv(input_path, sep="\t")
    return SNPScores(
        variant_ids=df["variant_id"].values,
        chrom=df["chrom"].values,
        pos=df["pos"].values,
        ref=df["ref"].values,
        alt=df["alt"].values,
        maf=df["maf"].values,
        pic=df["pic"].values,
        he=df["he"].values,
        pi_locus=df["pi_locus"].values,
        flank_uniqueness=df["flank_uniqueness"].values,
        flank_gc_content=df["flank_gc_content"].values,
        nearby_variants_count=df["nearby_variants_count"].values,
        repeat_overlap=df["repeat_overlap"].values,
        composite_score=df["composite_score"].values,
        n_samples=0,  # Not stored in table
        call_rates=df["call_rate"].values,
        het_rates=df["het_rate"].values,
    )


def rank_candidates(scores: SNPScores, top_n: int | None = None) -> pd.DataFrame:
    """Rank candidates by composite score."""
    df = pd.DataFrame({
        "variant_id": scores.variant_ids,
        "chrom": scores.chrom,
        "pos": scores.pos,
        "ref": scores.ref,
        "alt": scores.alt,
        "maf": scores.maf,
        "pic": scores.pic,
        "he": scores.he,
        "composite_score": scores.composite_score,
        "call_rate": scores.call_rates,
        "het_rate": scores.het_rates,
    })
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    if top_n:
        df = df.head(top_n)
    return df
