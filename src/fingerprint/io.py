"""
I/O utilities for genotype data: VCF, PLINK2, Parquet, metadata.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import allel
import numpy as np
import pandas as pd

# Optional imports - only needed for specific functions
try:
    import pysam
    _HAS_PYSAM = True
except ImportError:
    _HAS_PYSAM = False


@dataclass
class GenotypeData:
    """Container for genotype matrix and metadata."""
    gt: np.ndarray                    # Shape: (n_variants, n_samples), int8 (0/1/2, -1=missing)
    variant_ids: np.ndarray           # SNP IDs
    sample_ids: np.ndarray            # Sample/variety IDs
    chrom: np.ndarray                 # Chromosome per variant
    pos: np.ndarray                   # Position per variant
    ref: np.ndarray                   # Reference allele
    alt: np.ndarray                   # Alternate allele
    metadata: dict = field(default_factory=dict)

    @property
    def n_variants(self) -> int:
        return self.gt.shape[0]

    @property
    def n_samples(self) -> int:
        return self.gt.shape[1]

    def subset_variants(self, indices: np.ndarray) -> GenotypeData:
        """Return a new GenotypeData with only selected variants."""
        return GenotypeData(
            gt=self.gt[indices],
            variant_ids=self.variant_ids[indices],
            sample_ids=self.sample_ids,
            chrom=self.chrom[indices],
            pos=self.pos[indices],
            ref=self.ref[indices],
            alt=self.alt[indices],
            metadata=self.metadata.copy(),
        )

    def subset_samples(self, indices: np.ndarray) -> GenotypeData:
        """Return a new GenotypeData with only selected samples."""
        return GenotypeData(
            gt=self.gt[:, indices],
            variant_ids=self.variant_ids,
            sample_ids=self.sample_ids[indices],
            chrom=self.chrom,
            pos=self.pos,
            ref=self.ref,
            alt=self.alt,
            metadata=self.metadata.copy(),
        )


def read_vcf_to_genotype(
    vcf_path: str | Path,
    region: str | None = None,
    samples: list[str] | None = None,
    max_missing: float = 1.0,
) -> GenotypeData:
    """
    Read VCF/BCF to GenotypeData using scikit-allel.

    Parameters
    ----------
    vcf_path : Path to VCF/BCF file (bgzipped and indexed)
    region : Optional genomic region (e.g., "chr1:1000000-2000000")
    samples : Optional list of sample names to include
    max_missing : Maximum missing rate per variant (1.0 = no filter)

    Returns
    -------
    GenotypeData object
    """
    vcf_path = str(vcf_path)

    # Open with allel
    callset = allel.read_vcf(
        vcf_path,
        region=region,
        samples=samples,
        fields=[
            "samples",
            "variants/ID",
            "variants/CHROM",
            "variants/POS",
            "variants/REF",
            "variants/ALT",
            "calldata/GT",
        ],
    )

    gt = allel.GenotypeArray(callset["calldata/GT"])  # (n_var, n_samp, 2)
    n_called = (gt[:, :, 0] >= 0).sum(axis=1)
    miss_rate = 1 - n_called / gt.shape[1]
    keep = miss_rate <= max_missing

    # Convert to dosage (0, 1, 2, -1)
    dos = gt.to_n_alt(fill=-1).astype(np.int8)  # (n_var, n_samp)

    return GenotypeData(
        gt=dos[keep],
        variant_ids=callset["variants/ID"][keep],
        sample_ids=callset["samples"],
        chrom=callset["variants/CHROM"][keep],
        pos=callset["variants/POS"][keep],
        ref=callset["variants/REF"][keep],
        alt=np.array([a[0] for a in callset["variants/ALT"][keep]]),
        metadata={"source": vcf_path, "region": region, "max_missing": max_missing},
    )


def read_plink2(
    prefix: str | Path,
    variant_subset: np.ndarray | None = None,
    sample_subset: np.ndarray | None = None,
) -> GenotypeData:
    """
    Read PLINK2 binary files (.pgen, .pvar, .psam) to GenotypeData.

    Uses pysam for efficient reading.
    """
    prefix = Path(prefix)
    pgen_path = prefix.with_suffix(".pgen")
    pvar_path = prefix.with_suffix(".pvar")
    psam_path = prefix.with_suffix(".psam")

    # Read variant info
    pvar = pd.read_csv(pvar_path, sep="\t", comment="#")
    if variant_subset is not None:
        pvar = pvar.iloc[variant_subset].reset_index(drop=True)

    # Read sample info
    psam = pd.read_csv(psam_path, sep="\t", comment="#")
    if sample_subset is not None:
        psam = psam.iloc[sample_subset].reset_index(drop=True)

    # Read genotype data using pysam
    # Note: pysam doesn't directly support .pgen, use plink2 --export or allel
    # For now, use allel's read_plink (requires .bed/.bim/.fam) or convert
    # This is a placeholder - actual implementation depends on format
    raise NotImplementedError("Use plink2 --export vcf and read_vcf_to_genotype for now")


def write_genotype_parquet(data: GenotypeData, path: str | Path) -> None:
    """Write GenotypeData to Parquet for fast random access."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Create variant table
    var_df = pd.DataFrame({
        "variant_id": data.variant_ids,
        "chrom": data.chrom,
        "pos": data.pos,
        "ref": data.ref,
        "alt": data.alt,
    })

    # Create genotype matrix as DataFrame (samples as columns)
    gt_df = pd.DataFrame(
        data.gt.T,  # transpose to (samples x variants)
        index=data.sample_ids,
        columns=data.variant_ids,
    )

    # Write both
    var_df.to_parquet(path.with_suffix(".variants.parquet"), index=False)
    gt_df.to_parquet(path.with_suffix(".genotypes.parquet"))

    # Write metadata
    with open(path.with_suffix(".meta.json"), "w") as f:
        json.dump(data.metadata, f, indent=2)


def read_genotype_parquet(path: str | Path) -> GenotypeData:
    """Read GenotypeData from Parquet files."""
    path = Path(path)
    var_df = pd.read_parquet(path.with_suffix(".variants.parquet"))
    gt_df = pd.read_parquet(path.with_suffix(".genotypes.parquet"))

    with open(path.with_suffix(".meta.json")) as f:
        metadata = json.load(f)

    return GenotypeData(
        gt=gt_df.T.values.astype(np.int8),  # back to (variants x samples)
        variant_ids=var_df["variant_id"].values,
        sample_ids=gt_df.index.values,
        chrom=var_df["chrom"].values,
        pos=var_df["pos"].values,
        ref=var_df["ref"].values,
        alt=var_df["alt"].values,
        metadata=metadata,
    )


def load_reference_panel(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load a reference fingerprint panel for matching.

    Returns
    -------
    (genotypes, variant_ids, sample_ids)
    """
    data = read_genotype_parquet(path)
    return data.gt, data.variant_ids, data.sample_ids


def read_manifest(manifest_path: str | Path) -> pd.DataFrame:
    """Read data manifest TSV."""
    return pd.read_csv(manifest_path, sep="\t")


def write_manifest(df: pd.DataFrame, manifest_path: str | Path) -> None:
    """Write data manifest TSV."""
    df.to_csv(manifest_path, sep="\t", index=False)


def compute_file_sha256(filepath: str | Path) -> str:
    """Compute SHA-256 checksum of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def harmonise_sample_ids(
    sample_ids: np.ndarray,
    alias_table: pd.DataFrame,
    id_col: str = "alias",
    canonical_col: str = "canonical_id",
) -> tuple[np.ndarray, dict[str, str]]:
    """
    Harmonise sample IDs using an alias table.

    Parameters
    ----------
    sample_ids : Array of sample IDs to harmonise
    alias_table : DataFrame with columns [alias, canonical_id]
    id_col : Column name for aliases
    canonical_col : Column name for canonical IDs

    Returns
    -------
    (harmonised_ids, mapping_dict)
    """
    alias_map = dict(zip(alias_table[id_col], alias_table[canonical_col]))
    harmonised = np.array([alias_map.get(s, s) for s in sample_ids])
    return harmonised, alias_map


def build_alias_table(
    sample_lists: list[np.ndarray],
    similarity_threshold: float = 0.95,
) -> pd.DataFrame:
    """
    Build an alias table from multiple sample lists using fuzzy matching.

    This is a simplified version - in practice use a dedicated entity resolution tool.
    """
    from difflib import SequenceMatcher

    all_samples = set()
    for sl in sample_lists:
        all_samples.update(sl)

    all_samples = sorted(all_samples)
    canonical = {}
    used = set()

    for s in all_samples:
        if s in used:
            continue
        # Find similar names
        matches = [s]
        for t in all_samples:
            if t in used or t == s:
                continue
            if SequenceMatcher(None, s.lower(), t.lower()).ratio() >= similarity_threshold:
                matches.append(t)
        for m in matches:
            canonical[m] = s
            used.add(m)

    return pd.DataFrame({
        "alias": list(canonical.keys()),
        "canonical_id": list(canonical.values()),
    })
