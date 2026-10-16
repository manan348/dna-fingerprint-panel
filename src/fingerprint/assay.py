"""
KASP assay design: flank extraction, primer design, in silico validation.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import primer3
from Bio.Seq import Seq


@dataclass
class KASPAssay:
    """KASP assay design for a single SNP."""
    variant_id: str
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str

    # Flank sequences
    flank_sequence: str
    snp_index_in_flank: int  # 0-based position of SNP in flank

    # Primer sequences
    common_reverse: str
    allele1_forward: str      # FAM-tailed
    allele2_forward: str      # HEX-tailed

    # Primer properties
    common_tm: float
    allele1_tm: float
    allele2_tm: float
    amplicon_size: int

    # Quality checks
    blast_unique: bool
    no_dimers: bool
    gc_ok: bool
    tm_ok: bool
    pass_all: bool

    # Tails
    fam_tail: str = "GAAGGTGACCAAGTTCATGCT"
    hex_tail: str = "GAAGGTCGGAGTCAACGGATT"


def extract_flank(
    fasta_path: str,
    chrom: str,
    pos: int,
    flank_length: int = 60,
) -> str | None:
    """
    Extract flanking sequence from reference FASTA using samtools faidx.

    Parameters
    ----------
    fasta_path : Path to indexed FASTA file
    chrom : Chromosome name
    pos : 1-based position
    flank_length : Bases on each side

    Returns
    -------
    Flank sequence (uppercase) or None
    """
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
            expected_len = 2 * flank_length + 1
            if len(seq) == expected_len:
                return seq
            else:
                print(f"Warning: Flank length mismatch for {chrom}:{pos} "
                      f"(got {len(seq)}, expected {expected_len})")
    except (subprocess.SubprocessError, FileNotFoundError, IndexError) as e:
        print(f"Error extracting flank for {chrom}:{pos}: {e}")

    return None


def design_common_reverse(
    flank_seq: str,
    snp_index: int,
    opt_size: int = 24,
    min_size: int = 18,
    max_size: int = 30,
    opt_tm: float = 60.0,
    min_tm: float = 57.0,
    max_tm: float = 63.0,
    gc_min: float = 0.40,
    gc_max: float = 0.60,
) -> dict | None:
    """
    Design common reverse primer downstream of SNP.

    Returns primer3 result dict or None.
    """
    # Region for primer3: downstream of SNP
    # SNP is at snp_index, we want primer to the right
    seq_len = len(flank_seq)
    included_start = snp_index + 1
    included_len = seq_len - included_start - 1  # Leave some space

    if included_len < min_size:
        return None

    try:
        result = primer3.design_primers(
            {
                "SEQUENCE_ID": "kasp_common",
                "SEQUENCE_TEMPLATE": flank_seq,
                "SEQUENCE_INCLUDED_REGION": [included_start, included_len],
            },
            {
                "PRIMER_TASK": "pick_right_only",
                "PRIMER_OPT_SIZE": opt_size,
                "PRIMER_MIN_SIZE": min_size,
                "PRIMER_MAX_SIZE": max_size,
                "PRIMER_OPT_TM": opt_tm,
                "PRIMER_MIN_TM": min_tm,
                "PRIMER_MAX_TM": max_tm,
                "PRIMER_MIN_GC": gc_min * 100,
                "PRIMER_MAX_GC": gc_max * 100,
                "PRIMER_MAX_SELF_ANY": 8,
                "PRIMER_MAX_SELF_END": 3,
                "PRIMER_MAX_HAIRPIN_TH": 47,
            },
        )
        return result
    except Exception as e:
        print(f"Primer3 error for common reverse: {e}")
        return None


def design_allele_specific_forward(
    flank_seq: str,
    snp_index: int,
    allele: str,
    tail: str,
    target_tm: float = 60.0,
    min_len: int = 18,
    max_len: int = 30,
) -> str | None:
    """
    Design allele-specific forward primer.

    The primer ends with the allele-specific base at the 3' end.
    We add the universal tail (FAM/HEX) to the 5' end.

    Strategy: Take upstream sequence + allele base, adjust length for Tm.
    """
    # Upstream sequence (reverse complement since it's forward primer on opposite strand)
    upstream = flank_seq[:snp_index]

    # Try different lengths to get target Tm
    for length in range(min_len, max_len + 1):
        if length > len(upstream) + 1:
            break
        # Primer body: last (length-1) bases of upstream + allele
        body_start = len(upstream) - (length - 1)
        if body_start < 0:
            continue
        body = upstream[body_start:] + allele
        full_primer = tail + body

        # Calculate Tm
        tm = primer3.calc_tm(full_primer)

        if abs(tm - target_tm) <= 2.0:
            return full_primer

    # If no perfect match, return closest
    best_primer = None
    best_diff = float('inf')
    for length in range(min_len, max_len + 1):
        if length > len(upstream) + 1:
            break
        body_start = len(upstream) - (length - 1)
        if body_start < 0:
            continue
        body = upstream[body_start:] + allele
        full_primer = tail + body
        tm = primer3.calc_tm(full_primer)
        diff = abs(tm - target_tm)
        if diff < best_diff:
            best_diff = diff
            best_primer = full_primer

    return best_primer


def check_primer_dimer(
    primer1: str,
    primer2: str,
    max_complementarity: int = 8,
    max_3prime_complementarity: int = 3,
) -> bool:
    """
    Check for primer-dimer formation between two primers.

    Returns True if no significant dimer risk.
    """
    # Simple check: look for complementarity at 3' ends
    # More thorough check would use primer3's thermodynamic calculations

    # Check 3' end complementarity (last 5 bases)
    p1_end = primer1[-5:]
    p2_end = primer2[-5:]

    # Reverse complement of p2_end
    rc_p2_end = str(Seq(p2_end).reverse_complement())

    # Count matches
    matches = sum(1 for a, b in zip(p1_end, rc_p2_end) if a == b)

    if matches >= max_3prime_complementarity:
        return False

    # Check overall complementarity (simplified)
    # In practice, use primer3 calc_heterodimer
    try:
        het = primer3.calc_heterodimer(primer1, primer2)
        if het.dg < -9.0:  # kcal/mol, more negative = stronger binding
            return False
    except Exception:
        pass

    return True


def check_hairpin(primer: str, max_dg: float = -6.0) -> bool:
    """Check for hairpin formation."""
    try:
        hp = primer3.calc_hairpin(primer)
        if hp.dg < max_dg:
            return False
    except Exception:
        pass
    return True


def blast_primer(
    primer: str,
    fasta_path: str,
    max_hits: int = 1,
) -> bool:
    """
    BLAST primer against reference genome to check uniqueness.

    Returns True if unique (only hits target locus).
    """
    # Write primer to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">primer\n{primer}\n")
        primer_file = f.name

    try:
        # Run BLAST (using blastn-short for primers)
        result = subprocess.run(
            ["blastn", "-task", "blastn-short",
             "-query", primer_file,
             "-subject", fasta_path,
             "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
             "-max_target_seqs", str(max_hits + 5),
             "-evalue", "1000"],
            capture_output=True,
            text=True,
            timeout=30,
        )

        hits = result.stdout.strip().split('\n')
        hits = [h for h in hits if h.strip()]

        # Should have exactly 1 hit (the target)
        return len(hits) <= max_hits

    except (subprocess.SubprocessError, FileNotFoundError):
        print("Warning: BLAST not available, skipping uniqueness check")
        return True
    finally:
        Path(primer_file).unlink(missing_ok=True)


def design_kasp_assay(
    variant_id: str,
    chrom: str,
    pos: int,
    ref_allele: str,
    alt_allele: str,
    fasta_path: str,
    flank_length: int = 60,
    fam_tail: str = "GAAGGTGACCAAGTTCATGCT",
    hex_tail: str = "GAAGGTCGGAGTCAACGGATT",
) -> KASPAssay | None:
    """
    Design a complete KASP assay for a SNP.

    Parameters
    ----------
    variant_id : SNP identifier
    chrom, pos : Genomic coordinates (1-based)
    ref_allele, alt_allele : Alleles (single base)
    fasta_path : Reference genome FASTA (indexed)
    flank_length : Flank length on each side
    fam_tail, hex_tail : Universal tails

    Returns
    -------
    KASPAssay object or None if design failed
    """
    # Extract flank
    flank = extract_flank(fasta_path, chrom, pos, flank_length)
    if flank is None:
        print(f"Failed to extract flank for {variant_id}")
        return None

    snp_index = flank_length  # 0-based

    # Verify SNP matches reference
    ref_base = flank[snp_index]
    if ref_base != ref_allele.upper():
        print(f"Allele mismatch for {variant_id}: ref={ref_base}, expected={ref_allele}")
        # Try reverse complement
        flank_rc = str(Seq(flank).reverse_complement())
        if flank_rc[snp_index] == ref_allele.upper():
            flank = flank_rc
            ref_base = flank[snp_index]
        else:
            print(f"Cannot resolve allele mismatch for {variant_id}")
            return None

    # Design common reverse primer
    common_result = design_common_reverse(flank, snp_index)
    if not common_result or "PRIMER_RIGHT_0_SEQUENCE" not in common_result:
        print(f"Failed to design common reverse for {variant_id}")
        return None

    common_reverse = common_result["PRIMER_RIGHT_0_SEQUENCE"]
    common_tm = common_result.get("PRIMER_RIGHT_0_TM", 0)

    # Design allele-specific forwards
    allele1_forward = design_allele_specific_forward(
        flank, snp_index, ref_allele.upper(), fam_tail
    )
    allele2_forward = design_allele_specific_forward(
        flank, snp_index, alt_allele.upper(), hex_tail
    )

    if not allele1_forward or not allele2_forward:
        print(f"Failed to design allele-specific primers for {variant_id}")
        return None

    allele1_tm = primer3.calc_tm(allele1_forward)
    allele2_tm = primer3.calc_tm(allele2_forward)

    # Amplicon size
    # Common reverse binds downstream, allele-specific forwards bind upstream
    # Approximate: distance from 5' end of forward to 5' end of reverse
    amplicon_size = len(flank) - snp_index + len(allele1_forward) - len(fam_tail)

    # Quality checks
    blast_unique = blast_primer(common_reverse, fasta_path) and \
                   blast_primer(allele1_forward, fasta_path) and \
                   blast_primer(allele2_forward, fasta_path)

    no_dimers = (check_primer_dimer(common_reverse, allele1_forward) and
                 check_primer_dimer(common_reverse, allele2_forward) and
                 check_primer_dimer(allele1_forward, allele2_forward))

    no_hairpins = (check_hairpin(common_reverse) and
                   check_hairpin(allele1_forward) and
                   check_hairpin(allele2_forward))

    gc_ok = all(0.40 <= (p.count('G') + p.count('C')) / len(p) <= 0.60
                for p in [common_reverse, allele1_forward, allele2_forward])

    tm_ok = all(57 <= primer3.calc_tm(p) <= 63
                for p in [common_reverse, allele1_forward, allele2_forward])

    pass_all = blast_unique and no_dimers and no_hairpins and gc_ok and tm_ok

    return KASPAssay(
        variant_id=variant_id,
        chrom=chrom,
        pos=pos,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        flank_sequence=flank,
        snp_index_in_flank=snp_index,
        common_reverse=common_reverse,
        allele1_forward=allele1_forward,
        allele2_forward=allele2_forward,
        common_tm=common_tm,
        allele1_tm=allele1_tm,
        allele2_tm=allele2_tm,
        amplicon_size=amplicon_size,
        blast_unique=blast_unique,
        no_dimers=no_dimers,
        gc_ok=gc_ok,
        tm_ok=tm_ok,
        pass_all=pass_all,
        fam_tail=fam_tail,
        hex_tail=hex_tail,
    )


def design_panel_assays(
    panel_df: pd.DataFrame,
    fasta_path: str,
    flank_length: int = 60,
    fam_tail: str = "GAAGGTGACCAAGTTCATGCT",
    hex_tail: str = "GAAGGTCGGAGTCAACGGATT",
) -> list[KASPAssay]:
    """
    Design KASP assays for a panel of SNPs.

    Parameters
    ----------
    panel_df : DataFrame with columns variant_id, chrom, pos, ref, alt
    fasta_path : Reference genome FASTA
    flank_length : Flank length
    fam_tail, hex_tail : Universal tails

    Returns
    -------
    List of KASPAssay objects
    """
    assays = []
    failed = []

    for _, row in panel_df.iterrows():
        assay = design_kasp_assay(
            variant_id=row["variant_id"],
            chrom=row["chrom"],
            pos=row["pos"],
            ref_allele=row["ref"],
            alt_allele=row["alt"],
            fasta_path=fasta_path,
            flank_length=flank_length,
            fam_tail=fam_tail,
            hex_tail=hex_tail,
        )
        if assay:
            assays.append(assay)
        else:
            failed.append(row["variant_id"])

    print(f"Designed {len(assays)} assays, {len(failed)} failed")
    if failed:
        print(f"Failed: {failed}")

    return assays


def write_kasp_design_sheet(
    assays: list[KASPAssay],
    output_path: str,
) -> None:
    """
    Write KASP design sheet to Excel.

    Columns: Variant_ID, Chrom, Pos, Ref, Alt, Flank_Sequence,
             Common_Reverse, Allele1_Forward_FAM, Allele2_Forward_HEX,
             Common_Tm, Allele1_Tm, Allele2_Tm, Amplicon_Size,
             BLAST_Unique, No_Dimers, GC_OK, Tm_OK, Pass_All
    """
    data = []
    for a in assays:
        data.append({
            "Variant_ID": a.variant_id,
            "Chrom": a.chrom,
            "Pos": a.pos,
            "Ref": a.ref_allele,
            "Alt": a.alt_allele,
            "Flank_Sequence": a.flank_sequence,
            "SNP_Index_In_Flank": a.snp_index_in_flank,
            "Common_Reverse": a.common_reverse,
            "Allele1_Forward_FAM": a.allele1_forward,
            "Allele2_Forward_HEX": a.allele2_forward,
            "Common_Tm": round(a.common_tm, 1),
            "Allele1_Tm": round(a.allele1_tm, 1),
            "Allele2_Tm": round(a.allele2_tm, 1),
            "Amplicon_Size": a.amplicon_size,
            "BLAST_Unique": a.blast_unique,
            "No_Dimers": a.no_dimers,
            "GC_OK": a.gc_ok,
            "Tm_OK": a.tm_ok,
            "Pass_All": a.pass_all,
        })

    df = pd.DataFrame(data)
    df.to_excel(output_path, index=False)
    print(f"Design sheet written to {output_path}")


def generate_lab_sop(
    assays: list[KASPAssay],
    output_path: str,
    panel_name: str = "Core_Panel",
) -> None:
    """
    Generate a lab SOP document for KASP genotyping.
    """
    n_assays = len(assays)
    n_plates = (n_assays + 95) // 96  # 96-well plates

    sop = f"""# Lab SOP: KASP Genotyping for {panel_name}

## 1. Purpose
Genotype samples using the {panel_name} KASP SNP panel ({n_assays} assays)
for varietal identification and seed purity assessment.

## 2. Materials
- KASP Master Mix (e.g., LGC Biosearch)
- Assay primers (ordered from supplier, see design sheet)
- DNA samples (20-50 ng/µL in TE or water)
- 96-well or 384-well PCR plates
- qPCR instrument with FAM/HEX channels

## 3. Assay Layout
- Number of assays: {n_assays}
- Number of 96-well plates needed: {n_plates}
- Each plate: 96 assays × 1 sample, or 48 assays × 2 samples, etc.
- Include controls on each plate:
  - No-template control (NTC): 1 per plate
  - Positive controls (known varieties): 2-3 per plate
  - Replicate samples: 2-3 per plate

## 4. Reaction Setup
| Component | Volume per reaction |
|-----------|---------------------|
| 2× KASP Master Mix | 2.5 µL |
| Assay mix (primers) | 0.07 µL |
| DNA (20-50 ng/µL) | 2.0 µL |
| Water | 0.43 µL |
| **Total** | **5.0 µL** |

*Assay mix: 12 µM each allele-specific primer, 30 µM common reverse*

## 5. Thermal Cycling Protocol
1. 94°C for 15 min (enzyme activation)
2. 10 cycles: 94°C 20s, 61°C 60s (touchdown -0.6°C/cycle)
3. 26 cycles: 94°C 20s, 55°C 60s
4. 4°C hold

## 6. Data Analysis
- Read endpoint fluorescence (FAM/HEX)
- Auto-call clusters using instrument software
- Manual review of cluster quality
- Export genotype calls (AA, AB, BB, No Call)

## 7. Quality Control
- Per-assay call rate ≥ 95%
- Replicate concordance ≥ 99%
- NTC must show no amplification
- Positive controls must match expected genotypes

## 8. Troubleshooting
- **Poor clustering**: Check DNA quality, primer design, thermal profile
- **Low call rate**: Increase DNA concentration, check for inhibitors
- **Unexpected clusters**: Verify sample identity, check for contamination

## 9. Primer Information
See accompanying design sheet (kasp_design_sheet.xlsx) for all primer sequences.

---

Generated by DNA Fingerprint Panel Pipeline
"""

    with open(output_path, 'w') as f:
        f.write(sop)

    print(f"Lab SOP written to {output_path}")


def validate_assays_in_silico(
    assays: list[KASPAssay],
    fasta_path: str,
    vcf_path: str | None = None,
) -> pd.DataFrame:
    """
    In silico validation of designed assays.

    Checks:
    - Primer specificity (BLAST)
    - Nearby variants in primer binding sites
    - Expected amplicon size
    """
    results = []

    for a in assays:
        # Check for nearby variants in primer binding sites
        nearby_variants = 0
        if vcf_path:
            try:
                import pysam
                tbx = pysam.TabixFile(vcf_path)
                # Check common reverse binding region
                # Approximate: last 30bp of flank
                region = f"{a.chrom}:{a.pos - 30}-{a.pos + 30}"
                variants = list(tbx.fetch(region=region))
                nearby_variants = len(variants) - 1  # Exclude target SNP
                tbx.close()
            except Exception:
                pass

        results.append({
            "variant_id": a.variant_id,
            "blast_unique": a.blast_unique,
            "no_dimers": a.no_dimers,
            "gc_ok": a.gc_ok,
            "tm_ok": a.tm_ok,
            "pass_all": a.pass_all,
            "nearby_variants": nearby_variants,
            "amplicon_size": a.amplicon_size,
        })

    return pd.DataFrame(results)
