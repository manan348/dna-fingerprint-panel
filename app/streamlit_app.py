"""
Streamlit web interface for the DNA Fingerprint Panel toolkit.
"""

import json
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from fingerprint.io import load_reference_panel
from fingerprint.matcher import (
    Thresholds,
    classify_query,
    verify_hybrid,
    seed_lot_purity_test,
)
from fingerprint.assay import design_panel_assays, write_kasp_design_sheet


st.set_page_config(
    page_title="DNA Fingerprint Panel - Varietal Identification",
    page_icon=None,
    layout="wide",
)


@st.cache_resource
def load_reference(ref_path: str):
    """Load reference panel from Parquet."""
    return load_reference_panel(ref_path)


@st.cache_resource
def load_thresholds(thresh_path: str):
    """Load calibrated thresholds from JSON."""
    with open(thresh_path) as f:
        d = json.load(f)
    return Thresholds(**d)


def load_demo_data():
    """Load built-in demo data so users can immediately try the app."""
    demo_dir = Path("data/demo")
    ref_path = str(demo_dir / "reference_panel")
    thresh_path = str(demo_dir / "calibrated_thresholds.json")

    try:
        ref_gt, ref_var_ids, ref_sample_ids = load_reference(ref_path)
        st.session_state["ref_loaded"] = True
        st.session_state["ref_gt"] = ref_gt
        st.session_state["ref_var_ids"] = ref_var_ids
        st.session_state["ref_sample_ids"] = ref_sample_ids

        thresholds = load_thresholds(thresh_path)
        st.session_state["thresholds"] = thresholds
        st.session_state["ref_path"] = ref_path
        st.session_state["thresh_path"] = thresh_path

        return True, ref_path, thresh_path
    except Exception as e:
        return False, str(e), None


def main():
    st.title("DNA Fingerprint Panel - Varietal Identification")
    st.markdown(
        "Match unknown samples to reference varieties, verify hybrids, "
        "and design KASP assays for seed purity testing."
    )

    # Sidebar configuration
    with st.sidebar:
        st.header("Configuration")
        ref_path = st.text_input(
            "Reference Database",
            value="results/reference_panel.parquet",
            help="Path to Parquet reference panel",
        )
        thresh_path = st.text_input(
            "Thresholds JSON",
            value="results/tables/calibrated_thresholds.json",
            help="Path to calibrated thresholds",
        )

        st.divider()

        if st.button("Load Demo Data", type="primary", use_container_width=True):
            success, msg, _ = load_demo_data()
            if success:
                st.success("Demo data loaded successfully!")
                st.info(
                    "You can now try the Match Sample, Verify Hybrid, "
                    "and Seed Purity tabs with built-in data."
                )
            else:
                st.error(f"Failed to load demo data: {msg}")

        st.divider()

        if st.button("Load Reference", type="secondary"):
            try:
                ref_gt, ref_var_ids, ref_sample_ids = load_reference(ref_path)
                st.session_state["ref_loaded"] = True
                st.session_state["ref_gt"] = ref_gt
                st.session_state["ref_var_ids"] = ref_var_ids
                st.session_state["ref_sample_ids"] = ref_sample_ids
                st.success(
                    f"Loaded {ref_gt.shape[1]} varieties x "
                    f"{ref_gt.shape[0]} SNPs"
                )
            except Exception as e:
                st.error(f"Failed to load reference: {e}")

        if st.button("Load Thresholds"):
            try:
                thresholds = load_thresholds(thresh_path)
                st.session_state["thresholds"] = thresholds
                st.success(
                    f"Match threshold: {thresholds.match_threshold:.4f}"
                )
            except Exception as e:
                st.error(f"Failed to load thresholds: {e}")

        # Show demo file hints
        st.caption(
            "Demo files are in data/demo/. "
            "Generate them with: python scripts/generate_demo_data.py"
        )

    # Main tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Match Sample",
        "Verify Hybrid",
        "Seed Purity",
        "Design Assays",
        "Reference Panel",
    ])

    # Tab 1: Match Sample
    with tab1:
        st.header("Match Unknown Sample")
        st.markdown(
            "Upload a genotype file (CSV with columns: `snp_id`, `genotype`)"
        )

        col1, col2 = st.columns([3, 1])
        with col1:
            uploaded = st.file_uploader(
                "Query Genotype CSV", type=["csv"]
            )
        with col2:
            st.write("")
            st.write("")
            use_demo = st.button("Use Demo Query", key="demo_query")

        if use_demo:
            demo_query = Path("data/demo/query_sample.csv")
            if demo_query.exists():
                uploaded = demo_query
                st.info("Using demo query sample (matches variety V01)")
            else:
                st.error("Demo query file not found. Generate demo data first.")

        if uploaded:
            query_df = pd.read_csv(uploaded)
            st.write(f"Query has {len(query_df)} SNP calls")

            # Align to reference
            if "ref_loaded" not in st.session_state:
                st.error(
                    "Reference panel not loaded. "
                    "Please load demo data or your own reference panel first."
                )
            else:
                ref_var_ids = st.session_state["ref_var_ids"]
                query_gt = np.full(len(ref_var_ids), -1, dtype=np.int8)
                for _, row in query_df.iterrows():
                    idx = np.where(ref_var_ids == row["snp_id"])[0]
                    if len(idx) > 0:
                        query_gt[idx[0]] = row["genotype"]

                n_called = int((query_gt >= 0).sum())
                n_total = len(ref_var_ids)
                st.metric(
                    "Loci Called",
                    f"{n_called} / {n_total} ({n_called / n_total:.1%})",
                )

                if "thresholds" not in st.session_state:
                    st.error(
                        "Thresholds not loaded. "
                        "Please load demo data or your own thresholds first."
                    )
                else:
                    thresholds = st.session_state["thresholds"]
                    ref_gt = st.session_state["ref_gt"]
                    ref_sample_ids = st.session_state["ref_sample_ids"]

                    result = classify_query(
                        query_gt, ref_gt, ref_sample_ids, thresholds
                    )

                    # Display result
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Decision", result.decision)
                    with col2:
                        st.metric("Best Match", result.best_match)
                    with col3:
                        st.metric("Distance", f"{result.best_distance:.4f}")

                    st.metric("Loci Used", result.n_loci_used)
                    st.metric(
                        "Probability of Identity",
                        f"{result.probability_identity:.2e}",
                    )

                    if result.warnings:
                        for w in result.warnings:
                            st.warning(w)

                    # Top hits table
                    st.subheader("Top Matches")
                    hits_df = pd.DataFrame(
                        result.top_hits, columns=["Variety", "Distance"]
                    )
                    st.dataframe(hits_df, use_container_width=True)

    # Tab 2: Verify Hybrid
    with tab2:
        st.header("Verify F1 Hybrid")
        st.markdown(
            "Upload genotypes for the hybrid and both parents"
        )

        col1, col2, col3 = st.columns(3)
        with col1:
            h_file = st.file_uploader(
                "Hybrid Genotype", type=["csv"], key="hybrid"
            )
        with col2:
            p1_file = st.file_uploader(
                "Parent 1 Genotype", type=["csv"], key="parent1"
            )
        with col3:
            p2_file = st.file_uploader(
                "Parent 2 Genotype", type=["csv"], key="parent2"
            )

        st.write("")
        col_demo1, col_demo2, col_demo3 = st.columns(3)
        with col_demo1:
            use_demo_hybrid = st.button("Use Demo Files", key="demo_hybrid")

        if use_demo_hybrid:
            demo_dir = Path("data/demo")
            if all(
                (demo_dir / f).exists()
                for f in [
                    "hybrid_genotype.csv",
                    "parent1_genotype.csv",
                    "parent2_genotype.csv",
                ]
            ):
                h_file = demo_dir / "hybrid_genotype.csv"
                p1_file = demo_dir / "parent1_genotype.csv"
                p2_file = demo_dir / "parent2_genotype.csv"
                st.info("Using demo hybrid verification files")
            else:
                st.error("Demo files not found. Generate demo data first.")

        if h_file and p1_file and p2_file:
            h_df = pd.read_csv(h_file)
            p1_df = pd.read_csv(p1_file)
            p2_df = pd.read_csv(p2_file)

            # Merge on snp_id
            merged = h_df.merge(
                p1_df, on="snp_id", suffixes=("", "_p1")
            )
            merged = merged.merge(
                p2_df, on="snp_id", suffixes=("", "_p2")
            )

            h_gt = merged["genotype"].values
            p1_gt = merged["genotype_p1"].values
            p2_gt = merged["genotype_p2"].values

            result = verify_hybrid(h_gt, p1_gt, p2_gt)

            st.subheader("Verification Result")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Concordance", f"{result['concordance']:.2%}")
            with col2:
                st.metric(
                    "Selfed-like Fraction",
                    f"{result['selfed_like']:.2%}",
                )
            with col3:
                st.metric("Loci Compared", result["loci_compared"])

            if result["concordance"] > 0.95 and result["selfed_like"] < 0.05:
                st.success(
                    "Result is consistent with a true F1 hybrid"
                )
            elif result["selfed_like"] > 0.10:
                st.error(
                    "High selfed-like fraction - "
                    "possible selfing or contamination"
                )
            else:
                st.warning(
                    "Intermediate result - review manually"
                )

    # Tab 3: Seed Purity
    with tab3:
        st.header("Seed-Lot Purity Test")
        st.markdown(
            "Calculate off-type rate and confidence intervals "
            "from seed genotyping"
        )

        col1, col2, col3 = st.columns(3)
        with col1:
            n_seeds = st.number_input(
                "Seeds Tested", min_value=1, value=100
            )
        with col2:
            n_offtypes = st.number_input(
                "Off-types Detected",
                min_value=0,
                value=0,
                max_value=n_seeds,
            )
        with col3:
            confidence = st.selectbox(
                "Confidence", [0.90, 0.95, 0.99], index=1
            )

        result = seed_lot_purity_test(n_seeds, n_offtypes, confidence)

        st.subheader("Results")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric(
                "Observed Off-type Rate",
                f"{result['observed_rate']:.2%}",
            )
        with col2:
            st.metric(
                f"{confidence:.0%} Confidence Interval",
                f"[{result['ci_lower']:.2%}, {result['ci_upper']:.2%}]",
            )
        with col3:
            st.metric(
                "Seeds for 5% Detection",
                result["n_for_5pct"],
            )

        st.markdown(
            f"""
            **Interpretation:**

            - Tested {n_seeds} seeds, found {n_offtypes} off-types.
            - True off-type rate is between
              {result['ci_lower']:.2%} and {result['ci_upper']:.2%}
              with {confidence:.0%} confidence.
            - To detect a 5% off-type rate with 95% confidence,
              {result['n_for_5pct']} seeds are required.
            - To detect a 2% off-type rate with 95% confidence,
              {result['n_for_2pct']} seeds are required.
            """
        )

    # Tab 4: Design Assays
    with tab4:
        st.header("Design KASP Assays")
        st.markdown("Design primers for a panel of SNPs")

        panel_file = st.file_uploader(
            "Panel SNP List (TSV)", type=["tsv", "txt"]
        )
        fasta_file = st.text_input(
            "Reference FASTA Path", value="data/reference/B73_v5.fa"
        )

        if panel_file and st.button("Design Assays", type="primary"):
            panel_df = pd.read_csv(panel_file, sep="\t")
            st.write(
                f"Designing assays for {len(panel_df)} SNPs..."
            )

            with st.spinner("Designing primers..."):
                assays = design_panel_assays(
                    panel_df=panel_df,
                    fasta_path=fasta_file,
                    flank_length=60,
                )

            passed = sum(1 for a in assays if a.pass_all)
            st.success(
                f"Designed {len(assays)} assays, "
                f"{passed} passed all checks"
            )

            # Show results
            results_df = pd.DataFrame([
                {
                    "Variant": a.variant_id,
                    "Chrom": a.chrom,
                    "Pos": a.pos,
                    "Common Reverse": a.common_reverse,
                    "Allele 1 (FAM)": a.allele1_forward,
                    "Allele 2 (HEX)": a.allele2_forward,
                    "Tm (C)": f"{a.common_tm:.1f}/{a.allele1_tm:.1f}/{a.allele2_tm:.1f}",
                    "Amplicon": a.amplicon_size,
                    "Pass": "Passed" if a.pass_all else "Failed",
                }
                for a in assays
            ])

            st.dataframe(results_df, use_container_width=True)

            # Download buttons
            excel_buffer = BytesIO()
            write_kasp_design_sheet(assays, excel_buffer)
            st.download_button(
                "Download Design Sheet (Excel)",
                data=excel_buffer.getvalue(),
                file_name="kasp_design_sheet.xlsx",
                mime=(
                    "application/vnd.openxmlformats-"
                    "officedocument.spreadsheetml.sheet"
                ),
            )

    # Tab 5: Reference Panel
    with tab5:
        st.header("Reference Panel Overview")

        if "ref_loaded" in st.session_state:
            ref_gt = st.session_state["ref_gt"]
            ref_var_ids = st.session_state["ref_var_ids"]
            ref_sample_ids = st.session_state["ref_sample_ids"]

            st.metric("Varieties", len(np.unique(ref_sample_ids)))
            st.metric("SNPs in Panel", len(ref_var_ids))

            # Variety list
            var_counts = pd.Series(ref_sample_ids).value_counts()
            st.subheader("Varieties and Replicate Counts")
            st.dataframe(
                var_counts.reset_index().rename(
                    columns={"index": "Variety", "count": "Replicates"}
                ),
                use_container_width=True,
            )

            # SNP info
            st.subheader("Panel SNPs")
            snp_df = pd.DataFrame({
                "Variant_ID": ref_var_ids,
                "Index": range(len(ref_var_ids)),
            })
            st.dataframe(snp_df, use_container_width=True)
        else:
            st.info(
                "Click 'Load Demo Data' in the sidebar to see the "
                "reference panel, or load your own data."
            )


if __name__ == "__main__":
    main()