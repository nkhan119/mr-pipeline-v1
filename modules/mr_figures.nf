// ================================================================
//  modules/mr_figures.nf — Stage B5
//  Publication-quality MR figures (300 DPI PNG + PDF)
//  No grid lines · generous whitespace · Nature-style axes
// ================================================================
process MR_FIGURES {
    label 'medium'
    publishDir "${params.out_dir}/figures", mode: 'copy'

    input:
    path mr_tsvs
    path wald_tsvs
    path het_tsvs
    path loo_tsvs
    path snp_tsvs
    path mvmr_tsvs
    val  author
    val  institute
    val  affiliation

    output:
    path "F1_forest_uvmr.png",        optional: true, emit: f1_png
    path "F1_forest_uvmr.pdf",        optional: true, emit: f1_pdf
    path "F2_heterogeneity.png",      optional: true, emit: f2_png
    path "F2_heterogeneity.pdf",      optional: true, emit: f2_pdf
    path "F3_scatter_funnel.png",     optional: true, emit: f3_png
    path "F3_scatter_funnel.pdf",     optional: true, emit: f3_pdf
    path "F4_leave_one_out.png",      optional: true, emit: f4_png
    path "F4_leave_one_out.pdf",      optional: true, emit: f4_pdf
    path "F5_single_snp.png",         optional: true, emit: f5_png
    path "F5_single_snp.pdf",         optional: true, emit: f5_pdf
    path "F6_mvmr_forest.png",        optional: true, emit: f6_png
    path "F6_mvmr_forest.pdf",        optional: true, emit: f6_pdf

    script:
    """
    python3 ${projectDir}/bin/mr_figures.py \
        --mr_dir      . \
        --out_dir     . \
        --author      "${author}" \
        --institute   "${institute}" \
        --affiliation "${affiliation}"

    echo "[FIGURES] Complete:"
    ls -lh F*.png F*.pdf 2>/dev/null || echo "  No figure files produced"
    """
}
