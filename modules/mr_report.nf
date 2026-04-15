// ================================================================
//  modules/mr_report.nf — Stage B6
//  Interactive Plotly HTML MR report
// ================================================================
process MR_REPORT {
    label 'small'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path mr_tsvs
    path wald_tsvs
    path het_tsvs
    path loo_tsvs
    path snp_tsvs
    path mvmr_tsvs
    val  author
    val  affiliation
    val  institute
    val  github

    output:
    path "MR_Report.html", emit: report

    script:
    """
    python3 ${projectDir}/bin/render_mr_report.py \
        --template    ${projectDir}/assets/mr_report_template.html \
        --author      "${author}" \
        --affiliation "${affiliation}" \
        --institute   "${institute}" \
        --github      "${github}"
    """
}
