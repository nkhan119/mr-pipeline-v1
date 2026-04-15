#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================
//  MR Pipeline — CDC 1.0.0  Phase B
//  Author : Nadeem Khan, Bioinformatician, INRS-CAFSB
//  GitHub : github.com/nkhan119
//
//  Stages:
//    B1  BUILD_RSID_MAP   — build CHR:BP → rsID lookup from bim files
//    B2  UNIVARIABLE_MR   — IVW + Egger + Weighted Median per pair × direction
//    B3  HETEROGENEITY    — Cochran Q · Egger intercept · I² · Steiger · PRESSO · LOO · single-SNP
//    B4  MVMR             — MVMR-IVW + MVMR-Egger for multi-exposure pairs
//    B5  MR_FIGURES       — publication-quality figures (300 DPI PNG + PDF)
//    B6  MR_REPORT        — interactive Plotly HTML report
//
//  Change (2026-04-13):
//    params.ref_panel is no longer passed as a `val` channel argument
//    to UNIVARIABLE_MR / HETEROGENEITY / MVMR. Each module resolves it
//    directly from params inside the script block. This avoids Nextflow
//    treating the path string as a file object and staging it into the
//    work directory (which broke PLINK bfile resolution on compute nodes).
// ================================================================

include { BUILD_RSID_MAP   } from './modules/build_rsid_map'
include { UNIVARIABLE_MR   } from './modules/univariable_mr'
include { HETEROGENEITY    } from './modules/heterogeneity'
include { MVMR             } from './modules/mvmr'
include { MR_FIGURES       } from './modules/mr_figures'
include { MR_REPORT        } from './modules/mr_report'

// ── Helper: parse "exp:out:p:label" list ──────────────────────
def parse_uvmr_pairs(str) {
    str.toString().split(",").collect { item ->
        def parts = item.trim().split(":")
        [exposure: parts[0], outcome: parts[1],
         p_exp: parts[2], label: parts[3]]
    }
}

def parse_mvmr_pairs(str) {
    str.toString().split(",").collect { item ->
        def parts = item.trim().split(":")
        [exposures: parts[0], p_exps: parts[1],
         outcome: parts[2], label: parts[3]]
    }
}

new File("${params.out_dir}/logs").mkdirs()

log.info """
╔══════════════════════════════════════════════════════════════════╗
║   MR Pipeline · CDC 1.0.0 · Phase B                            ║
║   Author    : ${params.author.padRight(52)}║
║   Institute : ${(params.institute ?: "").take(52).padRight(52)}║
╠══════════════════════════════════════════════════════════════════╣
║   sumstat_dir : ${params.sumstat_dir.padRight(49)}║
║   ref_bim_dir : ${(params.ref_bim_dir ?: "").padRight(49)}║
║   out_dir     : ${params.out_dir.padRight(49)}║
╠══════════════════════════════════════════════════════════════════╣
║   skip_uvmr   = ${params.skip_uvmr.toString().padRight(47)}║
║   skip_het    = ${params.skip_het.toString().padRight(47)}║
║   skip_mvmr   = ${params.skip_mvmr.toString().padRight(47)}║
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()

workflow {

    def cohort_list = params.cohorts.toString().split(",").collect { it.trim() }
    def uvmr_list   = parse_uvmr_pairs(params.uvmr_pairs)
    def mvmr_list   = parse_mvmr_pairs(params.mvmr_pairs)

    def directions = [
        [exp_c: cohort_list[0], out_c: cohort_list[1]],
        [exp_c: cohort_list[1], out_c: cohort_list[0]]
    ]

    // ── B1  rsID map ──────────────────────────────────────────────
    def bim_dir  = params.ref_bim_dir ?: "${params.base_dir}"
    def bim_file = file("${bim_dir}/1000G_phase3_common_norel.bim")
    if (!bim_file.exists())
        error "rsID reference bim not found: ${bim_file}\nCheck params.ref_bim_dir"

    ch_rsid_map = Channel.value(bim_file)

    // ── B2  Univariable MR ────────────────────────────────────────
    // Note: ref_panel no longer passed as a channel arg — resolved via
    //       params.ref_panel inside the module script block.
    if (!params.skip_uvmr) {

        ch_uvmr = Channel.from(uvmr_list)
            .combine(Channel.from(directions))
            .map { pair, dir ->
                tuple(
                    "${pair.label}__${dir.exp_c}_${dir.out_c}",
                    pair.exposure, pair.outcome, pair.p_exp,
                    dir.exp_c, dir.out_c
                )
            }

        UNIVARIABLE_MR(ch_uvmr, ch_rsid_map, params.sumstat_dir)

        ch_mr_tsv   = UNIVARIABLE_MR.out.mr_tsv.collect().ifEmpty([])
        ch_wald_tsv = UNIVARIABLE_MR.out.wald_tsv.collect().ifEmpty([])

    } else {
        log.info "⏭  skip_uvmr — loading existing MR results from disk"
        ch_mr_tsv   = Channel.fromPath(
            "${params.out_dir}/uvmr/*_mr.tsv").collect().ifEmpty([])
        ch_wald_tsv = Channel.fromPath(
            "${params.out_dir}/uvmr/*_wald.tsv").collect().ifEmpty([])
    }

    // ── B3  Heterogeneity ─────────────────────────────────────────
    // Note: ref_panel no longer passed as a channel arg.
    if (!params.skip_het) {

        ch_het_input = Channel.from(uvmr_list)
            .combine(Channel.from(directions))
            .map { pair, dir ->
                tuple(
                    "${pair.label}__${dir.exp_c}_${dir.out_c}",
                    pair.exposure, pair.outcome, pair.p_exp,
                    dir.exp_c, dir.out_c
                )
            }

        HETEROGENEITY(ch_het_input, ch_rsid_map, params.sumstat_dir)

        ch_het_tsv = HETEROGENEITY.out.het_tsv.collect().ifEmpty([])
        ch_loo_tsv = HETEROGENEITY.out.loo_tsv.collect().ifEmpty([])
        ch_snp_tsv = HETEROGENEITY.out.snp_tsv.collect().ifEmpty([])

    } else {
        log.info "⏭  skip_het — loading existing heterogeneity results"
        ch_het_tsv = Channel.fromPath(
            "${params.out_dir}/heterogeneity/*_het.tsv").collect().ifEmpty([])
        ch_loo_tsv = Channel.fromPath(
            "${params.out_dir}/heterogeneity/*_loo.tsv").collect().ifEmpty([])
        ch_snp_tsv = Channel.fromPath(
            "${params.out_dir}/heterogeneity/*_single.tsv").collect().ifEmpty([])
    }

    // ── B4  MVMR ──────────────────────────────────────────────────
    // Note: ref_panel no longer passed as a channel arg.
    if (!params.skip_mvmr) {

        ch_mvmr = Channel.from(mvmr_list)
            .combine(Channel.from(directions))
            .map { pair, dir ->
                tuple(
                    "${pair.label}__${dir.exp_c}_${dir.out_c}",
                    pair.exposures, pair.p_exps, pair.outcome,
                    dir.exp_c, dir.out_c
                )
            }

        MVMR(ch_mvmr, ch_rsid_map, params.sumstat_dir)

        ch_mvmr_tsv = MVMR.out.mvmr_tsv.collect().ifEmpty([])

    } else {
        log.info "⏭  skip_mvmr — loading existing MVMR results"
        ch_mvmr_tsv = Channel.fromPath(
            "${params.out_dir}/mvmr/*_mvmr.tsv").collect().ifEmpty([])
    }

    // ── B5  Figures ───────────────────────────────────────────────
    if (!params.skip_report) {
        MR_FIGURES(
            ch_mr_tsv, ch_wald_tsv, ch_het_tsv,
            ch_loo_tsv, ch_snp_tsv, ch_mvmr_tsv,
            params.author, params.institute, params.affiliation
        )
    }

    // ── B6  HTML Report ───────────────────────────────────────────
    if (!params.skip_report) {
        MR_REPORT(
            ch_mr_tsv, ch_wald_tsv, ch_het_tsv,
            ch_loo_tsv, ch_snp_tsv, ch_mvmr_tsv,
            params.author, params.affiliation,
            params.institute, params.github
        )
    }
}

workflow.onComplete {
    def status = workflow.success ? "SUCCESS ✓" : "FAILED  ✗"
    log.info """
╔══════════════════════════════════════════════════════════════════╗
║  MR Pipeline   ${status.padRight(47)}║
║  Duration    : ${workflow.duration.toString().padRight(49)}║
║  HTML Report : ${(params.out_dir + '/MR_Report.html').padRight(49)}║
║  Figures     : ${(params.out_dir + '/figures/').padRight(49)}║
╚══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
