// ================================================================
//  modules/mvmr.nf — Multivariable MR
// ================================================================

process MVMR {
    label 'medium'
    tag   "${analysis_id}"
    publishDir "${params.out_dir}/mvmr", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(analysis_id),
          val(exposures_str), val(p_exps_str), val(outcome),
          val(exp_cohort), val(out_cohort)
    path  rsid_map
    val   sumstat_dir

    output:
    path "${analysis_id}_mvmr.tsv", emit: mvmr_tsv, optional: true
    path "${analysis_id}_mvmr.log", emit: mvmr_log

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
      library(TwoSampleMR)
      library(MVMR)
      library(MendelianRandomization)
      library(data.table)
    })

    if (nchar("${params.opengwas_token}") > 0) {
      Sys.setenv(OPENGWAS_JWT = "${params.opengwas_token}")
    }

    aid       <- "${analysis_id}"
    exposures <- strsplit("${exposures_str}", "+", fixed = TRUE)[[1]]
    p_exps    <- as.numeric(strsplit("${p_exps_str}", "+", fixed = TRUE)[[1]])
    outcome   <- "${outcome}"
    exp_coh   <- "${exp_cohort}"
    out_coh   <- "${out_cohort}"
    ss_dir    <- "${sumstat_dir}"
    ref_panel <- "${params.ref_panel}"
    rsid_f    <- "${rsid_map}"

    log_f <- file(paste0(aid, "_mvmr.log"), open = "wt")
    lg <- function(m) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), m)
      cat(msg, "\\n", file = log_f, append = TRUE)
      cat(msg, "\\n")
    }
    lg(paste("START MVMR:", aid))

    # ── rsID mapping ─────────────────────────────────────────────
    bim <- fread(rsid_f, header = FALSE, data.table = TRUE,
                 col.names = c("CHR","rsID","cM","BP","A1","A2"))
    bim[, CHR := as.integer(CHR)][, BP := as.integer(BP)]
    bim <- bim[grepl("^rs", rsID)]
    bim <- unique(bim, by = c("CHR", "BP"))

    map_rsid <- function(dat) {
      if (all(c("chr.exposure","pos.exposure") %in% names(dat))) {
        dat\$CHR_n <- as.integer(dat\$chr.exposure)
        dat\$BP_n  <- as.integer(dat\$pos.exposure)
      } else {
        sp <- strsplit(dat\$SNP, ":", fixed = TRUE)
        dat\$CHR_n <- suppressWarnings(as.integer(sapply(sp, `[`, 1)))
        dat\$BP_n  <- suppressWarnings(as.integer(sapply(sp, `[`, 2)))
      }
      dt <- as.data.table(dat)
      m  <- merge(dt, bim[, .(CHR, BP, rsID)],
                  by.x = c("CHR_n","BP_n"), by.y = c("CHR","BP"), all.x = TRUE)
      m  <- m[!is.na(rsID)]
      m[, SNP := rsID][, c("CHR_n","BP_n","rsID") := NULL]
      as.data.frame(m)
    }

    load_fmt <- function(cohort, trait, pth = NULL, is_exp = TRUE) {
      path <- file.path(ss_dir, cohort, "sumstats",
                        sprintf("%s_%s_sumstats.tsv.gz", cohort, trait))
      if (!file.exists(path)) {
        lg(paste("Missing sumstats:", path))
        return(NULL)
      }
      dt <- fread(path, data.table = FALSE)
      for (c in c("BETA","SE","P","EAF","N","CHR","BP"))
        if (c %in% names(dt)) dt[[c]] <- as.numeric(as.character(dt[[c]]))

      if (is_exp && !is.null(pth))
        dt <- dt[!is.na(dt\$P) & dt\$P < pth & !is.na(dt\$BETA) & !is.na(dt\$SE), ]

      if (nrow(dt) == 0) return(NULL)

      tp <- if (is_exp) "exposure" else "outcome"
      fd <- format_data(dat = dt, type = tp,
                        snp_col = "SNP", beta_col = "BETA", se_col = "SE",
                        effect_allele_col = "A1", other_allele_col = "A2",
                        eaf_col = "EAF", pval_col = "P", samplesize_col = "N",
                        chr_col = "CHR", pos_col = "BP",
                        phenotype_col = trait)
      fd[[tp]] <- trait
      fd
    }

    # ── Clumping (server → local PLINK) ─────────────────────────
    clump <- function(dat) {
      cl <- tryCatch(
        clump_data(dat, clump_r2 = ${params.clump_r2}, clump_kb = ${params.clump_kb}, pop = "EUR"),
        error = function(e) NULL
      )
      if (!is.null(cl) && nrow(cl) > 0) {
        lg(sprintf("  Server clumping: %d instruments retained", nrow(cl)))
        return(cl)
      }

      lg("  Server clumping failed — trying local PLINK...")
      pb <- trimws(system("which plink2 2>/dev/null || which plink", intern = TRUE)[1])
      if (length(pb) == 0 || is.na(pb)) { lg("  PLINK not found"); return(NULL) }

      tmp_prefix <- tempfile(pattern = "clump_mvmr_")
      inp_file   <- paste0(tmp_prefix, "_input.txt")
      write.table(data.frame(SNP = dat\$SNP, P = dat\$pval.exposure),
                  inp_file, row.names = FALSE, quote = FALSE, sep = "\\t")

      rc <- system2(pb, args = c("--bfile", ref_panel, "--clump", inp_file,
                                 "--clump-kb", ${params.clump_kb},
                                 "--clump-p1", 1, "--clump-r2", ${params.clump_r2},
                                 "--out", tmp_prefix), stdout = FALSE, stderr = FALSE)

      if (rc != 0) { lg(sprintf("  PLINK exited with %d", rc)); return(NULL) }

      res_file <- if (file.exists(paste0(tmp_prefix, ".clumps"))) paste0(tmp_prefix, ".clumps") else
                  if (file.exists(paste0(tmp_prefix, ".clumped"))) paste0(tmp_prefix, ".clumped") else NULL

      if (is.null(res_file)) { lg("  No clumping output"); return(NULL) }

      clumped <- tryCatch(fread(res_file, data.table = FALSE), error = function(e) NULL)
      if (is.null(clumped) || nrow(clumped) == 0) return(NULL)

      id_col <- if ("ID" %in% names(clumped)) "ID" else "SNP"
      keep <- dat[dat\$SNP %in% clumped[[id_col]], ]
      lg(sprintf("  Local clumping: %d instruments retained", nrow(keep)))
      keep
    }

    # ── Main MVMR pipeline (identical logic to your working script) ─────
    exp_clumped_list <- lapply(seq_along(exposures), function(i) {
      d <- load_fmt(exp_coh, exposures[i], p_exps[i], TRUE)
      if (is.null(d) || nrow(d) == 0) return(NULL)
      d <- map_rsid(d)
      if (nrow(d) == 0) return(NULL)
      cl <- clump(d)
      if (is.null(cl) || nrow(cl) == 0) return(NULL)
      cl\$f_stat <- (cl\$beta.exposure / cl\$se.exposure)^2
      cl[cl\$f_stat > ${params.f_stat_min}, ]
    })

    if (any(sapply(exp_clumped_list, is.null))) {
      lg("SKIP: missing instruments for one or more exposures")
      close(log_f); quit("no")
    }

    all_snps <- unique(unlist(lapply(exp_clumped_list, `[[`, "SNP")))
    lg(sprintf("  Union instrument set: %d SNPs", length(all_snps)))
    if (length(all_snps) < 3) {
      lg("SKIP: <3 union SNPs"); close(log_f); quit("no")
    }

    exp_full_list <- lapply(exposures, function(trait) {
      d <- load_fmt(exp_coh, trait, NULL, TRUE)
      if (is.null(d)) return(NULL)
      d <- map_rsid(d)
      d[d\$SNP %in% all_snps, ]
    })

    if (any(sapply(exp_full_list, function(x) is.null(x) || nrow(x) == 0))) {
      lg("SKIP: failed to extract full data at union SNPs")
      close(log_f); quit("no")
    }

    out_dat <- load_fmt(out_coh, outcome, NULL, FALSE)
    if (is.null(out_dat)) {
      lg("SKIP: no outcome data"); close(log_f); quit("no")
    }
    out_dat <- map_rsid(out_dat)
    out_sub <- out_dat[out_dat\$SNP %in% all_snps, ]

    common_snps <- Reduce(intersect, c(lapply(exp_full_list, `[[`, "SNP"), list(out_sub\$SNP)))
    lg(sprintf("  %d SNPs common across all datasets", length(common_snps)))
    if (length(common_snps) < 3) {
      lg("SKIP: <3 common SNPs"); close(log_f); quit("no")
    }

    align_snps <- function(df, snps) df[match(snps, df\$SNP), ]
    exp_al <- lapply(exp_full_list, align_snps, common_snps)
    out_al <- align_snps(out_sub, common_snps)

    betaX   <- do.call(cbind, lapply(exp_al, `[[`, "beta.exposure"))
    betaXse <- do.call(cbind, lapply(exp_al, `[[`, "se.exposure"))
    betaY   <- out_al\$beta.outcome
    betaYse <- out_al\$se.outcome

    colnames(betaX) <- colnames(betaXse) <- exposures

    keep <- complete.cases(betaX, betaXse, betaY, betaYse)
    if (sum(keep) < 3) {
      lg("SKIP: <3 complete cases"); close(log_f); quit("no")
    }

    betaX   <- betaX[keep, , drop = FALSE]
    betaXse <- betaXse[keep, , drop = FALSE]
    betaY   <- betaY[keep]
    betaYse <- betaYse[keep]
    snps_used <- common_snps[keep]

    lg(sprintf("  Running MVMR: %d SNPs x %d exposures", nrow(betaX), ncol(betaX)))

    lbl_dir <- sprintf("%s->%s", exp_coh, out_coh)
    lbl_exp <- paste(exposures, collapse = "+")

    # Format for MVMR package
    mvmr_fmt <- tryCatch({
      MVMR::format_mvmr(
        BXGs   = as.matrix(betaX),
        BYG    = as.numeric(betaY),
        seBXGs = as.matrix(betaXse),
        seBYG  = as.numeric(betaYse),
        RSID   = as.character(snps_used)
      )
    }, error = function(e) {
      lg(paste("  format_mvmr failed:", e\$message))
      NULL
    })

    if (is.null(mvmr_fmt)) {
      close(log_f); quit("no")
    }

    # MVMR-IVW using MVMR::mvmr()  (exactly as in your working script)
    mvmr_ivw <- tryCatch({
      res <- MVMR::mvmr(r_input = mvmr_fmt, gencov = 0)

      if (!inherits(res, "MVMRIVW"))
        stop("mvmr() did not return MVMRIVW object")

      coef_m <- res\$coef
      f_vals <- tryCatch(as.numeric(res\$Q_strength), error = function(e) rep(NA_real_, length(exposures)))

      data.frame(
        analysis_id = aid,
        label       = lbl_exp,
        outcome     = outcome,
        direction   = lbl_dir,
        method      = "MVMR-IVW",
        exposure    = exposures,
        beta        = as.numeric(coef_m[,1]),
        se          = as.numeric(coef_m[,2]),
        pval        = as.numeric(coef_m[,4]),
        n_snps      = nrow(betaX),
        cond_F      = f_vals,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      lg(paste("  MVMR-IVW failed:", e\$message))
      NULL
    })

    # MVMR-Egger using MendelianRandomization::mr_mvegger()
    mvmr_egger <- tryCatch({
      mv_input <- MendelianRandomization::mr_mvinput(
        bx = as.matrix(betaX), bxse = as.matrix(betaXse),
        by = as.numeric(betaY), byse = as.numeric(betaYse),
        snps = as.character(snps_used)
      )

      egger <- MendelianRandomization::mr_mvegger(mv_input)

      data.frame(
        analysis_id = aid,
        label       = lbl_exp,
        outcome     = outcome,
        direction   = lbl_dir,
        method      = "MVMR-Egger",
        exposure    = exposures,
        beta        = as.numeric(egger@Estimate),
        se          = as.numeric(egger@StdError.Est),
        pval        = as.numeric(egger@Pvalue.Est),
        n_snps      = nrow(betaX),
        cond_F      = NA_real_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      lg(paste("  MVMR-Egger failed:", e\$message))
      NULL
    })

    # Write results
    res_list <- Filter(Negate(is.null), list(mvmr_ivw, mvmr_egger))
    if (length(res_list) > 0) {
      final <- do.call(rbind, res_list)
      write.table(final, paste0(aid, "_mvmr.tsv"),
                  sep = "\\t", row.names = FALSE, quote = FALSE)
      lg(sprintf("  SUCCESS: Wrote %d MVMR estimates", nrow(final)))
    } else {
      lg("  WARNING: No MVMR results generated")
    }

    close(log_f)
    """
}
