// ================================================================
//  modules/heterogeneity.nf — Stage B3
//
//  Per-analysis heterogeneity & pleiotropy diagnostics:
//    1. Cochran Q            — instrument heterogeneity (IVW + Egger)
//    2. MR-Egger intercept   — directional pleiotropy test
//    3. I² statistic         — heterogeneity magnitude (Bowden formula)
//    4. Steiger filtering    — confirms exposure → outcome direction
//    5. MR-PRESSO            — outlier-robust estimate (if installed)
//    6. Leave-one-out        — IVW stability when each IV dropped
//    7. Single-SNP MR        — per-instrument Wald ratios
//
//  Change (2026-04-13):
//    val ref_panel removed from input block; resolved via params.ref_panel
//    inside the script. clump() rewritten to call PLINK directly via
//    system2() and handle both .clumps (PLINK2) and .clumped (PLINK 1.9).
// ================================================================
process HETEROGENEITY {
    label 'medium'
    tag   "${analysis_id}"
    publishDir "${params.out_dir}/heterogeneity", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(analysis_id),
          val(exposure), val(outcome), val(p_exp),
          val(exp_cohort), val(out_cohort)
    path  rsid_map
    val   sumstat_dir

    output:
    path "${analysis_id}_het.tsv",    emit: het_tsv,    optional: true
    path "${analysis_id}_loo.tsv",    emit: loo_tsv,    optional: true
    path "${analysis_id}_single.tsv", emit: snp_tsv,    optional: true
    path "${analysis_id}_het.log",    emit: het_log

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
      library(TwoSampleMR); library(data.table); library(ieugwasr)
    })
    presso_ok <- requireNamespace("MRPRESSO", quietly = TRUE)

    if (nchar("${params.opengwas_token}") > 0)
      Sys.setenv(OPENGWAS_JWT = "${params.opengwas_token}")

    aid       <- "${analysis_id}"
    exposure  <- "${exposure}"
    outcome   <- "${outcome}"
    p_exp     <- as.numeric("${p_exp}")
    exp_coh   <- "${exp_cohort}"
    out_coh   <- "${out_cohort}"
    ss_dir    <- "${sumstat_dir}"
    ref_panel <- "${params.ref_panel}"
    rsid_f    <- "${rsid_map}"

    log_f <- file(paste0(aid, "_het.log"), open = "wt")
    lg <- function(m) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), m)
      cat(msg, "\\n", file = log_f, append = TRUE); cat(msg, "\\n")
    }
    lg(paste("START HETEROGENEITY:", aid))

    # ── rsID lookup ────────────────────────────────────────────
    bim <- fread(rsid_f, header = FALSE, data.table = TRUE,
                 col.names = c("CHR","rsID","cM","BP","A1","A2"))
    bim[, CHR := as.integer(CHR)][, BP := as.integer(BP)]
    bim <- bim[grepl("^rs", rsID)]
    bim <- unique(bim, by = c("CHR","BP"))

    map_rsid <- function(dat, lbl = "") {
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
      m  <- m[!is.na(rsID)]; m[, SNP := rsID][, c("CHR_n","BP_n","rsID") := NULL]
      as.data.frame(m)
    }

    load_fmt <- function(cohort, trait, pth = NULL, is_exp = TRUE) {
      path <- file.path(ss_dir, cohort, "sumstats",
                        sprintf("%s_%s_sumstats.tsv.gz", cohort, trait))
      if (!file.exists(path)) return(NULL)
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
                        chr_col = "CHR", pos_col = "BP")
      fd[[tp]] <- trait; fd
    }

    # ── LD clumping (server → direct PLINK fallback) ───────────
    clump <- function(dat) {
      cl <- tryCatch(
        clump_data(dat, clump_r2 = ${params.clump_r2},
                   clump_kb = ${params.clump_kb}, pop = "EUR"),
        error = function(e) NULL)
      if (!is.null(cl) && nrow(cl) > 0) {
        lg(sprintf("  Server clumping: %d instruments retained", nrow(cl)))
        return(cl)
      }

      lg("  Server clumping failed — trying local PLINK...")

      pb <- tryCatch({
        candidates <- trimws(system(
          "which plink2 2>/dev/null; which plink 2>/dev/null", intern = TRUE))
        candidates <- candidates[nchar(candidates) > 0]
        if (length(candidates) == 0) stop("no plink binary found")
        candidates[1]
      }, error = function(e) {
        lg(paste("  PLINK not found:", e\$message)); return(NULL)
      })
      if (is.null(pb)) return(NULL)

      tmp_prefix <- tempfile(pattern = "clump_het_")
      inp_file   <- paste0(tmp_prefix, "_input.txt")
      write.table(
        data.frame(SNP = dat\$SNP, P = dat\$pval.exposure,
                   stringsAsFactors = FALSE),
        inp_file, row.names = FALSE, quote = FALSE, sep = "\\t")

      rc <- system2(pb, args = c(
        "--bfile",    ref_panel,
        "--clump",    inp_file,
        "--clump-kb", ${params.clump_kb},
        "--clump-p1", 1,
        "--clump-r2", ${params.clump_r2},
        "--out",      tmp_prefix
      ), stdout = FALSE, stderr = FALSE)

      if (rc != 0) { lg(sprintf("  PLINK exited with code %d", rc)); return(NULL) }

      out_clumps  <- paste0(tmp_prefix, ".clumps")
      out_clumped <- paste0(tmp_prefix, ".clumped")
      result_file <- if (file.exists(out_clumps))  out_clumps  else
                     if (file.exists(out_clumped)) out_clumped else NULL

      if (is.null(result_file) || file.size(result_file) == 0) {
        lg("  Local clump: no output file produced"); return(NULL)
      }

      clumped <- tryCatch(fread(result_file, data.table = FALSE),
                          error = function(e) NULL)
      if (is.null(clumped) || nrow(clumped) == 0) {
        lg("  Local clump: output empty or unreadable"); return(NULL)
      }

      id_col <- if ("ID"  %in% names(clumped)) "ID"  else
                if ("SNP" %in% names(clumped)) "SNP" else NULL
      if (is.null(id_col)) {
        lg(sprintf("  Local clump: unrecognised columns: %s",
                   paste(names(clumped), collapse = ", ")))
        return(NULL)
      }

      keep <- dat[dat\$SNP %in% clumped[[id_col]], ]
      lg(sprintf("  Local clumping: %d instruments retained", nrow(keep)))
      keep
    }

    # ── Load, map, clump, harmonise ────────────────────────────
    exp_dat <- load_fmt(exp_coh, exposure, p_exp, is_exp = TRUE)
    if (is.null(exp_dat)) { lg("SKIP: no exposure"); close(log_f); quit("no") }
    exp_dat <- map_rsid(exp_dat, exposure)
    if (nrow(exp_dat) == 0) { lg("SKIP: 0 rsID"); close(log_f); quit("no") }

    exp_cl <- clump(exp_dat)
    if (is.null(exp_cl) || nrow(exp_cl) == 0) {
      lg("SKIP: clumping = 0"); close(log_f); quit("no") }

    exp_cl\$f_stat <- (exp_cl\$beta.exposure / exp_cl\$se.exposure)^2
    exp_cl <- exp_cl[exp_cl\$f_stat > ${params.f_stat_min}, ]
    if (nrow(exp_cl) < 3) {
      lg(sprintf("SKIP: only %d F>${params.f_stat_min} instruments", nrow(exp_cl)))
      close(log_f); quit("no")
    }
    lg(sprintf("  %d instruments (F>${params.f_stat_min})", nrow(exp_cl)))

    out_dat <- load_fmt(out_coh, outcome, is_exp = FALSE)
    if (is.null(out_dat)) { lg("SKIP: no outcome"); close(log_f); quit("no") }
    out_dat <- map_rsid(out_dat, outcome)

    dat_h <- tryCatch(harmonise_data(exp_cl, out_dat, action = 2),
                      error = function(e) NULL)
    if (is.null(dat_h)) { lg("SKIP: harmonise failed"); close(log_f); quit("no") }
    dat_h <- dat_h[dat_h\$mr_keep, ]
    if (nrow(dat_h) < 3) {
      lg(sprintf("SKIP: %d SNPs post-harmonise (need >= 3)", nrow(dat_h)))
      close(log_f); quit("no")
    }
    lg(sprintf("  %d SNPs after harmonisation", nrow(dat_h)))

    n_snps   <- nrow(dat_h)
    lbl_pair <- paste(exposure, outcome, sep = "->")
    lbl_dir  <- sprintf("%s->%s", exp_coh, out_coh)

    het_rows <- list()

    # ── 1. Cochran Q ───────────────────────────────────────────
    q_res <- tryCatch(mr_heterogeneity(dat_h), error = function(e) NULL)
    if (!is.null(q_res)) {
      for (i in seq_len(nrow(q_res))) {
        r  <- q_res[i, ]
        i2 <- max(0, 100 * (r\$Q - r\$Q_df) / r\$Q)
        het_rows[[length(het_rows)+1]] <- data.frame(
          analysis_id = aid, causal_pair = lbl_pair, direction = lbl_dir,
          test = "Cochran_Q", method = r\$method,
          statistic = r\$Q, df = r\$Q_df, pvalue = r\$Q_pval,
          i2 = i2, n_snps = n_snps,
          note = "p<0.05 = heterogeneous instruments",
          stringsAsFactors = FALSE)
      }
      lg(sprintf("  Cochran Q: Q=%.3f df=%d p=%.4f I2=%.1f%%",
                 q_res\$Q[1], q_res\$Q_df[1], q_res\$Q_pval[1],
                 max(0, 100*(q_res\$Q[1]-q_res\$Q_df[1])/q_res\$Q[1])))
    }

    # ── 2. MR-Egger intercept ──────────────────────────────────
    eg_res <- tryCatch(mr_pleiotropy_test(dat_h), error = function(e) NULL)
    if (!is.null(eg_res)) {
      het_rows[[length(het_rows)+1]] <- data.frame(
        analysis_id = aid, causal_pair = lbl_pair, direction = lbl_dir,
        test = "Egger_intercept", method = "MR Egger",
        statistic = eg_res\$egger_intercept, df = NA_real_,
        pvalue = eg_res\$pval, i2 = NA_real_, n_snps = n_snps,
        note = sprintf("Intercept=%.4f SE=%.4f (p<0.05 = directional pleiotropy)",
                       eg_res\$egger_intercept, eg_res\$se),
        stringsAsFactors = FALSE)
      lg(sprintf("  Egger intercept: %.4f (SE=%.4f) p=%.4f",
                 eg_res\$egger_intercept, eg_res\$se, eg_res\$pval))
    }

    # ── 3. I² statistic (Bowden et al) ─────────────────────────
    b     <- dat_h\$beta.exposure
    se_b  <- dat_h\$se.exposure
    bhat  <- sum(b/se_b^2) / sum(1/se_b^2)
    Q_isq <- sum((b - bhat)^2 / se_b^2)
    df_i  <- n_snps - 1
    isq   <- max(0, 100 * (Q_isq - df_i) / Q_isq)
    het_rows[[length(het_rows)+1]] <- data.frame(
      analysis_id = aid, causal_pair = lbl_pair, direction = lbl_dir,
      test = "I2", method = "All instruments",
      statistic = isq, df = df_i,
      pvalue = pchisq(Q_isq, df_i, lower.tail = FALSE),
      i2 = isq, n_snps = n_snps,
      note = "0-24%: low; 25-49%: moderate; 50-74%: substantial; >=75%: high",
      stringsAsFactors = FALSE)
    lg(sprintf("  I² = %.1f%%", isq))

    # ── 4. Steiger filtering ───────────────────────────────────
    stg <- tryCatch({
      dat_stg <- steiger_filtering(dat_h)
      n_pass  <- sum(dat_stg\$steiger_dir == TRUE,  na.rm = TRUE)
      n_fail  <- sum(dat_stg\$steiger_dir == FALSE, na.rm = TRUE)
      data.frame(
        analysis_id = aid, causal_pair = lbl_pair, direction = lbl_dir,
        test = "Steiger", method = "Steiger filtering",
        statistic = n_pass, df = n_fail, pvalue = NA_real_,
        i2 = NA_real_, n_snps = n_snps,
        note = sprintf("Pass=%d Fail=%d | Direction confirmed=%s",
                       n_pass, n_fail,
                       ifelse(n_pass > n_fail, "YES", "NO — check reverse causation")),
        stringsAsFactors = FALSE)
    }, error = function(e) {
      lg(sprintf("  Steiger error: %s", e\$message)); NULL
    })
    if (!is.null(stg)) {
      het_rows[[length(het_rows)+1]] <- stg
      lg(sprintf("  Steiger: %s", stg\$note))
    }

    # ── 5. MR-PRESSO (optional) ────────────────────────────────
    if (presso_ok && n_snps >= 5) {
      presso <- tryCatch({
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome   = "se.outcome",   SdExposure   = "se.exposure",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
          data = dat_h,
          NbDistribution = ${params.presso_n_dist},
          SignifThreshold = 0.05)
      }, error = function(e) {
        lg(sprintf("  MR-PRESSO error: %s", e\$message)); NULL
      })
      if (!is.null(presso)) {
        gt <- presso\$`MR-PRESSO results`\$`Global Test`
        het_rows[[length(het_rows)+1]] <- data.frame(
          analysis_id = aid, causal_pair = lbl_pair, direction = lbl_dir,
          test = "MRPRESSO_Global", method = "MR-PRESSO",
          statistic = gt\$RSSobs, df = NA_real_, pvalue = gt\$Pvalue,
          i2 = NA_real_, n_snps = n_snps,
          note = "Global test: p<0.05 = outlier instruments present",
          stringsAsFactors = FALSE)
        lg(sprintf("  MR-PRESSO Global: RSSobs=%.3f p=%.4f",
                   gt\$RSSobs, gt\$Pvalue))
        mr_main <- presso\$`Main MR results`
        if (!is.null(mr_main) && nrow(mr_main) >= 2) {
          corr <- mr_main[2, ]
          het_rows[[length(het_rows)+1]] <- data.frame(
            analysis_id = aid, causal_pair = lbl_pair, direction = lbl_dir,
            test = "MRPRESSO_Corrected", method = "MR-PRESSO",
            statistic = as.numeric(corr[,"Causal Estimate"]),
            df = NA_real_,
            pvalue = as.numeric(corr[,"P-value"]),
            i2 = NA_real_, n_snps = n_snps,
            note = sprintf("Outlier-corrected beta=%.4f",
                           as.numeric(corr[,"Causal Estimate"])),
            stringsAsFactors = FALSE)
        }
      }
    } else if (!presso_ok) {
      lg("  MR-PRESSO not installed — skipping (optional)")
      lg("  Install: devtools::install_github('rondolab/MR-PRESSO')")
    }

    # ── Write heterogeneity TSV ────────────────────────────────
    if (length(het_rows) > 0) {
      het_df <- do.call(rbind, het_rows)
      write.table(het_df, paste0(aid, "_het.tsv"),
                  sep = "\\t", row.names = FALSE, quote = FALSE)
      lg(sprintf("  Heterogeneity: %d tests written", nrow(het_df)))
    }

    # ── 6. Leave-one-out ───────────────────────────────────────
    loo <- tryCatch(mr_leaveoneout(dat_h), error = function(e) NULL)
    if (!is.null(loo) && nrow(loo) > 0) {
      loo\$analysis_id <- aid
      loo\$causal_pair <- lbl_pair
      loo\$direction   <- lbl_dir
      write.table(loo, paste0(aid, "_loo.tsv"),
                  sep = "\\t", row.names = FALSE, quote = FALSE)
      lg(sprintf("  LOO: %d rows", nrow(loo)))
    }

    # ── 7. Single-SNP MR ──────────────────────────────────────
    snp_mr <- tryCatch(mr_singlesnp(dat_h), error = function(e) NULL)
    if (!is.null(snp_mr) && nrow(snp_mr) > 0) {
      snp_mr\$analysis_id <- aid
      snp_mr\$causal_pair <- lbl_pair
      snp_mr\$direction   <- lbl_dir
      f_vec <- c((dat_h\$beta.exposure/dat_h\$se.exposure)^2[
                    match(snp_mr\$SNP, dat_h\$SNP)])
      snp_mr\$f_stat <- f_vec
      write.table(snp_mr, paste0(aid, "_single.tsv"),
                  sep = "\\t", row.names = FALSE, quote = FALSE)
      lg(sprintf("  Single-SNP MR: %d rows", nrow(snp_mr)))
    }

    lg("DONE")
    close(log_f)
    """
}
