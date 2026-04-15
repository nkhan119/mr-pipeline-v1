// ================================================================
//  modules/univariable_mr.nf — Stage B2
//  Univariable MR: IVW + MR-Egger + Weighted Median
//  per causal pair × direction
//
//  Fix (2026-04-13): PLINK2 writes .clumps (not .clumped).
//  ld_clump() from TwoSampleMR hardcodes .clumped, so we bypass
//  it entirely and call PLINK directly via system2(), then read
//  the output ourselves — works for both PLINK 1.9 and PLINK2.
// ================================================================
process UNIVARIABLE_MR {
    label 'medium'
    tag   "${analysis_id}"
    publishDir "${params.out_dir}/uvmr", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(analysis_id),
          val(exposure), val(outcome), val(p_exp),
          val(exp_cohort), val(out_cohort)
    path  rsid_map
    val   sumstat_dir

    output:
    path "${analysis_id}_mr.tsv",   emit: mr_tsv,   optional: true
    path "${analysis_id}_wald.tsv", emit: wald_tsv, optional: true
    path "${analysis_id}_mr.log",   emit: log_file

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
      library(TwoSampleMR); library(data.table); library(ieugwasr)
    })

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

    log_f <- file(paste0(aid, "_mr.log"), open = "wt")
    lg <- function(m) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), m)
      cat(msg, "\\n", file = log_f, append = TRUE); cat(msg, "\\n")
    }
    lg(paste("START uvMR:", aid))

    # ── rsID lookup ────────────────────────────────────────────
    bim <- fread(rsid_f, header = FALSE, data.table = TRUE,
                 col.names = c("CHR","rsID","cM","BP","A1","A2"))
    bim[, CHR := as.integer(CHR)][, BP := as.integer(BP)]
    bim <- bim[grepl("^rs", rsID)]
    bim <- unique(bim, by = c("CHR","BP"))
    lg(sprintf("  rsID map: %s entries", format(nrow(bim), big.mark = ",")))

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
      lg(sprintf("  rsID [%s]: %d / %d mapped", lbl, sum(!is.na(m\$rsID)), nrow(dat)))
      m <- m[!is.na(rsID)]; m[, SNP := rsID][, c("CHR_n","BP_n","rsID") := NULL]
      as.data.frame(m)
    }

    # ── Load & format ──────────────────────────────────────────
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

    # ── LD clumping ────────────────────────────────────────────
    # Bypasses TwoSampleMR::ld_clump() entirely because it hardcodes
    # the .clumped extension which PLINK2 no longer writes (.clumps).
    # Strategy: server first, then direct PLINK system2() call with
    # manual output file handling for both PLINK 1.9 and PLINK2.
    clump <- function(dat) {

      # 1. Try OpenGWAS server clumping
      cl <- tryCatch(
        clump_data(dat, clump_r2 = ${params.clump_r2},
                   clump_kb = ${params.clump_kb}, pop = "EUR"),
        error = function(e) NULL)
      if (!is.null(cl) && nrow(cl) > 0) {
        lg(sprintf("  Server clumping: %d instruments retained", nrow(cl)))
        return(cl)
      }

      lg("  Server clumping failed — trying local PLINK...")

      # 2. Locate PLINK binary (prefer plink2, fall back to plink)
      pb <- tryCatch({
        candidates <- trimws(system(
          "which plink2 2>/dev/null; which plink 2>/dev/null",
          intern = TRUE))
        candidates <- candidates[nchar(candidates) > 0]
        if (length(candidates) == 0) stop("no plink binary found")
        candidates[1]
      }, error = function(e) {
        lg(paste("  PLINK not found:", e\$message)); return(NULL)
      })
      if (is.null(pb)) return(NULL)

      lg(sprintf("  Using PLINK binary: %s", pb))

      # 3. Write SNP / p-value input file
      tmp_prefix <- tempfile(pattern = "clump_mr_")
      inp_file   <- paste0(tmp_prefix, "_input.txt")
      write.table(
        data.frame(SNP = dat\$SNP, P = dat\$pval.exposure,
                   stringsAsFactors = FALSE),
        inp_file, row.names = FALSE, quote = FALSE, sep = "\\t")

      # 4. Run PLINK
      rc <- system2(pb, args = c(
        "--bfile",    ref_panel,
        "--clump",    inp_file,
        "--clump-kb", ${params.clump_kb},
        "--clump-p1", 1,
        "--clump-r2", ${params.clump_r2},
        "--out",      tmp_prefix
      ), stdout = FALSE, stderr = FALSE)

      if (rc != 0) {
        lg(sprintf("  PLINK exited with code %d", rc))
        return(NULL)
      }

      # 5. Handle .clumps (PLINK2) vs .clumped (PLINK 1.9)
      out_clumps  <- paste0(tmp_prefix, ".clumps")   # PLINK2
      out_clumped <- paste0(tmp_prefix, ".clumped")  # PLINK 1.9

      result_file <- if (file.exists(out_clumps))  out_clumps  else
                     if (file.exists(out_clumped)) out_clumped else NULL

      if (is.null(result_file) || file.size(result_file) == 0) {
        lg("  Local clump: no output file produced"); return(NULL)
      }

      clumped <- tryCatch(fread(result_file, data.table = FALSE),
                          error = function(e) NULL)
      if (is.null(clumped) || nrow(clumped) == 0) {
        lg("  Local clump: output file empty or unreadable"); return(NULL)
      }

      # PLINK2 uses 'ID'; PLINK 1.9 uses 'SNP'
      id_col <- if ("ID" %in% names(clumped)) "ID" else
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

    # ── Main ──────────────────────────────────────────────────
    exp_dat <- load_fmt(exp_coh, exposure, p_exp, is_exp = TRUE)
    if (is.null(exp_dat)) { lg("SKIP: no exposure"); close(log_f); quit("no") }
    exp_dat <- map_rsid(exp_dat, exposure)
    if (nrow(exp_dat) == 0) { lg("SKIP: 0 rsID"); close(log_f); quit("no") }

    lg(sprintf("  Clumping %d variants...", nrow(exp_dat)))
    exp_cl <- clump(exp_dat)
    if (is.null(exp_cl) || nrow(exp_cl) == 0) {
      lg("SKIP: clumping = 0"); close(log_f); quit("no") }

    exp_cl\$f_stat <- (exp_cl\$beta.exposure / exp_cl\$se.exposure)^2
    exp_cl <- exp_cl[exp_cl\$f_stat > ${params.f_stat_min}, ]
    if (nrow(exp_cl) == 0) {
      lg("SKIP: all F <= ${params.f_stat_min}"); close(log_f); quit("no") }
    lg(sprintf("  %d instruments (F > ${params.f_stat_min})", nrow(exp_cl)))

    out_dat <- load_fmt(out_coh, outcome, is_exp = FALSE)
    if (is.null(out_dat)) { lg("SKIP: no outcome"); close(log_f); quit("no") }
    out_dat <- map_rsid(out_dat, outcome)

    dat_h <- tryCatch(harmonise_data(exp_cl, out_dat, action = 2),
                      error = function(e) NULL)
    if (is.null(dat_h)) { lg("SKIP: harmonise failed"); close(log_f); quit("no") }
    dat_h <- dat_h[dat_h\$mr_keep, ]
    if (nrow(dat_h) < ${params.min_instruments}) {
      lg(sprintf("SKIP: %d SNPs post-harmonise", nrow(dat_h)))
      close(log_f); quit("no")
    }
    lg(sprintf("  %d SNPs after harmonisation", nrow(dat_h)))

    mr_res <- tryCatch(
      mr(dat_h, method_list = c("mr_ivw","mr_egger_regression",
                                "mr_weighted_median")),
      error = function(e) NULL)

    if (!is.null(mr_res) && nrow(mr_res) > 0) {
      mr_res\$causal_pair <- paste(exposure, outcome, sep = "->")
      mr_res\$exposure    <- exposure
      mr_res\$outcome     <- outcome
      mr_res\$direction   <- sprintf("%s->%s", exp_coh, out_coh)
      mr_res\$n_snps      <- nrow(dat_h)
      write.table(mr_res, paste0(aid, "_mr.tsv"),
                  sep = "\\t", row.names = FALSE, quote = FALSE)

      wald <- data.frame(
        analysis_id = aid,
        SNP         = dat_h\$SNP,
        exposure    = exposure,
        outcome     = outcome,
        direction   = sprintf("%s->%s", exp_coh, out_coh),
        beta_exp    = dat_h\$beta.exposure,
        beta_out    = dat_h\$beta.outcome,
        se_exp      = dat_h\$se.exposure,
        se_out      = dat_h\$se.outcome,
        wald_ratio  = dat_h\$beta.outcome / dat_h\$beta.exposure,
        f_stat      = (dat_h\$beta.exposure / dat_h\$se.exposure)^2,
        stringsAsFactors = FALSE)
      write.table(wald, paste0(aid, "_wald.tsv"),
                  sep = "\\t", row.names = FALSE, quote = FALSE)
      lg(sprintf("SUCCESS: %d estimates written", nrow(mr_res)))
    } else {
      lg("WARNING: mr() returned no results")
    }
    close(log_f)
    """
}
