# =====================================================================
# Courbe KM combinée EXPLORATOIRE : OS AT, LT vs LR avec résection vasculaire
# (IPD reconstruites poolées)
#   - Dong 2025    : bras LR_avecVR (Dong_LR+VR_OS_AT.csv) vs LT ; LR sans VR écarté
#   - Poletto 2026 : bras LR (Poletto_LR_OS_AT.csv) vs LT
# NB : pooling naïf -> descriptif ; HR = Cox stratifié par étude (k = 2).
# Sortie : results/meta/KM_combinee_OS_AT_VR.png
# =====================================================================

library(survival)

HORIZON <- 60                       # mois (la plupart des courbes s'arrêtent à 5 ans)
COUL    <- c(LR = "#2166AC", LT = "#B2182B")

km_combinee <- function(ipd, etudes, png_out, titre, lab_LR = "LR") {
  ipd$groupe <- factor(ipd$groupe, levels = c("LR", "LT"))
  lab <- c(LT = "LT", LR = lab_LR)

  # Troncature administrative à 60 mois
  ipd$status <- ifelse(ipd$time > HORIZON, 0, ipd$status)
  ipd$time   <- pmin(ipd$time, HORIZON)

  # La stratification exige la colonne etude, avec une valeur par étude incluse
  stopifnot("etude" %in% names(ipd), !anyNA(ipd$etude),
            setequal(unique(ipd$etude), names(etudes)))

  fit <- survfit(Surv(time, status) ~ groupe, data = ipd)
  # Log-rank non stratifié (patients mélangés) : console seulement, pour comparaison
  p_lr <- 1 - pchisq(survdiff(Surv(time, status) ~ groupe, data = ipd)$chisq, df = 1)
  # Log-rank et Cox stratifiés par étude (affichés dans la figure)
  p_lr_strat <- 1 - pchisq(survdiff(Surv(time, status) ~ groupe + strata(etude), data = ipd)$chisq, df = 1)
  cx  <- coxph(Surv(time, status) ~ groupe + strata(etude), data = ipd)
  hr  <- exp(c(coef(cx), confint(cx)))
  p_cox <- summary(cx)$coefficients[1, "Pr(>|z|)"]
  fmt_p <- function(p) if (p < 1e-4) "p < 0.0001" else sprintf("p = %s", format(signif(p, 2), scientific = FALSE))

  s <- summary(fit, times = c(12, 36, 60), extend = TRUE)
  tab <- data.frame(groupe = sub("groupe=", "", s$strata), mois = s$time,
                    survie = round(100 * s$surv, 1),
                    IC_inf = round(100 * s$lower, 1), IC_sup = round(100 * s$upper, 1),
                    n_risque = s$n.risk)

  n_lignes <- ceiling(length(etudes) / 2)
  png(png_out, width = 2640, height = (1750 + 75 * n_lignes) * 1.2, res = 300)
  layout(matrix(1:3), heights = c(4, 1.25, 0.3 + 0.19 * n_lignes))
  par(mar = c(4, 10, 3, 5))
  plot(fit, col = COUL[c("LR", "LT")], lwd = 2, conf.int = FALSE, lty = 1,
       mark.time = TRUE, cex = 0.6, xmax = HORIZON, xaxt = "n", las = 1,
       xlab = "Months", ylab = "Overall survival", main = titre, cex.main = 1)
  axis(1, at = seq(0, HORIZON, 12))
  abline(v = c(12, 36, 60), col = "grey85", lty = 3)
  ns <- table(ipd$groupe); ev <- tapply(ipd$status, ipd$groupe, sum)
  legend("topright", bty = "n", lwd = 2, col = COUL[c("LT", "LR")],
         legend = c(sprintf("%s (n = %d, %d deaths)", lab["LT"], ns["LT"], ev["LT"]),
                    sprintf("%s (n = %d, %d deaths)", lab["LR"], ns["LR"], ev["LR"])))
  legend("bottomleft", bty = "n", cex = 0.8,
         legend = c(sprintf("HR %s vs %s (Cox stratified by study): %.2f [%.2f–%.2f]; %s",
                            lab["LT"], lab["LR"], hr[1], hr[2], hr[3], fmt_p(p_cox)),
                    sprintf("Stratified log-rank: %s", fmt_p(p_lr_strat))))

  # Table des patients à risque : bordures noires, une cellule par temps
  t_ar <- seq(0, HORIZON, 12)
  nr <- summary(fit, times = t_ar, extend = TRUE)
  par(mar = c(1, 10, 1.5, 5))
  plot(NA, xlim = c(0, HORIZON), ylim = c(0.5, 2.5), axes = FALSE, xlab = "", ylab = "")
  title(main = "Number at risk", cex.main = 0.85, adj = 0)
  for (i in 1:2) {
    g <- c("LT", "LR")[i]
    rect(c(-12, t_ar - 6), 3 - i - 0.5, c(-6, t_ar + 6), 3 - i + 0.5,
         border = "black", lwd = 1, xpd = NA)
    text(-9, 3 - i, lab[g], col = COUL[g], font = 2, xpd = NA, cex = if (nchar(lab[g]) > 3) 0.8 else 1)
    text(t_ar, 3 - i, nr$n.risk[nr$strata == paste0("groupe=", g)], col = COUL[g], cex = 0.9)
  }

  # Études incluses (n par bras), sur deux colonnes
  n_et <- table(factor(ipd$etude, levels = names(etudes)), ipd$groupe)
  lib  <- sprintf("%s %d (%s %d / %s %d)", names(etudes), etudes,
                  lab["LT"], n_et[, "LT"], lab["LR"], n_et[, "LR"])
  par(mar = c(0.5, 10, 1.5, 5))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  title(main = sprintf("Included studies (k = %d)", length(etudes)), cex.main = 0.85, adj = 0)
  y <- seq(0.85, 0.16, length.out = max(n_lignes, 2))[seq_len(n_lignes)]
  text(rep(c(0, 0.5), each = n_lignes)[seq_along(lib)], rep(y, 2)[seq_along(lib)], lib, adj = 0, cex = 0.8)
  dev.off()

  list(survie = tab, HR = hr, p_logrank = p_lr, p_logrank_strat = p_lr_strat, p_cox = p_cox)
}

etudes <- c(Dong = 2025, Poletto = 2026)
ipd_vr <- do.call(rbind, lapply(names(etudes), function(e)
  read.csv(file.path("data/IPD", paste0(e, "_OS_AT.csv")))))
ipd_vr <- subset(ipd_vr, bras %in% c("LT", "LR", "LR_avecVR"))   # Dong : retire LR_sansVR
print(table(ipd_vr$etude, ipd_vr$bras))

res <- km_combinee(ipd_vr, etudes, "results/meta/KM_combinee_OS_AT_VR.png",
                   "Pooled OS (as-treated, exploratory) – LT vs LR with vascular resection (2 studies)",
                   lab_LR = "LR + VR")

print(res$survie); print(round(res$HR, 2))
cat("\nComparaison des p-values :\n",
    "  log-rank non stratifié (patients mélangés) :", format.pval(res$p_logrank, digits = 3), "\n",
    "  log-rank stratifié par étude              :", format.pval(res$p_logrank_strat, digits = 3), "\n",
    "  Cox stratifié par étude (Wald)            :", format.pval(res$p_cox, digits = 3), "\n")
