# =====================================================================
# Forest plot OS AT : LT vs LR avec résection vasculaire (LR + VR)
#   - Dong 2025    : bras LR_avecVR (Dong_LR+VR_OS_AT.csv) vs LT ; LR sans VR écarté
#   - Poletto 2026 : bras LR (Poletto_LR_OS_AT.csv) vs LT
# Même présentation que forest_harmonise() de meta_analyse.R : HR par étude (Cox),
# k = 2 -> estimation poolée one-stage (Cox stratifié par étude), hétérogénéité sous la figure
# Sortie : results/meta/OS_AT_LT_vs_LR_VR.png
# =====================================================================

library(survival)
library(metafor)

etudes <- c(Dong = 2025, Poletto = 2026)
d <- do.call(rbind, lapply(names(etudes), function(e)
  read.csv(file.path("data/IPD", paste0(e, "_OS_AT.csv")))))
d <- subset(d, bras %in% c("LT", "LR", "LR_avecVR"))   # Dong : retire LR_sansVR
d$groupe <- factor(d$groupe, levels = c("LR", "LT"))   # LR + VR = référence -> HR LT vs LR + VR
print(table(d$etude, d$bras))

# HR par étude (Cox) + test des risques proportionnels
par_etude <- do.call(rbind, lapply(names(etudes), function(e) {
  x  <- d[d$etude == e, ]
  cx <- coxph(Surv(time, status) ~ groupe, data = x)
  data.frame(etude = e, annee = etudes[[e]],
             n_LT = sum(x$groupe == "LT"), n_LR = sum(x$groupe == "LR"),
             ev_LT = sum(x$status[x$groupe == "LT"]), ev_LR = sum(x$status[x$groupe == "LR"]),
             logHR = unname(coef(cx)), se = sqrt(vcov(cx))[1, 1],
             p_PH = cox.zph(cx)$table["groupe", "p"])
}))
par_etude$HR <- exp(par_etude$logHR)
par_etude$IC_inf <- exp(par_etude$logHR - 1.96 * par_etude$se)
par_etude$IC_sup <- exp(par_etude$logHR + 1.96 * par_etude$se)

# One-stage : Cox stratifié par étude (estimation poolée affichée, k = 2)
os <- coxph(Surv(time, status) ~ groupe + strata(etude), data = d)
one_stage <- c(exp(coef(os)), exp(confint(os)), summary(os)$coefficients[, "Pr(>|z|)"])

# Two-stage (pour I² / τ² / Q seulement)
ma <- rma(yi = logHR, sei = se, data = par_etude, method = "REML", test = "knha",
          slab = paste(etude, annee))

AXE_AT <- log(c(0.05, 0.1, 0.25, 0.5, 1, 2, 4))
k <- nrow(par_etude)
png("results/meta/OS_AT_LT_vs_LR_VR.png", width = 2000, height = 550 + 150 * k, res = 220)
par(mar = c(7, 4, 1, 2))
fp <- forest(ma, atransf = exp, refline = 0, addfit = FALSE, xlab = "HR (LT vs LR + VR)",
             alim = range(AXE_AT), at = AXE_AT, xlim = c(-8, 5), cex.axis = 0.75,
             rows = (k + 2):3, ylim = c(0.5, k + 5),
             header = c("Étude", "HR [IC 95 %]"))
abline(h = 2)
addpoly(log(one_stage[1]), ci.lb = log(one_stage[2]), ci.ub = log(one_stage[3]),
        rows = 1, atransf = exp, mlab = "One-stage (Cox stratifié)", cex = fp$cex)
p_q <- if (ma$QEp < 0.001) "< 0.001" else sprintf("= %.3f", ma$QEp)
mtext(sprintf("Hétérogénéité : I² = %.0f %% ; τ² = %.3f ; p (test Q) %s", ma$I2, ma$tau2, p_q),
      side = 1, line = 5.5, cex = 0.9)
dev.off()

print(within(par_etude, { logHR <- NULL; se <- NULL }), digits = 3)
cat(sprintf("\nOne-stage (Cox stratifié) : HR = %.2f [%.2f-%.2f], p = %.2g\n",
            one_stage[1], one_stage[2], one_stage[3], one_stage[4]))
cat(sprintf("Two-stage (REML, HK, 1 ddl) : HR = %.2f [%.2f-%.2f] ; I2 = %.0f %%, tau2 = %.3f, p(Q) = %.3f\n",
            exp(ma$b), exp(ma$ci.lb), exp(ma$ci.ub), ma$I2, ma$tau2, ma$QEp))
