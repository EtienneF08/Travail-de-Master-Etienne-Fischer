# Hong et al. 2011, Arch Surg 146:683-689 - RFS des tumeurs hilaires (figure 3B)
# OLT n = 13 vs LR (radical resection) n = 7, temps en mois depuis la chirurgie.
# Endpoint de l'article : "tumor recurrence-free patient survival" (récidive ou décès).
#
# Fichiers (data/) :
#   Hong2011_fig3B_coords_OLT.csv / _LR.csv : courbes extraites du tracé vectoriel du PDF,
#     colonnes time_months, survival (0-1)
#   Hong2011_RFS_IPD_reconstruit.csv : IPD reconstruites (algorithme de Guyot 2012),
#     colonnes study, arm, time_months, event (1 = récidive/décès, 0 = censuré)
#
# Limites : pas de tableau des patients à risque ni de nombre total d'événements
# dans l'article, donc aucune censure n'est placée avant la fin du suivi.
# A lancer depuis la racine du projet : Rscript script/hong2011_RFS.R

library(survival)
library(IPDfromKM)

dir <- "data"
ipd <- read.csv(file.path(dir, "Hong2011_RFS_IPD_reconstruit.csv"))
ipd$arm <- factor(ipd$arm, levels = c("LR", "OLT"))
coords <- lapply(c(OLT = "OLT", LR = "LR"), function(a)
  read.csv(file.path(dir, paste0("Hong2011_fig3B_coords_", a, ".csv"))))

# 1. Courbes reconstruites
fit <- survfit(Surv(time_months, event) ~ arm, data = ipd)
print(summary(fit, times = c(12, 36, 60), extend = TRUE))
cat("Publié (lu sur les courbes) à 12 / 36 / 60 mois :\n")
for (a in names(coords)) {
  co <- coords[[a]]
  cat(a, sapply(c(12, 36, 60), function(m) round(co$survival[max(which(co$time_months <= m))], 3)), "\n")
}

# 2. Superposition avec les courbes extraites de l'article (contrôle visuel)
png("results/Hong2011_RFS_reconstruit_vs_publie.png", width = 1600, height = 1100, res = 200)
plot(fit, lty = c(2, 1), xlab = "Mois", ylab = "Survie sans récidive",
     xlim = c(0, 60), main = "Hong 2011, hilaires : reconstruit vs publié")
for (a in names(coords))
  lines(coords[[a]]$time_months, coords[[a]]$survival, type = "s", col = "red",
        lty = ifelse(a == "OLT", 1, 2))
legend("topright", c("LR reconstruit", "OLT reconstruit", "Publié (rouge)"),
       lty = c(2, 1, 1), col = c("black", "black", "red"), bty = "n")
invisible(dev.off())

# 3. Log-rank (publié : p = 0.62) et HR de Cox OLT vs LR
print(survdiff(Surv(time_months, event) ~ arm, data = ipd))
cox <- coxph(Surv(time_months, event) ~ arm, data = ipd)
print(summary(cox))

# 4. Variante : reconstruction indépendante avec IPDfromKM (contrôle)
n <- c(OLT = 13, LR = 7)
ipd2 <- do.call(rbind, lapply(names(n), function(a) {
  pre <- preprocess(dat = coords[[a]], trisk = 0, nrisk = n[[a]], totalpts = n[[a]], maxy = 1)
  rec <- getIPD(prep = pre, armID = match(a, names(n)))$IPD
  data.frame(arm = a, time_months = rec[[1]], event = rec[[2]])
}))
ipd2$arm <- factor(ipd2$arm, levels = c("LR", "OLT"))
print(table(ipd2$arm, ipd2$event))
print(survdiff(Surv(time_months, event) ~ arm, data = ipd2))
cox2 <- coxph(Surv(time_months, event) ~ arm, data = ipd2)

# 5. Résumé pour la méta-analyse (HR OLT vs LR)
res <- function(m, lab) { s <- summary(m)
  cat(sprintf("%-22s HR = %.2f [IC95 %.2f-%.2f], p = %.3f ; logHR = %.4f, seTE = %.4f\n", lab,
      s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5], coef(m), sqrt(vcov(m)))) }
res(cox, "IPD fournies :"); res(cox2, "IPDfromKM :")
