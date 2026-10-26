# =====================================================================
# Fonctions de reconstruction des IPD (Guyot, IPDfromKM)
# Chargées par reconstruction_IPD.R et reconstruction_toutes.R
# =====================================================================

library(IPDfromKM)
library(readxl)
library(survival)
library(survminer)

# ---------------------------------------------------------------------
# Reconstruit UN bras à partir du CSV automeris (temps ; survie 0-1, décimale ",")
# et du fichier _atrisk.xlsx (colonnes time, nrisk). Tout est rendu en mois.
# unite = "mois" ou "annees" ; n_total : par défaut le N à risque au temps 0
# ---------------------------------------------------------------------
reconstruire <- function(csv, atrisk_xlsx, etude, critere, analyse, bras,
                         unite = "mois", n_total = NULL) {
  facteur <- if (unite == "annees") 12 else 1

  # 1. Points numérisés : nettoyage (temps >= 0, survie dans [0,1] et non croissante)
  pts <- read.table(csv, sep = ";", dec = ",", strip.white = TRUE, col.names = c("time", "surv"))
  pts$time <- pmax(pts$time, 0)
  pts <- pts[order(pts$time), ]
  pts$surv <- cummin(pmin(pmax(pts$surv, 0), 1))
  pts <- rbind(data.frame(time = 0, surv = 1), pts)

  # 2. Effectifs à risque (on garde ceux couverts par la courbe numérisée ;
  #    si la courbe s'arrête juste avant un temps du tableau, ex. 59.9 vs 60, on la prolonge)
  ar <- as.data.frame(read_excel(atrisk_xlsx))
  ar <- ar[!is.na(ar$nrisk), ]
  t_fin <- max(ar$time[ar$time <= max(pts$time) + 1])
  if (t_fin > max(pts$time)) pts <- rbind(pts, data.frame(time = t_fin, surv = tail(pts$surv, 1)))
  ar <- ar[ar$time <= t_fin, ]
  if (is.null(n_total)) n_total <- ar$nrisk[1]

  # 3. Algorithme de Guyot
  prep <- preprocess(dat = pts, trisk = ar$time, nrisk = ar$nrisk,
                     totalpts = n_total, maxy = 1)
  rec  <- getIPD(prep = prep, armID = 1)

  ipd <- data.frame(etude   = etude,
                    critere = critere,
                    analyse = analyse,
                    bras    = bras,
                    groupe  = if (grepl("^LT", bras)) "LT" else "LR",
                    time    = rec$IPD$time * facteur,   # tout en mois
                    status  = rec$IPD$status)
  # N à risque reconstruit aux temps du tableau, pour contrôle
  attr(ipd, "atrisk") <- data.frame(time = ar$time * facteur, publie = ar$nrisk,
                                    reconstruit = sapply(ar$time * facteur, function(t) sum(ipd$time >= t)))
  ipd
}

# ---------------------------------------------------------------------
# Valide une comparaison LT vs LR (à comparer à l'article)
# ---------------------------------------------------------------------
valider <- function(ipd, titre = "", fichier_png = NULL) {
  ipd$groupe <- factor(ipd$groupe, levels = c("LR", "LT"))
  fit <- survfit(Surv(time, status) ~ groupe, data = ipd)
  print(fit)                                                   # N, événements, médianes
  s <- summary(fit, times = c(12, 36, 60), extend = TRUE)
  print(data.frame(groupe = s$strata, mois = s$time, n_risque = s$n.risk,
                   survie = round(s$surv, 3)))
  cx <- coxph(Surv(time, status) ~ groupe, data = ipd)
  cat("HR LT vs LR :", round(exp(coef(cx)), 2),
      " IC95% :", round(exp(confint(cx)), 2), "\n")
  g <- ggsurvplot(fit, data = ipd, risk.table = TRUE, title = titre,
                  break.time.by = 12, xlab = "Mois")
  if (is.null(fichier_png)) print(g) else {
    png(fichier_png, width = 1800, height = 1400, res = 200); print(g); invisible(dev.off())
  }
  invisible(cx)
}
