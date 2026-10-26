# Reconstruction des IPD (méthode de Guyot, package IPDfromKM) – OS transplantation (LT) vs résection (LR)
# Source : figure 1 (Kaplan-Meier, axe X en années 0-10, Y en %), capture data/Capture d'écran 2026-10-12 à 17.01.09.png
library(png); library(IPDfromKM); library(survival)

img <- list.files("data", pattern = "17\\.01\\.09\\.png$", full.names = TRUE)
im  <- readPNG(img); D <- (im[,,1] + im[,,2] + im[,,3]) / 3 < 0.5

# Calibration des axes (graduations) : 0 an = col 269.5, 1 an = 78.75 px ; 0 % = ligne 569.5, 1 % = 5.04 px
X0 <- 269.5; PXY <- 78.75; Y0 <- 569.5; PXP <- 5.04

runs <- function(c) {
  r <- which(D[55:566, c]) + 54
  if (!length(r)) return(NULL)
  b <- cumsum(c(1, diff(r) > 1)); t(sapply(split(r, b), range))
}
# Suivi d'une courbe de droite à gauche (la survie ne peut que remonter en allant vers t = 0)
trace <- function(c_start, row_start, c_end = 271, cap_before = 0) {
  prev <- row_start; out <- NULL
  for (c in c_start:c_end) {
    R <- runs(c)
    if (!is.null(R)) {
      touch <- which(R[, 1] <= prev + 3 & R[, 2] >= prev - 9)
      if (length(touch)) {
        r <- R[touch[1], ]
        lvl <- if (diff(r) <= 8) mean(r) else r[1] + 2.5
        if (c < cap_before) lvl <- max(lvl, prev - 6)  # zone de chevauchement : remontée limitée à 6 px/colonne
        prev <- min(prev, lvl)
      }
    }
    out <- rbind(out, c(c, prev))
  }
  d <- data.frame(t = (out[, 1] - X0) / PXY * 12, s = (Y0 - out[, 2]) / PXP)
  d <- d[order(d$t), ]; d$s <- pmin(100, cummin(d$s))
  rbind(data.frame(t = 0, s = 100), d, data.frame(t = 120, s = tail(d$s, 1)))
}
lt <- trace(round(X0 + 8 * PXY), 492)  # trait-point (LT), départ à 8 ans / 15.4 %
lr <- trace(1055, 492, cap_before = X0 + 1.7 * PXY)  # trait plein épais (résection), départ à 10 ans ;
# avant 1.7 an les deux courbes se chevauchent : on limite la remontée pour ne pas suivre les marches de LT

# Effectifs à risque (tableau sous la figure), temps en mois
trisk <- seq(0, 120, 12)
nr_lr <- c(136, 93, 60, 36, 31, 26, 22, 18, 9, 6, 4)
nr_lt <- c(19, 15, 8, 5, 4, 4, 3, 3, 3, 3, 2)

ipd <- function(d, nr, arm) {
  p <- preprocess(dat = d, trisk = trisk, nrisk = nr, totalpts = nr[1], maxy = 100)
  x <- getIPD(prep = p, armID = arm)$IPD; names(x) <- c("time", "status", "arm"); x
}
IPD <- rbind(ipd(lr, nr_lr, 0), ipd(lt, nr_lt, 1))
IPD$groupe <- factor(IPD$arm, 0:1, c("Resection", "Transplantation"))
write.csv(IPD, "results/IPD_OS_LT_vs_LR.csv", row.names = FALSE)

# Contrôles : survie reconstruite à 12 / 36 / 60 mois et effectifs à risque
fit <- survfit(Surv(time, status) ~ groupe, data = IPD)
print(summary(fit, times = c(12, 36, 60)))
print(summary(fit, times = trisk)$n.risk)

# HR (Cox) LT vs résection + test du log-rank
cox <- coxph(Surv(time, status) ~ groupe, data = IPD)
print(summary(cox))
print(survdiff(Surv(time, status) ~ groupe, data = IPD))
print(cox.zph(cox))
s <- summary(cox)
cat(sprintf("\nHR (LT vs LR) = %.2f [IC95 %.2f-%.2f], p = %.3f ; logHR = %.4f, seTE = %.4f\n",
    s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5], coef(cox), sqrt(vcov(cox))))

png("results/KM_reconstruit_OS_LT_vs_LR.png", width = 1600, height = 1100, res = 200)
plot(fit, lty = c(1, 4), lwd = c(3, 2), xlab = "Mois", ylab = "Survie globale", xmax = 120, mark.time = TRUE)
lines(lr$t, lr$s / 100, type = "s", col = "grey60"); lines(lt$t, lt$s / 100, type = "s", col = "grey60")
legend("topright", c("Résection (IPD reconstruits)", "Transplantation (IPD reconstruits)", "Courbes numérisées"),
       lty = c(1, 4, 1), lwd = c(3, 2, 1), col = c(1, 1, "grey60"), bty = "n")
invisible(dev.off())
