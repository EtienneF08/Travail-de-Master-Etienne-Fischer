# =====================================================================
# Méta-analyse LT vs LR (pCCA) à partir des IPD reconstruites
# Entrée : data/IPD_tous.csv (etude, critere, analyse, bras, groupe, time, status)
# Sortie : results/meta/ (tableaux + forest plots)
# =====================================================================

# install.packages(c("survival", "metafor", "survRM2", "coxme", "dplyr"))
library(survival)
library(metafor)
library(survRM2)
library(coxme)
library(dplyr)

dir.create("results/meta", recursive = TRUE, showWarnings = FALSE)

ipd <- read.csv("data/IPD_tous.csv")
# groupe = LT / LR (les 2 courbes LR de Dong AT sont déjà regroupées en LR)
if (!"groupe" %in% names(ipd)) ipd$groupe <- ifelse(ipd$bras == "LT", "LT", "LR")
ipd$groupe <- factor(ipd$groupe, levels = c("LR", "LT"))   # LR = référence -> HR LT vs LR

# Année de publication (CodeBook) : tri et étiquettes « Étude année » des forest plots
cb <- as.data.frame(suppressMessages(readxl::read_excel("data/CodeBook.1.def..xlsx", "Extraction",
                                                        col_names = FALSE, col_types = "text")))
cb <- setNames(cb[3:nrow(cb), ], unlist(cb[1, ]))
nom_cb <- c(Giuseppucci = "A.Giuseppuci", Poletto = "E.Poletto", Hidalgo = "E.Hidalgo",
            Ethun = "Ethnun CG.", Dong = "Y.Dong", Breuer = "E.Breuer", Croome = "P.Croome",
            Iwatsuki = "S.Iwatsuki", Robles = "R.Robles", Pichlmayr = "R.Pichlmayr")
annee <- setNames(as.integer(cb$year[match(nom_cb, cb$Study)]), names(nom_cb))

# ---------------------------------------------------------------------
# Les analyses (principales et de sensibilité)
# Croome chevauche Dong (Mayo Rochester) : jamais les deux ensemble
# ---------------------------------------------------------------------
AT_tous <- c("Dong", "Ethun", "Poletto", "Hidalgo", "Breuer", "Giuseppucci",
             "Iwatsuki", "Robles", "Pichlmayr")
analyses <- list(
  "OS ITT - principal (Dong + Ethun Fig 4A)" = list(critere = "OS",  analyse = "ITT", etudes = c("Dong", "Ethun")),
  "OS ITT - sensibilité sans Ethun"         = list(critere = "OS",  analyse = "ITT", etudes = c("Dong")),
  "OS ITT - sensibilité Croome au lieu de Dong" = list(critere = "OS", analyse = "ITT", etudes = c("Croome", "Ethun")),
  "OS AT - principal"                       = list(critere = "OS",  analyse = "AT",  etudes = AT_tous),
  "OS AT - sensibilité Croome au lieu de Dong" = list(critere = "OS", analyse = "AT", etudes = c(setdiff(AT_tous, "Dong"), "Croome")),
  "OS AT - sensibilité sans Robles"         = list(critere = "OS",  analyse = "AT",  etudes = setdiff(AT_tous, "Robles")),
  "RFS AT - principal"                      = list(critere = "RFS", analyse = "AT",  etudes = c("Dong", "Poletto", "Breuer", "Giuseppucci")),
  "RFS AT - sensibilité Croome au lieu de Dong" = list(critere = "RFS", analyse = "AT", etudes = c("Croome", "Poletto", "Breuer", "Giuseppucci"))
)

# ---------------------------------------------------------------------
# Étape 1 : HR par étude (Cox) + test des risques proportionnels
# ---------------------------------------------------------------------
hr_par_etude <- function(d) {
  d %>% group_by(etude) %>% group_modify(function(x, k) {
    cx  <- coxph(Surv(time, status) ~ groupe, data = x)
    zph <- cox.zph(cx)$table["groupe", "p"]
    data.frame(n_LT = sum(x$groupe == "LT"), n_LR = sum(x$groupe == "LR"),
               ev_LT = sum(x$status[x$groupe == "LT"]), ev_LR = sum(x$status[x$groupe == "LR"]),
               logHR = coef(cx), se = sqrt(vcov(cx))[1, 1], p_PH = zph)
  }) %>% ungroup() %>% mutate(HR = exp(logHR), IC_inf = exp(logHR - 1.96 * se),
                              IC_sup = exp(logHR + 1.96 * se))
}

# ---------------------------------------------------------------------
# Forest plot harmonisé : études triées par année, même axe X partout,
# une seule estimation poolée : k = 2 -> one-stage Cox stratifié (HK à 1 ddl inutilisable),
# k >= 3 -> two-stage effets aléatoires (REML, HK) ; hétérogénéité sous la figure
# ---------------------------------------------------------------------
AXE_AT <- log(c(0.05, 0.1, 0.25, 0.5, 1, 2, 4))
forest_harmonise <- function(par_etude, one_stage, fichier) {
  k  <- nrow(par_etude)
  pe <- par_etude %>% mutate(annee = annee[etude], lab = paste(etude, annee)) %>% arrange(annee, etude)
  ma <- rma(yi = logHR, sei = se, data = pe, method = "REML", test = "knha", slab = lab)

  png(fichier, width = 2000, height = 550 + 150 * k, res = 220)
  par(mar = c(7, 4, 1, 2))
  fp <- forest(ma, atransf = exp, refline = 0, addfit = FALSE, xlab = "HR (LT vs LR)",
               alim = range(AXE_AT), at = AXE_AT, xlim = c(-8, 5), cex.axis = 0.75,
               rows = (k + 2):3, ylim = c(0.5, k + 5),
               header = c("Étude", "HR [IC 95 %]"))
  abline(h = 2)
  if (k == 2) {
    addpoly(log(one_stage[1]), ci.lb = log(one_stage[2]), ci.ub = log(one_stage[3]),
            rows = 1, atransf = exp, mlab = "One-stage (Cox stratifié)", cex = fp$cex)
  } else {
    addpoly(ma, row = 1, atransf = exp, mlab = "Two-stage effets aléatoires (REML, HK)", cex = fp$cex)
  }
  p_q <- if (ma$QEp < 0.001) "< 0.001" else sprintf("= %.3f", ma$QEp)
  mtext(sprintf("Hétérogénéité : I² = %.0f %% ; τ² = %.3f ; p (test Q) %s", ma$I2, ma$tau2, p_q),
        side = 1, line = 5.5, cex = 0.9)
  dev.off()
  ma
}

# ---------------------------------------------------------------------
# Boucle sur toutes les analyses
# ---------------------------------------------------------------------
resume <- list()
hr_etudes <- list()
for (nom in names(analyses)) {
  a <- analyses[[nom]]
  d <- ipd %>% filter(critere == a$critere, analyse == a$analyse, etude %in% a$etudes)
  etudes_presentes <- unique(d$etude)
  manquantes <- setdiff(a$etudes, etudes_presentes)
  cat("\n=====================================================\n", nom,
      "\n=====================================================\n")
  if (length(manquantes)) cat("Absentes des IPD :", paste(manquantes, collapse = ", "), "\n")

  # 1. HR par étude
  par_etude <- hr_par_etude(d)
  print(as.data.frame(par_etude %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
  k <- nrow(par_etude)
  hr_etudes[[nom]] <- data.frame(analyse = nom, par_etude)

  # 2. One-stage : Cox stratifié par étude (calculé d'abord : utilisé sous le forest plot si k = 2)
  os <- coxph(Surv(time, status) ~ groupe + strata(etude), data = d)
  one_stage <- c(exp(coef(os)), exp(confint(os)), summary(os)$coefficients[, "Pr(>|z|)"])

  # 3. Two-stage : méta-analyse des log(HR), REML + Hartung-Knapp
  #    Le tableau garde les deux méthodes pour toutes les analyses ; pour k = 2, HK n'a
  #    qu'un degré de liberté (IC très larges) : le forest plot montre alors le one-stage
  if (k >= 2) {
    ma <- forest_harmonise(par_etude, one_stage,
      file.path("results/meta", paste0(gsub("[^A-Za-z0-9]+", "_", chartr("éèêàç", "eeeac", nom)), ".png")))
    print(ma)
    two_stage <- c(exp(ma$b), exp(ma$ci.lb), exp(ma$ci.ub), ma$pval, ma$I2, ma$tau2)
  } else {
    cat("Une seule étude : pas de pooling, on rapporte son HR.\n")
    two_stage <- c(par_etude$HR, par_etude$IC_inf, par_etude$IC_sup, NA, NA, NA)
  }

  # 4. One-stage avec effet aléatoire de l'effet traitement si k >= 3
  # Pente aléatoire sur l'indicateur numérique LT (avec le facteur groupe, coxme crée
  # 2 effets aléatoires corrélés LR/LT non identifiables avec strata(etude) -> échec)
  one_stage_re <- rep(NA, 4)
  if (k >= 3) {
    d$lt <- as.numeric(d$groupe == "LT")
    me <- tryCatch(coxme(Surv(time, status) ~ lt + strata(etude) + (0 + lt | etude), data = d),
                   error = function(e) { cat("coxme a échoué :", conditionMessage(e), "\n"); NULL })
    if (!is.null(me)) {
      b <- fixef(me); se_b <- sqrt(vcov(me))[1, 1]
      one_stage_re <- c(exp(b), exp(b - 1.96 * se_b), exp(b + 1.96 * se_b), VarCorr(me)$etude[[1]])
      cat("One-stage coxme (effet aléatoire LT vs LR par étude) : HR =", round(one_stage_re[1], 2),
          "[", round(one_stage_re[2], 2), ";", round(one_stage_re[3], 2), "], tau2 =",
          round(one_stage_re[4], 3), "\n")
    }
  }

  resume[[nom]] <- data.frame(
    analyse = nom, k = k, N = nrow(d),
    HR_2stage = two_stage[1], IC_inf_2stage = two_stage[2], IC_sup_2stage = two_stage[3],
    p_2stage = two_stage[4], I2 = two_stage[5], tau2 = two_stage[6],
    HR_1stage = one_stage[1], IC_inf_1stage = one_stage[2], IC_sup_1stage = one_stage[3],
    p_1stage = one_stage[4],
    HR_1stage_RE = one_stage_re[1], IC_inf_1stage_RE = one_stage_re[2],
    IC_sup_1stage_RE = one_stage_re[3], tau2_1stage_RE = one_stage_re[4],
    PH_violee = paste(par_etude$etude[par_etude$p_PH < 0.05], collapse = ", "))
}

tableau <- bind_rows(resume)
tableau_aff <- tableau %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(tableau_aff)
write.csv(tableau_aff, "results/meta/resume_meta_analyses.csv", row.names = FALSE)

# HR de chaque étude, pour chaque analyse (à comparer aux HR publiés)
hr_tab <- bind_rows(hr_etudes) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(hr_tab, "results/meta/HR_par_etude.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# Étape 4 : RMST à 36 et 60 mois (utile si les risques ne sont pas proportionnels)
# Différence de RMST LT - LR par étude, puis poolée
# Une étude dont le bras le plus court s'arrête à <= TOL_RMST mois de l'horizon est gardée,
# avec sa RMST calculée jusqu'à sa fin de suivi (tau = min(horizon, fin de suivi))
# ---------------------------------------------------------------------
TOL_RMST <- 1.5
rmst_meta <- function(critere_, analyse_, etudes, tau) {
  d <- ipd %>% filter(critere == critere_, analyse == analyse_, etude %in% etudes)
  res <- d %>% group_by(etude) %>% group_modify(function(x, k) {
    fin <- min(tapply(x$time, x$groupe, max))            # fin de suivi du bras le plus court
    if (fin < tau - TOL_RMST) return(data.frame())       # suivi trop court
    tau_e <- min(tau, fin)
    r <- rmst2(x$time, x$status, as.integer(x$groupe == "LT"), tau = tau_e)
    u <- r$unadjusted.result[1, ]   # RMST (LT - LR)
    data.frame(tau_utilise = tau_e, diff = u[1], se = (u[3] - u[2]) / (2 * 1.96))
  }) %>% ungroup()
  exclues <- setdiff(etudes, res$etude)
  tronquees <- res %>% filter(tau_utilise < tau)
  if (nrow(res) < 2) return(NULL)
  ma <- rma(yi = diff, sei = se, data = res, method = "REML", test = "knha", slab = etude)
  cat("\nRMST", critere_, analyse_, "à", tau, "mois : différence LT - LR =",
      round(ma$b, 1), "mois [", round(ma$ci.lb, 1), ";", round(ma$ci.ub, 1), "], k =", nrow(res), "\n")
  if (nrow(tronquees)) cat("  Gardées avec tau = fin de suivi :",
      paste0(tronquees$etude, " (", round(tronquees$tau_utilise, 1), " m)", collapse = ", "), "\n")
  if (length(exclues)) cat("  Exclues (suivi trop court) :", paste(exclues, collapse = ", "), "\n")
  print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
  data.frame(critere = critere_, analyse = analyse_, tau = tau, k = nrow(res),
             diff_RMST = ma$b[1], IC_inf = ma$ci.lb, IC_sup = ma$ci.ub, p = ma$pval, I2 = ma$I2,
             etudes = paste(res$etude, collapse = ", "),
             tau_raccourci = paste0(tronquees$etude, if (nrow(tronquees)) " (" else "",
                                    if (nrow(tronquees)) round(tronquees$tau_utilise, 1) else "",
                                    if (nrow(tronquees)) " m)" else "", collapse = ", "),
             exclues_suivi_court = paste(exclues, collapse = ", "))
}
rmst_tab <- list()
for (tau in c(36, 60)) {
  rmst_tab[[length(rmst_tab) + 1]] <- rmst_meta("OS",  "ITT", c("Dong", "Ethun"), tau)
  rmst_tab[[length(rmst_tab) + 1]] <- rmst_meta("OS",  "AT",  AT_tous, tau)
  rmst_tab[[length(rmst_tab) + 1]] <- rmst_meta("RFS", "AT",  c("Dong", "Poletto", "Breuer", "Giuseppucci"), tau)
}
rmst_tab <- bind_rows(rmst_tab) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(rmst_tab, "results/meta/RMST.csv", row.names = FALSE)
