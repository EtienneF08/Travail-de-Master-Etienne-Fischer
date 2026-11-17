# =====================================================================
# Forest plot OS ITT ALL : analyse OS ITT principale (Dong + Ethun Fig 4A) + Vibert
# Même procédé que meta_analyse.R (HR par étude, one-stage, two-stage, coxme,
# forest_harmonise), avec le one-stage Cox stratifié comme estimation poolée
# pour k <= 3 (même règle que meta_analyse_sans_Poletto.R)
# Entrées : data/IPD_tous.csv + data/IPD/Vibert_OS_ITT.csv (script/reconstruction_Vibert.R)
# Sortie : results/meta/OS_ITT_ALL.png (n'écrase aucun fichier existant)
# =====================================================================

# Préparation de meta_analyse.R (données, années, hr_par_etude, forest_harmonise),
# sans lancer sa boucle d'analyses ni réécrire ses tableaux
expr <- parse("script/meta_analyse.R", keep.source = FALSE)
fin_prep <- which(sapply(expr, function(e) identical(e, quote(resume <- list())))) - 1
for (e in expr[seq_len(fin_prep)]) eval(e)

code_fp <- deparse(forest_harmonise)
stopifnot(sum(grepl("k == 2", code_fp, fixed = TRUE)) == 1)
forest_harmonise <- eval(parse(text = sub("k == 2", "k <= 3", code_fp, fixed = TRUE)))

# Vibert (TRANSPHIL, Milan 2024) : absent du CodeBook, année de la présentation
vibert <- read.csv("data/IPD/Vibert_OS_ITT.csv")
vibert$role <- "principal"
vibert$source <- "TRANSPHIL, diapositive OS ITT (Milan 2024)"
vibert$groupe <- factor(vibert$groupe, levels = c("LR", "LT"))
ipd <- bind_rows(ipd, vibert)
annee["Vibert"] <- 2024L

nom <- "OS ITT - ALL"
etudes <- c(analyses[["OS ITT - principal (Dong + Ethun Fig 4A)"]]$etudes, "Vibert")
d <- ipd %>% filter(critere == "OS", analyse == "ITT", etude %in% etudes)
cat("\n=====================================================\n", nom,
    "\n=====================================================\n")
manquantes <- setdiff(etudes, unique(d$etude))
if (length(manquantes)) cat("Absentes des IPD :", paste(manquantes, collapse = ", "), "\n")

# 1. HR par étude
par_etude <- hr_par_etude(d)
print(as.data.frame(par_etude %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))

# 2. One-stage : Cox stratifié par étude (estimation montrée sur le forest plot)
os <- coxph(Surv(time, status) ~ groupe + strata(etude), data = d)
one_stage <- c(exp(coef(os)), exp(confint(os)), summary(os)$coefficients[, "Pr(>|z|)"])
cat("One-stage Cox stratifié : HR =", round(one_stage[1], 2), "[", round(one_stage[2], 2), ";",
    round(one_stage[3], 2), "], p =", signif(one_stage[4], 3), "\n")

# 3. Two-stage REML + HK (pour information) et forest plot
ma <- forest_harmonise(par_etude, one_stage, "results/meta/OS_ITT_ALL.png")
print(ma)

# 4. One-stage avec effet aléatoire (coxme), pour information
d$lt <- as.numeric(d$groupe == "LT")
me <- tryCatch(coxme(Surv(time, status) ~ lt + strata(etude) + (0 + lt | etude), data = d),
               error = function(e) { cat("coxme a échoué :", conditionMessage(e), "\n"); NULL })
if (!is.null(me)) {
  b <- fixef(me); se_b <- sqrt(vcov(me))[1, 1]
  cat("One-stage coxme (effet aléatoire LT vs LR par étude) : HR =", round(exp(b), 2),
      "[", round(exp(b - 1.96 * se_b), 2), ";", round(exp(b + 1.96 * se_b), 2), "], tau2 =",
      round(VarCorr(me)$etude[[1]], 3), "\n")
}
