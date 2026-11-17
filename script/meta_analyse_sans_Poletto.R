# =====================================================================
# Forest plots OS AT et RFS AT principaux, sans l'étude Poletto
# Même procédé que meta_analyse.R (HR par étude, one-stage, forest_harmonise),
# avec le one-stage Cox stratifié comme estimation poolée pour k <= 3
# Sortie : results/meta/OS_AT_principal_2.png et RFS_AT_principal_2.png
# (n'écrase aucun fichier existant)
# =====================================================================

# Préparation de meta_analyse.R (données, années, hr_par_etude, forest_harmonise),
# sans lancer sa boucle d'analyses ni réécrire ses tableaux
expr <- parse("script/meta_analyse.R", keep.source = FALSE)
fin_prep <- which(sapply(expr, function(e) identical(e, quote(resume <- list())))) - 1
for (e in expr[seq_len(fin_prep)]) eval(e)

# Même forest plot, mais le one-stage (Cox stratifié) est aussi utilisé pour k = 3 :
# avec 3 études, HK n'a que 2 ddl (IC très larges), même logique que pour k = 2 en ITT
code_fp <- deparse(forest_harmonise)
stopifnot(sum(grepl("k == 2", code_fp, fixed = TRUE)) == 1)
forest_harmonise <- eval(parse(text = sub("k == 2", "k <= 3", code_fp, fixed = TRUE)))

analyses_2 <- list(
  "OS AT - principal 2"  = list(critere = "OS",  analyse = "AT", etudes = setdiff(AT_tous, "Poletto")),
  "RFS AT - principal 2" = list(critere = "RFS", analyse = "AT", etudes = c("Dong", "Breuer", "Giuseppucci"))
)

for (nom in names(analyses_2)) {
  a <- analyses_2[[nom]]
  d <- ipd %>% filter(critere == a$critere, analyse == a$analyse, etude %in% a$etudes)
  cat("\n=====================================================\n", nom,
      "\n=====================================================\n")
  manquantes <- setdiff(a$etudes, unique(d$etude))
  if (length(manquantes)) cat("Absentes des IPD :", paste(manquantes, collapse = ", "), "\n")

  par_etude <- hr_par_etude(d)
  print(as.data.frame(par_etude %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))

  os <- coxph(Surv(time, status) ~ groupe + strata(etude), data = d)
  one_stage <- c(exp(coef(os)), exp(confint(os)), summary(os)$coefficients[, "Pr(>|z|)"])
  cat("One-stage Cox stratifié : HR =", round(one_stage[1], 2), "[", round(one_stage[2], 2), ";",
      round(one_stage[3], 2), "], p =", signif(one_stage[4], 3), "\n")

  ma <- forest_harmonise(par_etude, one_stage,
                         file.path("results/meta", paste0(gsub("[^A-Za-z0-9]+", "_", nom), ".png")))
  print(ma)
}
