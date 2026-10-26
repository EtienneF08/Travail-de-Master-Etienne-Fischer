# =====================================================================
# Reconstruction des IPD de toutes les comparaisons LT vs LR
# (une étude à la fois), contrôle contre le CodeBook, puis fusion
# Sorties : data/IPD/Etude_Critere_Analyse.csv, data/IPD_tous.csv,
#           results/recap_IPD.csv, results/KM_IPD/*.png
# =====================================================================

source("script/fonctions_IPD.R")
dir.create("data/IPD", showWarnings = FALSE)
dir.create("results/KM_IPD", showWarnings = FALSE, recursive = TRUE)

SEUIL <- 0.03                  # écart toléré entre survie reconstruite et publiée
TEMPS <- c(12, 36, 60)         # mois

# ---------------------------------------------------------------------
# Valeurs publiées : data/CodeBook.1.def..xlsx (feuille Extraction)
# ---------------------------------------------------------------------
cb_brut <- as.data.frame(read_excel("data/CodeBook.1.def..xlsx", "Extraction",
                                    col_names = FALSE, col_types = "text"))
cb <- cb_brut[3:nrow(cb_brut), ]
names(cb) <- ifelse(is.na(cb_brut[1, ]), paste0("col", seq_along(cb_brut)), unlist(cb_brut[1, ]))
nom_cb <- c(Giuseppucci = "A.Giuseppuci", Poletto = "E.Poletto", Hidalgo = "E.Hidalgo",
            Ethun = "Ethnun CG.", Dong = "Y.Dong", Breuer = "E.Breuer", Croome = "P.Croome",
            Iwatsuki = "S.Iwatsuki", Robles = "R.Robles", Pichlmayr = "R.Pichlmayr")
num <- function(x) suppressWarnings(as.numeric(x))

cb_val <- function(etude, variable) {
  ligne <- cb[cb$Study == nom_cb[[etude]], ]
  if (!variable %in% names(ligne)) return(NA)
  num(ligne[[variable]])
}

# Survie publiée à 1/3/5 ans. Giuseppucci : colonnes AT vides, valeurs lues sur les courbes
# KM et saisies à droite de "Comments" (col86-91 = LR OS 1/3/5, RFS 1/3/5 ; col93-98 = LT)
survie_publiee <- function(etude, critere, analyse, groupe) {
  if (etude == "Giuseppucci") {
    debut <- if (groupe == "LR") 86 else 93
    if (critere == "RFS") debut <- debut + 3
    return(num(unlist(cb[cb$Study == nom_cb[[etude]], debut:(debut + 2)])))
  }
  sapply(c(1, 3, 5), function(an)
    cb_val(etude, sprintf("%s_%dyr_perc_%s_%s", tolower(critere), an, groupe, analyse)))
}

n_codebook <- function(etude, analyse, groupe) {
  if (analyse == "ITT") cb_val(etude, if (groupe == "LR") "n_LR_ITT" else "n_LT_ITT")
  else cb_val(etude, if (groupe == "LR") "n_resected" else "n_transplated")
}

hr_publie <- function(etude, critere) {
  v <- cb[cb$Study == nom_cb[[etude]], paste0("hr_", tolower(critere))]
  if (is.na(v) || v %in% c("NA", "NO")) "" else v
}

# ---------------------------------------------------------------------
# Liste des comparaisons, dans l'ordre demandé
# ---------------------------------------------------------------------
B <- function(bras, csv, atrisk) list(bras = bras, csv = csv, atrisk = atrisk)
C <- function(etude, critere, analyse, ..., role = "principal", source = "")
  list(etude = etude, critere = critere, analyse = analyse, bras = list(...), role = role, source = source)

d <- function(etude, fichier) file.path("data", etude, fichier)

comparaisons <- list(
  # Ethun ITT = Fig 4A : LR <3 cm N0 (n = 58) vs LT tous inscrits (n = 66). Analyse principale ;
  # les survies ITT du CodeBook portent sur toute la cohorte (LR n = 234), donc non comparables.
  # Sensibilité : méta-analyse sans Ethun (data/IPD_tous_sansEthun.csv)
  C("Ethun", "OS", "ITT", source = "Fig 4A (LR <3 cm N0 vs LT tous inscrits)",
    B("LT", d("Ethun", "Ethun_LT_OS_ITT.csv"), d("Ethun", "Ethun_OS_ITT_LT_atrisk.xlsx")),
    B("LR", d("Ethun", "Ethun_LR_OS_ITT.csv"), d("Ethun", "Ethun_OS_ITT_LR_atrisk.xlsx"))),
  C("Dong", "OS", "AT",
    B("LT",        d("Dong", "Dong_LT_OS_AT.csv"),      d("Dong", "Dong_OS_AT_LT_atrisk.xlsx")),
    B("LR_sansVR", d("Dong", "Dong_LRw_oVR_OS_AT.csv"), d("Dong", "Dong_OS_AT_LRsansVR_atrisk.xlsx")),
    B("LR_avecVR", d("Dong", "Dong_LR+VR_OS_AT.csv"),   d("Dong", "Dong_OS_AT_LRavecVR_atrisk.xlsx"))),
  C("Dong", "OS", "ITT",
    B("LT", d("Dong", "Dong_LT_OS_ITT.csv"), d("Dong", "Dong_OS_ITT_LT_atrisk.xlsx")),
    B("LR", d("Dong", "Dong_LR_OS_ITT.csv"), d("Dong", "Dong_OS_ITT_LR_atrisk.xlsx"))),
  C("Ethun", "OS", "AT",
    B("LT", d("Ethun", "Ethun_LT_OS_AT.csv"), d("Ethun", "Ethun_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Ethun", "Ethun_LR_OS_AT.csv"), d("Ethun", "Ethun_OS_AT_LR_atrisk.xlsx"))),
  C("Poletto", "OS", "AT",
    B("LT", d("Poletto", "Poletto_LT_OS_AT.csv"), d("Poletto", "Poletto_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Poletto", "Poletto_LR_OS_AT.csv"), d("Poletto", "Poletto_OS_AT_LR_atrisk.xlsx"))),
  C("Poletto", "RFS", "AT",
    B("LT", d("Poletto", "Poletto_LT_RFS_AT.csv"), d("Poletto", "Poletto_RFS_AT_LT_atrisk.xlsx")),
    B("LR", d("Poletto", "Poletto_LR_RFS_AT.csv"), d("Poletto", "Poletto_RFS_AT_LR_atrisk.xlsx"))),
  C("Hidalgo", "OS", "AT",
    B("LT", d("Hidalgo", "Hidalgo_LT_OS_AT.csv"), d("Hidalgo", "Hidalgo_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Hidalgo", "Hidalgo_LR_OS_AT.csv"), d("Hidalgo", "Hidalgo_OS_AT_LR_atrisk.xlsx"))),
  C("Breuer", "OS", "AT",
    B("LT", d("Breuer", "Breuer_LT_OS_AT.csv"), d("Breuer", "Breuer_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Breuer", "Breuer_LR_OS_AT.csv"), d("Breuer", "Breuer_OS_AT_LR_atrisk.xlsx"))),
  C("Breuer", "RFS", "AT",
    B("LT", d("Breuer", "Breuer_LT_RFS_AT.csv"), d("Breuer", "Breuer_RFS_AT_LT_atrisk copie.xlsx")),
    B("LR", d("Breuer", "Breuer_LR_RFS_AT.csv"), d("Breuer", "Breuer_RFS_AT_LR_atrisk.xlsx"))),
  C("Giuseppucci", "OS", "AT",
    B("LT", d("Giuseppucci", "Giuseppuci_LT_OS_AT.csv"), d("Giuseppucci", "Giuseppucci_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Giuseppucci", "Giuseppuci_LR_OS_AT.csv"), d("Giuseppucci", "Giuseppucci_OS_AT_LR_atrisk.xlsx"))),
  C("Giuseppucci", "RFS", "AT",
    B("LT", d("Giuseppucci", "Giuseppuci_LT_RFS_AT.csv"), d("Giuseppucci", "Giuseppucci_RFS_AT_LT_atrisk.xlsx")),
    B("LR", d("Giuseppucci", "Giuseppuci_LR_RFS_AT.csv"), d("Giuseppucci", "Giuseppucci_RFS_AT_LR_atrisk.xlsx"))),
  C("Iwatsuki", "OS", "AT",
    B("LT", d("Iwatsuki", "Iwatsuki_LT_OS_AT.csv"), d("Iwatsuki", "Iwatsuki_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Iwatsuki", "Iwatsuki_LR_OS_AT.csv"), d("Iwatsuki", "Iwatsuki_OS_AT_LR_atrisk.xlsx"))),
  C("Robles", "OS", "AT",
    B("LT", d("Robles", "Robles_LT_OS_AT.csv"), d("Robles", "Robles_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Robles", "Robles_LR_OS_AT.csv"), d("Robles", "Robles_OS_AT_LR_atrisk.xlsx"))),
  C("Pichlmayr", "OS", "AT",
    B("LT", d("Pichlmayr", "Pichlmayr_LT_OS_AT.csv"), d("Pichlmayr", "Pichlmayr_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Pichlmayr", "Pichlmayr_LR_OS_AT.csv"), d("Pichlmayr", "Pichlmayr_OS_AT_LR_atrisk.xlsx"))),
  C("Dong", "RFS", "AT",
    B("LT",        d("Dong", "Dong_LT_RFS_AT.csv"),      d("Dong", "Dong_RFS_AT_LT_atrisk.xlsx")),
    B("LR_sansVR", d("Dong", "Dong_LRw_oVR_RFS_AT.csv"), d("Dong", "Dong_RFS_AT_LRsansVR_atrisk.xlsx")),
    B("LR_avecVR", d("Dong", "Dong_LR+VR_RFS_AT.csv"),   d("Dong", "Dong_RFS_AT_LRavecVR_atrisk.xlsx"))),
  # Sensibilité
  C("Croome", "OS", "ITT", role = "sensibilite",
    B("LT", d("Croome", "Croome_LT_OS_ITT.csv"), d("Croome", "Croome_OS_ITT_LT_atrisk.xlsx")),
    B("LR", d("Croome", "Croome_LR_OS_ITT.csv"), d("Croome", "Croome_OS_ITT_LR_atrisk.xlsx"))),
  C("Croome", "OS", "AT", role = "sensibilite",
    B("LT", d("Croome", "Croome_LT_OS_AT.csv"), d("Croome", "Croome_OS_AT_LT_atrisk.xlsx")),
    B("LR", d("Croome", "Croome_LR_OS_AT.csv"), d("Croome", "Croome_OS_AT_LR_atrisk.xlsx"))),
  C("Croome", "RFS", "AT", role = "sensibilite",
    B("LT", d("Croome", "Croome_LT_RFS_AT.csv"), d("Croome", "Croome_RFS_AT_LT_atrisk.xlsx")),
    B("LR", d("Croome", "Croome_LR_RFS_AT.csv"), d("Croome", "Croome_RFS_AT_LR_atrisk.xlsx")))
  # Dong apparié (ITTapp) : fichiers atrisk présents mais aucun CSV de courbe -> non reconstruit
)

# ---------------------------------------------------------------------
# Résumé d'un sous-ensemble d'IPD (un bras ou un groupe)
# ---------------------------------------------------------------------
survie_a <- function(ipd, t) {
  if (t > max(ipd$time) + 1) return(NA)   # suivi trop court (tolérance 1 mois, ex. 35.9 pour 36)
  fit <- survfit(Surv(time, status) ~ 1, data = ipd)
  summary(fit, times = t, extend = TRUE)$surv
}

resumer <- function(ipd, etude, critere, analyse, bras, groupe, pub, n_cb, atrisk = NULL) {
  fit <- survfit(Surv(time, status) ~ 1, data = ipd)
  rec <- sapply(TEMPS, function(t) survie_a(ipd, t))
  alertes <- character(0)
  ecart <- rec - pub
  for (k in which(abs(ecart) > SEUIL))
    alertes <- c(alertes, sprintf("S%d : %+.1f pts", TEMPS[k], 100 * ecart[k]))
  if (!is.null(atrisk)) {
    diff_ar <- atrisk[atrisk$reconstruit != atrisk$publie, ]
    if (nrow(diff_ar)) alertes <- c(alertes, paste0("N à risque ≠ tableau (",
      paste(sprintf("%g m : %d vs %d", diff_ar$time, diff_ar$reconstruit, diff_ar$publie), collapse = ", "), ")"))
  }
  if (!is.na(n_cb) && n_cb != nrow(ipd))
    alertes <- c(alertes, sprintf("N courbe %d ≠ CodeBook %d", nrow(ipd), n_cb))
  data.frame(etude = etude, critere = critere, analyse = analyse, bras = bras,
             N = nrow(ipd), N_codebook = n_cb, evenements = sum(ipd$status),
             mediane = unname(summary(fit)$table["median"]),
             S12_rec = rec[1], S12_pub = pub[1], S36_rec = rec[2], S36_pub = pub[2],
             S60_rec = rec[3], S60_pub = pub[3],
             nrisk_tableau = if (is.null(atrisk)) NA else paste(atrisk$publie, collapse = "/"),
             nrisk_reconstruit = if (is.null(atrisk)) NA else paste(atrisk$reconstruit, collapse = "/"),
             alerte = paste(alertes, collapse = " ; "))
}

# ---------------------------------------------------------------------
# Boucle : une étude à la fois
# ---------------------------------------------------------------------
recap <- list(); ipd_tous <- list()
for (cmp in comparaisons) {
  nom <- paste(cmp$etude, cmp$critere, cmp$analyse, sep = "_")
  cat("\n=====================", nom, "=====================\n")

  bras <- lapply(cmp$bras, function(b)
    reconstruire(b$csv, b$atrisk, cmp$etude, cmp$critere, cmp$analyse, b$bras))
  ipd <- do.call(rbind, bras)
  write.csv(ipd, file.path("data/IPD", paste0(nom, ".csv")), row.names = FALSE)
  ipd_tous[[nom]] <- cbind(ipd, role = cmp$role, source = cmp$source)

  cx <- valider(ipd, gsub("_", " ", nom), file.path("results/KM_IPD", paste0(nom, ".png")))
  hr <- sprintf("%.2f (%.2f-%.2f)", exp(coef(cx)), exp(confint(cx))[1], exp(confint(cx))[2])

  lignes <- list()
  # une ligne par courbe numérisée
  for (k in seq_along(bras)) {
    b <- bras[[k]]; g <- b$groupe[1]; plusieurs <- sum(sapply(bras, function(x) x$groupe[1]) == g) > 1
    lignes[[length(lignes) + 1]] <- resumer(
      b, cmp$etude, cmp$critere, cmp$analyse, b$bras[1], g,
      pub  = if (plusieurs) rep(NA, 3) else survie_publiee(cmp$etude, cmp$critere, cmp$analyse, g),
      n_cb = if (plusieurs) NA else n_codebook(cmp$etude, cmp$analyse, g),
      atrisk = attr(b, "atrisk"))
  }
  # bras LR composé de plusieurs courbes (Dong) : ligne LR totale
  if (sum(ipd$groupe == "LR") > 0 && sum(sapply(bras, function(x) x$groupe[1]) == "LR") > 1)
    lignes[[length(lignes) + 1]] <- resumer(
      ipd[ipd$groupe == "LR", ], cmp$etude, cmp$critere, cmp$analyse, "LR (total)", "LR",
      survie_publiee(cmp$etude, cmp$critere, cmp$analyse, "LR"), n_codebook(cmp$etude, cmp$analyse, "LR"))

  r <- do.call(rbind, lignes)
  r$HR_LT_vs_LR <- hr
  r$HR_publie <- hr_publie(cmp$etude, cmp$critere)
  r$role <- cmp$role
  r$source <- cmp$source
  recap[[nom]] <- r
}

recap <- do.call(rbind, recap); rownames(recap) <- NULL
write.csv(recap, "results/recap_IPD.csv", row.names = FALSE)

ipd_tous <- do.call(rbind, ipd_tous); rownames(ipd_tous) <- NULL
write.csv(ipd_tous, "data/IPD_tous.csv", row.names = FALSE)
# Analyse de sensibilité : toutes les études sauf Ethun
write.csv(ipd_tous[ipd_tous$etude != "Ethun", ], "data/IPD_tous_sansEthun.csv", row.names = FALSE)

cat("\n\n##### Récapitulatif\n")
print(recap[, c("etude", "critere", "analyse", "bras", "N", "evenements", "mediane",
                "S12_rec", "S12_pub", "S36_rec", "S36_pub", "S60_rec", "S60_pub",
                "HR_LT_vs_LR", "alerte")], digits = 3)
cat("\nIPD_tous :", nrow(ipd_tous), "patients,", length(unique(paste(ipd_tous$etude, ipd_tous$critere, ipd_tous$analyse))), "comparaisons\n")
