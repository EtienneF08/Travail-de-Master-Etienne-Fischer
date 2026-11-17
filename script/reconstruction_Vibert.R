# =====================================================================
# Reconstruction des IPD de Vibert (essai TRANSPHIL, World Consensus Meeting
# Milan 2024, Lectures/2024_Vibert_ConsensusConf.pdf) : OS ITT, LT vs LR
# Même procédé que reconstruction_toutes.R (reconstruire, valider, resumer)
# Sorties : data/IPD/Vibert_OS_ITT.csv, results/KM_IPD/Vibert_OS_ITT.png,
#           results/recap_IPD_Vibert.csv (n'écrase aucun fichier existant)
# =====================================================================

# Préparation de reconstruction_toutes.R (fonctions, seuils), sans lancer sa boucle
expr <- parse("script/reconstruction_toutes.R", keep.source = FALSE)
fin_prep <- which(sapply(expr, function(e) identical(e, quote(recap <- list())))) - 1
for (e in expr[seq_len(fin_prep)]) eval(e)

# Vibert absent du CodeBook : valeurs publiées prises sur la diapositive (figure OS ITT)
# N = 20 par bras ; survie à 5 ans = 30 % dans les deux bras (encadré KM Est) ;
# 1 et 3 ans lus sur la courbe (aucune censure avant ~57 mois)
pub_vibert <- list(LT = c(0.80, 0.30, 0.30), LR = c(0.80, 0.45, 0.30))
n_vibert   <- c(LT = 20, LR = 20)

cmp <- C("Vibert", "OS", "ITT", source = "TRANSPHIL, diapositive OS ITT (Milan 2024)",
         B("LT", d("Vibert", "Vibert_LT_OS_ITT.csv"), d("Vibert", "Vibert_LT_OS_ITT.xlsx")),
         B("LR", d("Vibert", "Vibert_LR_OS_ITT.csv"), d("Vibert", "Vibert_LR_OS_ITT.xlsx")))

nom <- paste(cmp$etude, cmp$critere, cmp$analyse, sep = "_")
cat("\n=====================", nom, "=====================\n")

bras <- lapply(cmp$bras, function(b)
  reconstruire(b$csv, b$atrisk, cmp$etude, cmp$critere, cmp$analyse, b$bras))
ipd <- do.call(rbind, bras)
write.csv(ipd, file.path("data/IPD", paste0(nom, ".csv")), row.names = FALSE)

cx <- valider(ipd, gsub("_", " ", nom), file.path("results/KM_IPD", paste0(nom, ".png")))
hr <- sprintf("%.2f (%.2f-%.2f)", exp(coef(cx)), exp(confint(cx))[1], exp(confint(cx))[2])

r <- do.call(rbind, lapply(bras, function(b) {
  g <- b$groupe[1]
  resumer(b, cmp$etude, cmp$critere, cmp$analyse, b$bras[1], g,
          pub = pub_vibert[[g]], n_cb = n_vibert[[g]], atrisk = attr(b, "atrisk"))
}))
r$HR_LT_vs_LR <- hr
r$HR_publie <- ""
r$role <- cmp$role
r$source <- cmp$source
write.csv(r, "results/recap_IPD_Vibert.csv", row.names = FALSE)

cat("\n##### Récapitulatif\n")
print(r[, c("etude", "critere", "analyse", "bras", "N", "evenements", "mediane",
            "S12_rec", "S12_pub", "S36_rec", "S36_pub", "S60_rec", "S60_pub",
            "HR_LT_vs_LR", "alerte")], digits = 3)
