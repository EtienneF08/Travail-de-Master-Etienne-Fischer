# =====================================================================
# Survie à 1, 3 et 5 ans par étude et par bras, lue sur les courbes KM
# publiées (points numérisés automeris), avec correction de l'échelle
# temporelle quand la numérisation est décalée par rapport à la figure
# Entrées : data/<Etude>/*.csv, results/recap_IPD.csv (valeurs publiées)
# Sortie  : results/survie_1_3_5_ans.csv (colonnes S12/S36/S60 = valeurs retenues)
# R de base uniquement
# =====================================================================

TEMPS <- c(12, 36, 60)   # mois

# ---------------------------------------------------------------------
# Facteur k par courbe : t_figure = t_csv / k
# k estimé en superposant les points numérisés sur la figure originale
# (graduations repérées au pixel) et en cherchant le k qui maximise le
# recouvrement avec le tracé. k = 1 : numérisation alignée sur la figure.
#   Breuer RFS      k = 1.050 (axe étiré : marche LT à 62.5 m au lieu de ~59.5 m)
#   Breuer OS       k = 0.965 / 0.960 (axe comprimé)
#   Giuseppucci OS  k = 0.955 (axe comprimé)
# Toutes les autres courbes : recouvrement maximal à k = 1
# ---------------------------------------------------------------------
d <- function(etude, f) file.path("data", etude, f)
courbes <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
etude       critere analyse bras      fichier                     k
Breuer      OS      AT      LT        Breuer_LT_OS_AT.csv         0.965
Breuer      OS      AT      LR        Breuer_LR_OS_AT.csv         0.960
Breuer      RFS     AT      LT        Breuer_LT_RFS_AT.csv        1.050
Breuer      RFS     AT      LR        Breuer_LR_RFS_AT.csv        1.050
Croome      OS      AT      LT        Croome_LT_OS_AT.csv         1
Croome      OS      AT      LR        Croome_LR_OS_AT.csv         1
Croome      OS      ITT     LT        Croome_LT_OS_ITT.csv        1
Croome      OS      ITT     LR        Croome_LR_OS_ITT.csv        1
Croome      RFS     AT      LT        Croome_LT_RFS_AT.csv        1
Croome      RFS     AT      LR        Croome_LR_RFS_AT.csv        1
Dong        OS      AT      LT        Dong_LT_OS_AT.csv           1
Dong        OS      AT      LR_sansVR Dong_LRw_oVR_OS_AT.csv      1
Dong        OS      AT      LR_avecVR Dong_LR+VR_OS_AT.csv        1
Dong        OS      ITT     LT        Dong_LT_OS_ITT.csv          1
Dong        OS      ITT     LR        Dong_LR_OS_ITT.csv          1
Dong        RFS     AT      LT        Dong_LT_RFS_AT.csv          1
Dong        RFS     AT      LR_sansVR Dong_LRw_oVR_RFS_AT.csv     1
Dong        RFS     AT      LR_avecVR Dong_LR+VR_RFS_AT.csv       1
Ethun       OS      AT      LT        Ethun_LT_OS_AT.csv          1
Ethun       OS      AT      LR        Ethun_LR_OS_AT.csv          1
Ethun       OS      ITT     LT        Ethun_LT_OS_ITT.csv         1
Ethun       OS      ITT     LR        Ethun_LR_OS_ITT.csv         1
Giuseppucci OS      AT      LT        Giuseppuci_LT_OS_AT.csv     0.955
Giuseppucci OS      AT      LR        Giuseppuci_LR_OS_AT.csv     0.955
Giuseppucci RFS     AT      LT        Giuseppuci_LT_RFS_AT.csv    1
Giuseppucci RFS     AT      LR        Giuseppuci_LR_RFS_AT.csv    1
Hidalgo     OS      AT      LT        Hidalgo_LT_OS_AT.csv        1
Hidalgo     OS      AT      LR        Hidalgo_LR_OS_AT.csv        1
Iwatsuki    OS      AT      LT        Iwatsuki_LT_OS_AT.csv       1
Iwatsuki    OS      AT      LR        Iwatsuki_LR_OS_AT.csv       1
Pichlmayr   OS      AT      LT        Pichlmayr_LT_OS_AT.csv      1
Pichlmayr   OS      AT      LR        Pichlmayr_LR_OS_AT.csv      1
Poletto     OS      AT      LT        Poletto_LT_OS_AT.csv        1
Poletto     OS      AT      LR        Poletto_LR_OS_AT.csv        1
Poletto     RFS     AT      LT        Poletto_LT_RFS_AT.csv       1
Poletto     RFS     AT      LR        Poletto_LR_RFS_AT.csv       1
Robles      OS      AT      LT        Robles_LT_OS_AT.csv         1
Robles      OS      AT      LR        Robles_LR_OS_AT.csv         1
")

# Survie lue sur la courbe numérisée au temps t de la figure (NA au-delà du suivi)
lire_courbe <- function(csv, k, t) {
  pts <- read.table(csv, sep = ";", dec = ",", strip.white = TRUE, col.names = c("time", "surv"))
  pts <- pts[order(pts$time), ]
  pts$surv <- cummin(pmin(pmax(pts$surv, 0), 1))
  tc <- t * k
  sapply(tc, function(x) {
    if (x > max(pts$time) + 0.5) return(NA)        # tolérance : courbes arrêtées à 59.9 m
    av <- pts$surv[pts$time <= x]
    if (length(av)) tail(av, 1) else 1
  })
}

recap <- read.csv("results/recap_IPD.csv", stringsAsFactors = FALSE)
pub_de <- function(e, c, a, b) {
  r <- recap[recap$etude == e & recap$critere == c & recap$analyse == a & recap$bras == b, ]
  if (!nrow(r)) return(rep(NA, 3))
  as.numeric(r[1, c("S12_pub", "S36_pub", "S60_pub")])
}
rec_de <- function(e, c, a, b) {
  r <- recap[recap$etude == e & recap$critere == c & recap$analyse == a & recap$bras == b, ]
  round(as.numeric(r[1, c("S12_rec", "S36_rec", "S60_rec")]), 3)
}
# Comparaisons dont les valeurs retenues viennent des IPD reconstruits (choix de l'auteur)
reconstruit <- function(etude, analyse) etude == "Dong" || (etude == "Croome" && analyse == "ITT")

lignes <- lapply(seq_len(nrow(courbes)), function(i) {
  cr <- courbes[i, ]
  fig <- round(lire_courbe(d(cr$etude, cr$fichier), cr$k, TEMPS), 3)
  pub <- pub_de(cr$etude, cr$critere, cr$analyse, cr$bras)
  # Valeur retenue : publiée si disponible, sinon lue sur la figure. On privilégie la cohorte totale.
  # Exceptions : Robles (courbe KM retenue, chiffres du texte incohérents) et Giuseppucci
  # (rien de publié : courbe corrigée k = 0.955 plutôt que les lectures manuelles du CodeBook) ;
  # Dong (tout) et Croome ITT : IPD reconstruits (Croome AT : chiffres de l'article).
  # Ethun OS ITT LR : valeurs du CodeBook = cohorte LR complète (n=234), retenues ;
  # la courbe appariée (n=58) ne sert qu'à la reconstruction des IPD
  if (cr$etude %in% c("Robles", "Giuseppucci")) {
    retenu <- fig; source <- "figure (prioritaire)"
  } else if (reconstruit(cr$etude, cr$analyse)) {
    retenu <- rec_de(cr$etude, cr$critere, cr$analyse, cr$bras); source <- "IPD reconstruits"
  } else {
    retenu <- ifelse(is.na(pub), fig, pub)
    source <- paste(ifelse(is.na(pub), ifelse(is.na(fig), "NA", "figure"), "publié"), collapse = "/")
  }
  data.frame(cr[, c("etude", "critere", "analyse", "bras")], k = cr$k,
             S12_fig = fig[1], S36_fig = fig[2], S60_fig = fig[3],
             S12_pub = pub[1], S36_pub = pub[2], S60_pub = pub[3],
             S12 = retenu[1], S36 = retenu[2], S60 = retenu[3], source = source)
})
res <- do.call(rbind, lignes)

# Dong LR total (sans VR + avec VR) : pas de courbe publiée -> IPD reconstruits regroupés
for (c in c("OS", "RFS")) {
  s <- rec_de("Dong", c, "AT", "LR (total)")
  res <- rbind(res, data.frame(etude = "Dong", critere = c, analyse = "AT", bras = "LR (total)", k = NA,
                               S12_fig = NA, S36_fig = NA, S60_fig = NA,
                               S12_pub = NA, S36_pub = NA, S60_pub = NA,
                               S12 = s[1], S36 = s[2], S60 = s[3], source = "IPD reconstruits regroupés"))
}

res <- res[order(res$etude, res$critere, res$analyse, res$bras), ]
write.csv(res, "results/survie_1_3_5_ans.csv", row.names = FALSE)
print(res, row.names = FALSE)
