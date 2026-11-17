# =====================================================================
# Survie à 1, 3 et 5 ans de Vibert (TRANSPHIL, OS ITT), même procédé que
# survie_1_3_5_ans.R : lecture sur la courbe numérisée, valeur publiée si disponible
# Sortie : results/survie_1_3_5_ans_Vibert.csv (n'écrase aucun fichier existant)
# =====================================================================

# Préparation de survie_1_3_5_ans.R (TEMPS, d, lire_courbe), sans lire recap_IPD.csv
expr <- parse("script/survie_1_3_5_ans.R", keep.source = FALSE)
fin_prep <- which(sapply(expr, function(e) identical(e[[1]], as.name("<-")) &&
                           identical(e[[2]], as.name("recap")))) - 1
for (e in expr[seq_len(fin_prep)]) eval(e)

# k = 1 : la reconstruction retrouve exactement le tableau des patients à risque (pas de décalage)
# Publié (diapositive, encadré KM Est) : 5 ans = 30 % dans les deux bras ; rien à 1 et 3 ans
courbes_v <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
etude   critere analyse bras fichier              k
Vibert  OS      ITT     LT   Vibert_LT_OS_ITT.csv 1
Vibert  OS      ITT     LR   Vibert_LR_OS_ITT.csv 1
")
pub_v <- c(NA, NA, 0.30)

res <- do.call(rbind, lapply(seq_len(nrow(courbes_v)), function(i) {
  cr <- courbes_v[i, ]
  fig <- round(lire_courbe(d(cr$etude, cr$fichier), cr$k, TEMPS), 3)
  pub <- pub_v
  retenu <- ifelse(is.na(pub), fig, pub)
  source <- paste(ifelse(is.na(pub), ifelse(is.na(fig), "NA", "figure"), "publié"), collapse = "/")
  data.frame(cr[, c("etude", "critere", "analyse", "bras")], k = cr$k,
             S12_fig = fig[1], S36_fig = fig[2], S60_fig = fig[3],
             S12_pub = pub[1], S36_pub = pub[2], S60_pub = pub[3],
             S12 = retenu[1], S36 = retenu[2], S60 = retenu[3], source = source)
}))
write.csv(res, "results/survie_1_3_5_ans_Vibert.csv", row.names = FALSE)
print(res, row.names = FALSE)
