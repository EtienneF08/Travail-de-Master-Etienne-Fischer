# =====================================================================
# Reconstruction des IPD (Guyot, IPDfromKM), une courbe (un bras) à la fois
# Pour chaque bras : le CSV automeris (temps ; survie 0-1, décimale ",")
#                    + le fichier _atrisk.xlsx (colonnes time, nrisk)
# =====================================================================

# install.packages(c("IPDfromKM", "readxl", "survival", "survminer"))
source("script/fonctions_IPD.R")               # reconstruire() et valider()
dir.create("data/IPD", showWarnings = FALSE)   # chemins relatifs à la racine du projet (.Rproj)

# Toutes les études d'un coup : script/reconstruction_toutes.R

# =====================================================================
# EXEMPLE : Dong, OS en ITT
# =====================================================================
lt <- reconstruire("data/Dong/Dong_LT_OS_ITT.csv", "data/Dong/Dong_OS_ITT_LT_atrisk.xlsx",
                   "Dong", "OS", "ITT", "LT")
lr <- reconstruire("data/Dong/Dong_LR_OS_ITT.csv", "data/Dong/Dong_OS_ITT_LR_atrisk.xlsx",
                   "Dong", "OS", "ITT", "LR")
dong_os_itt <- rbind(lt, lr)
valider(dong_os_itt, "Dong OS ITT")
write.csv(dong_os_itt, "data/IPD/Dong_OS_ITT.csv", row.names = FALSE)

# Dong AT : le bras LR = 2 courbes (sans et avec résection vasculaire) -> on les empile
# (la colonne groupe regroupe tout ce qui ne commence pas par "LT" en LR)
# lt  <- reconstruire("data/Dong/Dong_LT_OS_AT.csv",      "data/Dong/Dong_OS_AT_LT_atrisk.xlsx",       "Dong", "OS", "AT", "LT")
# lr1 <- reconstruire("data/Dong/Dong_LRw_oVR_OS_AT.csv", "data/Dong/Dong_OS_AT_LRsansVR_atrisk.xlsx", "Dong", "OS", "AT", "LR_sansVR")
# lr2 <- reconstruire("data/Dong/Dong_LR+VR_OS_AT.csv",   "data/Dong/Dong_OS_AT_LRavecVR_atrisk.xlsx", "Dong", "OS", "AT", "LR_avecVR")
# dong_os_at <- rbind(lt, lr1, lr2); valider(dong_os_at, "Dong OS AT")
# write.csv(dong_os_at, "data/IPD/Dong_OS_AT.csv", row.names = FALSE)

# =====================================================================
# À LA FIN : tout rassembler en un seul fichier pour la méta-analyse
# =====================================================================
# ipd_tous <- do.call(rbind, lapply(list.files("data/IPD", full.names = TRUE), read.csv))
# write.csv(ipd_tous, "data/IPD_tous.csv", row.names = FALSE)
