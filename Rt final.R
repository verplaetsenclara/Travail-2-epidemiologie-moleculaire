if(!require(EpiEstim))install.packages("EpiEstim")
library(EpiEstim)

# 1. Chargement et nettoyage de vos données
tab <- read.csv("H5N1_genomic_analyses_simulated_dataset_4-1.csv", header = TRUE)
tab$date <- as.Date(tab$collection_date, format = "%Y-%m-%d")

# 2. Section 1 manquante / tronquée du PDF : Calcul de l'incidence quotidienne
# Détermination du jour de l'épidémie pour chaque ligne (de 1 à N)
days <- interval(min(ymd(tab$date)), ymd(tab$date)) %/% days(1) + 1

# Calcul du nombre total de jours (correction du '-1' erroné du PDF en '+1')
total_number_of_days <- interval(min(ymd(tab$date)), max(ymd(tab$date))) %/% days(1) + 1

# Remplissage sécurisé avec des 0 
daily_cases <- rep(0, total_number_of_days)

for (i in 1:length(daily_cases)) {
  daily_cases[i] <- sum(days == i)
}

n <- 1000  # Nombre d'itérations
mean_range <- c(3.0, 5.0)  # Intervalle uniforme pour la moyenne du SI (adapté au H5N1 aviaire)
sd_range <- c(1.0, 2.0)    # Intervalle uniforme pour l'écart-type du SI (adapté au H5N1 aviaire)

# Configuration de la fenêtre glissante de 7 jours
t_start <- seq(2, length(daily_cases) - 6) 
t_end <- seq(8, length(daily_cases))

# Matrice pour stocker toutes les estimations de chaque itération
all_Rt <- matrix(NA, nrow = length(t_start), ncol = n)

# Boucle officielle de propagation
for (i in 1:n) {
  # Tirage stochastique des paramètres du SI pour l'itération i
  mean_si_i <- runif(1, mean_range[1], mean_range[2])
  sd_si_i   <- runif(1, sd_range[1], sd_range[2])
  
  # Estimation du R pour l'itération i
  res_i <- estimate_R(incid = daily_cases, 
                      method = "parametric_si", 
                      config = make_config(list(mean_si = mean_si_i, 
                                                std_si = sd_si_i, 
                                                t_start = t_start, 
                                                t_end = t_end)))
  all_Rt[, i] <- res_i$R$`Mean(R)`
}

# 4. Calcul des statistiques globales (Médiane et Intervalles à 95%)
R_median <- apply(all_Rt, 1, median, na.rm = TRUE)
R_lower  <- apply(all_Rt, 1, quantile, probs = 0.025, na.rm = TRUE)
R_upper  <- apply(all_Rt, 1, quantile, probs = 0.975, na.rm = TRUE)

# Axe du temps : calcul du milieu de la fenêtre glissante
Rt_days <- (t_start + t_end) / 2
R_dates <- decimal_date(min(ymd(tab$date)) + Rt_days)

# 5. Dessin du graphique au format "ruban d'incertitude" (Style publication du PDF)
# Configuration des marges et paramètres graphiques de base
par(mar = c(3.5, 3.5, 1.5, 1), lwd = 0.5, bty = "o", col.axis = "gray30", fg = "gray30")

# Création du graphique vide ajusté à vos valeurs
plot(R_dates, R_median, type = "n", xlab = "", ylab = "",
     ylim = c(0, max(R_upper, na.rm = TRUE)), xlim = c(min(R_dates), max(R_dates)), axes = TRUE)

# Ajout du ruban gris (Intervalle de confiance à 95%)
xx_1 <- c(R_dates, rev(R_dates))
yy_1 <- c(R_lower, rev(R_upper))
polygon(xx_1, yy_1, col = rgb(187/255, 187/255, 187/255, 0.4), border = NA)

# Ajout de la ligne de la médiane (ligne rouge ou noire selon préférence)
lines(R_dates, R_median, lwd = 1.2, col = "black")

# Ajout de la ligne pointillée critique au seuil de Rt = 1
abline(h = 1, lty = 2, lwd = 0.8, col = "firebrick")

# Titre de l'axe Y
mtext("Nombre de reproduction effectif (Rt)", side = 2, line = 2.2, col = "gray30", cex = 0.9)
title(main = "Évolution du Rt avec propagation de l'incertitude du SI", line = 0.5, cex.main = 1)
