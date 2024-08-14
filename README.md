# Étude des propriétés de stabilité d'un modèle climatique simple
## Résolution numérique du double puit de potentiel
### 1
Dans le fichier "fonctions_pot_2_etats_stables", on trouve les fonctions utiles pour la résolution du double puits de potentiel et la fonction de fit utilisée également dans la résolution du modèle de Stommel. 
L'objectif, les arguments et les retours de chaque fonctions sont détaillés dans le code de celles-ci.
### 2
Dans le fichier "résolution_eq_double_puits", on résout l'équation différentielle stochastique décrivant le double puits de potentiel puis on affiche une trajectoire en fonction du temps et l'histogramme correspondant à cette trajectoire.
## Résolution numérique du modèle de Stommel
### 1
Dans le fichier "integration_stommel", on défini une fonction pour le compte des transitions adapté à un potentiel asymétrique puis on résout le modèle de Stommel et on affiche une trajectoire de $\psi$ e fonction du temps et $\Delta S$ en fonction de $\Delta T$.
### 2
Dans le fichier "diag_biffurcation_stommel", on résout les équations statiques du modèle de Stommel et on affiche le diagremme de biffurcation correspondant.
### 3
Dans les fichiers "tau_moy_de_delta_V" et "tau_moy_de_sigma", on résout le système d'équation différentielles stochastique du modèle de Stommel pour plusieurs valeurs de $\Delta V$ / $\sigma$ puis on trace les résultats obtenus.
