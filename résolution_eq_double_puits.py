from fonctions_pot_2_etats_stables import *
import matplotlib.pyplot as plt

# Paramètres
al = 0.45
sigma = 0.20
T = 50
nb_T = 5000

#potentiel
x, V_values = values_pot(-2, 2, 30, al)
plt.plot(x, V_values, label = "potentiel en fonction de l'état du système")
plt.legend(loc = 'upper right', fontsize=12)
plt.title("V en fonction de X")
plt.xlabel('X', fontsize=16)
plt.ylabel('Potentiel', fontsize=16)
plt.show()

#etats stables
etats = find_min(V_values, al)

#résolution équation avec euler-maruyama
tvalues, xvalues = euler_pour_pot(etats, T, nb_T, sigma, al)

#affichage de la trajectoire obtenue 
plt.plot(tvalues, xvalues, label = "sigma =" + str(sigma)+", alpha = "+str(al)+", période d'intégration = "+str(T)+", pas d'intégration = "+str(T/nb_T))
plt.axvline(0, c='black', linewidth=0.5) 
plt.axhline(0, c='black', linewidth=0.5)  
plt.legend(loc = 'upper right', fontsize=12)
plt.title("X en fonction du temps")
plt.show()

# hitogramme sur les valeurs obtenues 
plt.hist(xvalues, bins=700, alpha=0.7, color='green', edgecolor='black', label = "sigma =" + str(sigma)+", alpha = "+str(al)+", période d'intégration = "+str(T)+", pas d'intégration = "+str(T/nb_T)) 
plt.legend(loc = 'upper right', fontsize=12)
plt.title("Densité de probabilité")
plt.xlabel('X', fontsize=16)
plt.ylabel('Fréquence', fontsize=16)
plt.show()