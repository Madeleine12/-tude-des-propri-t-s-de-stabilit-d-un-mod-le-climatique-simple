from fonctions_pot_2_etats_stables import *
import matplotlib.pyplot as plt

# Paramètres
min_a = 2*np.sqrt(0.2)
max_a = 2*np.sqrt(1)
T = 100
nb_T = 10000
sigma = 0.6

# Calcul de tau en fonction de D
D, tau = tau_values_in_function_of_D(min_a, max_a, T, nb_T, sigma)

# Affichage du résultat
plt.plot(D, tau, 'o', label=r"$\sigma = $"+str(sigma)+", période d'intégration = "+str(T)+", pas intégration = "+str(T/nb_T)+r", $min_a = $"+str(min_a)+r", $\max_a = $"+str(max_a))
D_fit, tau_fit, (a, b), var = fit(D, tau, "exponential D")
plt.plot(D_fit, tau_fit, label="fit, a="+str(a)+", b="+str(b))
plt.grid()
plt.legend(loc = 'upper right', fontsize=13)
plt.xlabel(r"$\Delta V$", fontsize=16)
plt.ylabel(r"$\tau_{moyen}$", fontsize=16)
plt.show()

# log de tau en fonction de D
tau_log = np.log(tau)
D_fit, tau_log_fit, (a, b), variance_log = fit(D, tau_log, "affine")
print("matrice de covariance fit affine")
print(variance_log)

# Affichage du résultat
plt.plot(D, tau_log, 'o', label=r"simulation numérique avec $\sigma = $"+str(sigma)+r", $T=$"+str(T)+", pas d'intégration = "+str(T/nb_T))
plt.plot(D_fit, tau_log_fit, label = "fit fonction affine, a="+str(a)+", b = "+str(b))
plt.legend(loc = 'upper left', fontsize=13)
plt.xlabel(r"$\Delta V$", fontsize=16)
plt.ylabel(r"$log(\tau_{moyen})$", fontsize=16)
plt.show()

