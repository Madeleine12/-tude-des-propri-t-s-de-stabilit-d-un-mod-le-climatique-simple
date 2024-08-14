from fonctions_pot_2_etats_stables import *
import matplotlib.pyplot as plt

# paramètres
T = 100
nb_T = 10000
min_s = 0.125 
un_peu_moins = 0.25
bcp_moins = 0.5
max_s = 0.85 
nb_val_s_1 = 7 
nb_val_s_2 = 3
al = 0.45

# Calcul des états stables 
x, V_values = values_pot(-2, 2, 30, al)
etats_stables = find_min(V_values, al)

# calcul de tau en fonction de sigma 
tau, sig = tau_values_in_function_of_sigma(min_s, un_peu_moins, bcp_moins, max_s, nb_val_s_1, nb_val_s_2, etats_stables, T, nb_T)
sig_fit, tau_fit, (a, b), var = fit(sig, tau, 'exponential sig')

# affichage des résultats
plt.plot(sig, tau, 'o', label=r"$\alpha$ = "+str(al)+", période d'intégration ="+str(T)+", pas intégration = "+str(T/nb_T))
plt.plot(sig_fit, tau_fit, label = "fit exponentiel, a="+str(a)+",b="+str(b))
plt.legend(loc = 'upper right', fontsize=12)
plt.xlabel(r"$\sigma$", fontsize=16)
plt.ylabel(r"$\tau_{moyen}$", fontsize=16)
plt.show()


# tracé du log de tau en fonction de l'inverse du carré de sigma
sig_inverse = 1/(sig**2)
tau_log = np.log(tau)
plt.plot(sig_inverse, tau_log, 'o', label=r"$\alpha$ = "+str(al)+", période d'intégration ="+str(T)+", pas intégration = "+str(T/nb_T))
sig_inv_fit, tau_log_fit, (a, b), var = fit(sig_inverse, tau_log, "affine")
plt.plot(sig_inv_fit, tau_log_fit, label = "fit affine, a="+str(a)+", b="+str(b))

plt.legend(loc = 'upper right', fontsize=12)
plt.xlabel(r"1/$\sigma^2$", fontsize=16)
plt.ylabel(r"$log(\tau_{moyen})$", fontsize=16)
plt.show()
