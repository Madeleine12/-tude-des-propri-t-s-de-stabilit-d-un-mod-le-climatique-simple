import numpy as np
import sdeint as sd
import matplotlib.pyplot as plt
import itertools as it
from poub.lissage import lissage
from fonctions_pot_2_etats_stables import fit

# Paramètres 
mu = 7.5; theta = 1.; sig = 0.08; sig_0 = 0.9 # Stommel params
ep_a = 0.34 # timescales

# Etats d'équilibres stables 
u0 = np.array([0.57, 0.4]) #ON STATE
u1 =np.array([0.72,0.74]) #OFF state

# Définition des équations
def Stommel(u, t):
    T, S = u[0], u[1]
    dT = - (T - theta)/ep_a - T - mu * np.abs(T - S) * T
    dS = sig_0 - S - mu * np.abs(T - S) * S
    return np.array([dT,dS])

def count_pour_stommel(psivalues, et_1 =0.16999999999999993, et_2=-0.020000000000000018, ep = 0.03):
    """
    Compte le temps moyen passé dans le premier état d'équilibre une fois que le système y est.

    Paramètres : 
    psivalues : Valeurs de psi sur lesquelles on souhaite compter le nombre de transitions
    et_1 : premier état stable
    et_2 : second état stable
    ep : epaisseur du cache 

    Retours : 
    float : temps moyen passé dans le premier état stable
    """
    psi_lisse = lissage(psivalues, 80)
    mask = np.logical_or(psi_lisse < et_1+ep, psi_lisse > et_2-ep)
    psi_filtré = psi_lisse[mask]>0.09
    temps_passe = list(psi_filtré).count(True)*(T/nb_T)
    nb_ds_et_1 = len([key for key, gp in it.groupby(psi_filtré) if key == True])
    # print(temps_passe)
    # print(nb_ds_et_1)
    # return temps_passe/nb_ds_et_1
    return temps_passe/nb_ds_et_1

# Paramétrisation du temps d'intégration 
T = 1000; nb_T = 100000
tvalues = np.linspace(0, T, nb_T)

# Calcul et affichage d'une trajectoire pour psi
psivalues = np.array([T-S for [T, S] in sd.itoint(Stommel, lambda u, t: np.array([[sig, 0],[0, sig]]), u0, tvalues, np.random)])
plt.plot(tvalues, psivalues, label = r"$\psi$ en fonction du temps avec $\mu$ = "+str(mu)+ r", $\theta$ = "+str(theta)+r", amplitude bruit stochastique $\sigma$ = "+str(sig)+", bruit constant "+r'$\sigma_0 = $'+str(sig_0)+r", $\epsilon_a$ = "+str(ep_a))
plt.axhline(0.16999999999999993, c='black', linewidth=0.5) 
plt.axhline(-0.020000000000000018, c='black', linewidth=0.5)  
plt.xlabel('t', fontsize=16)
plt.ylabel(r'$\psi$', fontsize=16)
plt.legend(loc = 'upper right', fontsize=12)
plt.show()

# Calcul et affichage de tau en fonction de sigma
s_min = 0.03
s_max = 0.045
s_values = np.concatenate((np.linspace(s_min, 0.035, 8, endpoint=False), np.linspace(0.035, s_max, 5)))
tau_values = np.array([count_pour_stommel(np.array([T-S for [T, S] in sd.itoint(Stommel, lambda u, t: np.array([[s, 0],[0, s]]), u0, tvalues, np.random)])) for s in s_values])
plt.plot(s_values, tau_values, 'o', label = "période d'intégration ="+str(T)+", pas d'intégration = "+str(T/nb_T)+r", $\mu$ = "+str(mu)+ r", $\theta$ = "+str(theta)+", bruit constant "+r' $\sigma_0 = $'+str(sig_0)+r", $\epsilon_a$ = "+str(ep_a))

s_fit, tau_fit, (a, b), var = fit(s_values, np.array(tau_values), "exponential sig")
plt.plot(s_fit, tau_fit, label = "fit exponentiel, a = "+str(a)+", b = "+str(b)) 

plt.xlabel(r'$\sigma$', fontsize=16)
plt.ylabel(r'$\tau$', fontsize=16)
plt.legend(loc = 'upper right', fontsize=13)
plt.show()


# Calcul et affichage de S en fonction de T
result_int = sd.itoint(Stommel, lambda u, t: np.array([[0, 0],[0, sig]]), u0, tvalues, np.random)
T_values = np.array(result_int[:, 0])
S_values = np.array(result_int[:, 1])
plt.plot(T_values, S_values, label = "période d'intégration ="+str(T)+", pas d'intégration = "+str(T/nb_T)+r", $\mu$ = "+str(mu)+ r", $\theta$ = "+str(theta)+r', $\sigma_0 = $'+str(sig_0)+r', $\sigma = $'+str(sig)+r", $\epsilon_a$ = "+str(ep_a))
plt.xlabel(r'$\Delta T$', fontsize=16)
plt.ylabel(r'$\Delta S$', fontsize=16)
plt.legend(loc = 'upper right', fontsize=13)
plt.show()
