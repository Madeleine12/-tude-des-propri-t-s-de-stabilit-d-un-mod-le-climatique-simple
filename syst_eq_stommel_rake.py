import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Paramètres
mu = 7.5; theta = 1; ep_a = 0.34

def résolution_système_statique(sig_0, mu = 7.5, theta = 1, ep_a = 0.34, seuil = 10**-10):
    """
    Résout le système d'équations statiques décrivant le système de stommel pour les états stables.

    Paramètres : 
    sig_0 : flux d'eau douce
    mu, theta, ep_a : paramètres du modèle
    seuil : seuil de tolérance pour la partie imaginaire des solutions

    Retours :
    list : solutions obtenues
    """

    # Définir les variables symboliques
    T, S = sp.symbols('T S')

    # si T-S >= 0
    eq1 = - (T - theta)/ep_a - T - mu * (T - S) * T
    eq2 = sig_0 - S - mu * (T - S) * S

    solution_pos = sp.solve((eq1, eq2), (T, S))
    real_sol_pos = [(complex(x).real, complex(y).real) for (x, y) in solution_pos if all(abs(complex(el).imag) <= seuil for el in (x, y))]
    psi_pos = [x-y for (x,y) in real_sol_pos if x-y>=0]

    # si x-y < 0
    eq1 = - (T - theta)/ep_a - T - mu * (S - T) * T
    eq2 = sig_0 - S - mu * (S - T) * S

    solution_neg = sp.solve((eq1, eq2), (T, S))
    real_sol_neg = [(complex(x).real, complex(y).real) for (x,y) in solution_neg if all(abs(complex(el).imag <= seuil) for el in (x, y))]
    psi_neg = [x-y for (x, y) in real_sol_neg if x - y < 0]
    return psi_pos + psi_neg

# intervalle de valeurs pour sigma
sigma = np.linspace(0.6, 1.1, 100)

# Calcul des solutions pour psi
abs_sig = []
psi = []
for s in sigma:
    sols = résolution_système_statique(s, mu, theta, ep_a)
    psi.append(sols)
    # abs_sig.append(s for i in range(len(sols)))
    for i in range(len(sols)):
        abs_sig.append(s) 

ord_psi = [ele for sol in psi for ele in sol]
    
# Affichage du résultat
plt.plot(abs_sig, ord_psi, 'o', label = r"solutions au système d'équations pour $\mu$ =" + str(mu) + r", $\theta$ =" + str(theta) + r", $\epsilon_a$ = "+str(ep_a))
plt.xlabel(r"$\sigma$", fontsize=16)
plt.ylabel(r"$\psi$", fontsize=16)
plt.axhline(0, c='black', linewidth=0.5) 
plt.legend(loc = 'upper right', fontsize=12)
plt.show()
