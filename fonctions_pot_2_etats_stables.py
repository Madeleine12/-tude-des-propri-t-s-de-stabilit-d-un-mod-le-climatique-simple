import numpy as np
from scipy.signal import argrelextrema
import sdeint as sd
import itertools as it
from scipy.optimize import curve_fit


# fonction potentiel
def V(x, al = 0.45):
    return x**4-al*x**2

def values_pot(min, max, nb, al = 0.45): 
    """
    Calcul des valeurs du potentiel sur un intervalle 

    Paramètres : 
    min : minimum de l'intervalle
    max : maximum de l'intervalle 
    nb : nombres de valeurs calculées
    al : paramètre pour la profondeur des puits de potentiel

    Retours : 
    array : valeurs pour x réparties uniformément sur l'intervalle 
    list : images des valeurs de x par la fonction potentiel
    """
    x_values = np.linspace(min, max, nb)
    return x_values, [V(a, al) for a in x_values]

def find_min(V_values, al = 0.45):
    """
    Trouve les minimas locaux du potentiel

    Paramètres : 
    V_values : liste contenant les valeurs du potentiel 
    al : paramètre pour la profondeur des puits de potentiel

    Ratours : 
    list : liste des deux états d'équilibres stables 
    """
    x, _ = values_pot(-2, 2, 30, al)
    etats_stables = [x[indices] for indices in argrelextrema(np.array(V_values), np.less)[0]]
    return etats_stables

def euler_pour_pot(etats_stables, T, nb_T, sigma, al = 0.45):
    """
    Résout l'équation différentielle stochastique avec la méthode d'Euler-Maruyama

    Paramètres : 
    etats_stables : liste des états stables du potentiel
    T : temps d'intégration 
    nb_T : nombres de valeurs calculées (réparties uniformément entre 0 et T)
    sigma : écart-type du bruit stochastique ajouté
    al : paramètre pour la profondeur des puits de potentiel

    Retours : 
    array : valeurs du temps utilisé pour la résolution (intervalle de nb_T valeurs allant de 0 à T)
    array : valeurs trouvées pour x lors de la résolution
    """
    tvalues = np.linspace(0, T, nb_T)
    xvalues = sd.itoint(lambda x, t: -4 * x**3 + 2 * al * x, lambda x, t: sigma, etats_stables[0], tvalues, np.random)
    return tvalues, xvalues

def count_nb_transitions(x_values, etats_stables = "list", ep=0.1):
    """"
    Compte le nombre de transtitions entre deux états d'équilibre 

    Paramètres : 
    x_values : données
    etats_stables : liste des états stables du potentiel
    ep : largeur du masque appliqué aux données 

    Retours : 
    int : nombre de transitions entre les deux états stables 
    """
    mask = np.logical_or(x_values < etats_stables[0]+ep, x_values > etats_stables[1]-ep)
    return len([key for key, _ in it.groupby(x_values[mask]>0)])

def tau_values_in_function_of_sigma(min_s, un_peu_moins, bcp_moins, max_s, nb_val_s_1, nb_val_s_2, etats_stables, T, nb_val_T, ep=0.1, al = 0.45):
    """
    Calcul des valeurs de tau pour sigma dans un intervalle de valeurs

    Paramètres : 
    min_s : valeur minimale de l'intervalle de sigma
    un_peu_moins : seuil de début de la fin de la pente 
    bcp_moins : seuil de fin de la pente 
    max_s : valeur maximale de l'intervalle de sigma 
    nb_val_s_1 : nombre de valeurs avant le premier seuil
    nb_val_s_2 : nombre de valeurs entre les deux seuils 
    etats_stables : liste des états stables du potentiel
    T : temps d'intégration 
    nb_T : nombres de valeurs calculées (réparties uniformément entre 0 et T)
    ep : largeur du masque appliqué aux données pour le comptes du nombre de transitions
    al : paramètre pour la profondeur des puits de potentiel

    Retours : 
    array : valeurs trouvées pour le temps moyen passé dans un état d'équilibre stable 
    array : valeurs de sigma utilisées 
    """
    t_values = np.linspace(0, T, nb_val_T)
    s_values = np.concatenate((np.linspace(min_s, un_peu_moins, nb_val_s_1, endpoint=False), np.linspace(un_peu_moins, bcp_moins, nb_val_s_2, endpoint=False), np.linspace(bcp_moins, max_s, 3)))
    tau = [T/count_nb_transitions(sd.itoint(lambda x, t: -4 * x**3 + 2 * al * x, lambda x, t : s, etats_stables[0], t_values, np.random), etats_stables, ep) for s in s_values]
    return np.array(tau), np.array(s_values)

def tau_values_in_function_of_D(min_a, max_a, T, nb_T, sigma, ep = 0.1, al = 0.45):
    """
    Calcul des valeurs de tau pour alpha dans un intervalle de valeurs

    Paramètres : 
    min_a : valeur minimale de l'intervalle de alpha 
    max_a : valeur maximale de l'intervalle de alpha 
    T : temps d'intégration 
    nb_T : nombres de valeurs calculées (réparties uniformément entre 0 et T)
    sigma : écart type du bruit 
    ep : largeur du masque appliqué aux données pour le comptes du nombre de transitions
    al : paramètre pour la profondeur des puits de potentiel

    Retours : 
    array : valeurs de deltaV correspondant aux valeurs de alpha utilisées 
    array : valeurs trouvées pour le temps moyen passé dans un état d'équilibre stable 
    """
    seuil = 2*np.sqrt(0.7)
    al_values = np.concatenate((np.linspace(min_a, seuil, 4, endpoint=False), np.linspace(seuil, max_a, 6)))
    tau = []
    D = []
    for al in al_values:
        etats_stables = find_min(values_pot(-2,2, 30, al)[1], al)
        if len(etats_stables) < 2:
            print("Le potentiel n'a pas deux états stbles")
        tau.append(T/count_nb_transitions(euler_pour_pot(etats_stables, T, nb_T, sigma, al)[1], etats_stables, ep))
        D.append(abs(V(etats_stables[0], al)))
    return np.array(D), np.array(tau)

def fit(x_data, y_data, fit):
    """
    Calcul de régressions sur des données

    Paramètres : 
    x_data : abscisse pour la régression
    y_data : ordonée pour la régression
    fit : paramètre de choix du type de régression

    Retours : 
    array : abscisse des données après fit (avec plus de valeurs)
    array : ordonnée des données après fit
    tuple : parametters ajustés
    list : variances associées à l'estimation des paramètres
    """

    # Vérifie la dimension des données
    if x_data.shape != y_data.shape:
        raise ValueError("Les données x_data et y_data doivent avoir la même forme.")
    
    # Vérifie les valeurs x égales à 0 pour éviter la division par zéro
    if np.any(x_data == 0) and fit == "inverse square":
        raise ValueError("x_data contains zero values, which would cause division by zero in the model.")

    abs_fit = np.linspace(np.min(x_data), np.max(x_data), (np.max(x_data)-np.min(x_data))*10000)

    # Ajustement avec curve_fit
    if fit == "exponential sig": 
        def exponential_function_sig(x, a, b):
            return a * np.exp(b / (x**2))
        (a, b), mat_cov = curve_fit(exponential_function_sig, x_data, y_data, p0=[150.0, 0.0], maxfev = 10000000)
        y_fit = exponential_function_sig(abs_fit, a, b)
        return np.array(abs_fit), np.array(y_fit), (a,b), np.sqrt(np.diag(mat_cov))

    # if fit == "inverse square":
    #     def inverse_square_function(x, a, b):
    #         return a / ((x+b)**2)
    #     (a, b), mat_cov = curve_fit(inverse_square_function, x_data, y_data, p0=(1.0, -1.0), maxfev = 10000)
    #     y_fit = inverse_square_function(abs_fit, a, b)
    #     return abs_fit, y_fit, (a, b), np.sqrt(np.diag(mat_cov))
    
    if fit == "exponential D": 
        def exponential_function_D(x, a, b):
            return a * np.exp(b*x)
        (a, b), mat_cov = curve_fit(exponential_function_D, x_data, y_data, p0=(1.0, 3.3), maxfev = 10000)
        y_fit = exponential_function_D(abs_fit, a, b)
        return np.array(abs_fit), np.array(y_fit), (a, b), np.sqrt(np.diag(mat_cov))
    
    if fit == "affine":
        def lineaire(x, a, b):
            return a*x + b
        (a, b), mat_cov = curve_fit(lineaire, x_data, y_data, p0=(1.0, 1.0), maxfev = 10000)
        y_fit = lineaire(abs_fit, a, b)
        return np.array(abs_fit), np.array(y_fit), (a, b), np.sqrt(np.diag(mat_cov))

    # if fit == "linéaire":
    #     def lineaire(x, a):
    #         return a*x
    #     a, pcov = curve_fit(lineaire, x_data, y_data, p0=(1.0), maxfev = 10000)
    #     y_fit = lineaire(abs_fit, a)
    #     return abs_fit, y_fit, a, np.sqrt(np.diag(pcov))


def lissage(signal_brut,L):
    """
    Lisse un signal grâce à une moyenne glissante 

    Paramètres : 
    signal_brut : données à lisser 
    L : nombre de données prises en compte pour la moyenne 

    Retours : 
    array : signal lissé
    """
    res = np.copy(signal_brut)
    for i in range (1,len(signal_brut)-1): 
        L_g = min(i,L) 
        L_d = min(len(signal_brut)-i-1,L) 
        Li=min(L_g,L_d)
        res[i]=np.sum(signal_brut[i-Li:i+Li+1])/(2*Li+1)
    return res