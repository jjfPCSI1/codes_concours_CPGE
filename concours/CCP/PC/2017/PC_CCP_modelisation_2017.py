# coding: latin1

"""
Auteur: JJ Fleck, PCSI1, Kléber

Proposition d'implémentation pour l'épreuve de modélisation de
Physique-Chimie, filière PC, concours CCP concernant la résolution de
l'équation de poisson dans un plan (x,y) pour une distribution volumique de
charges donnée.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

#plt.style.use('classic')

# Exemples que l'on veut exécuter:
cylindre_infini = True


# Q6
def nouveau_potentiel(V,rhos,frontiere,i,j):
    if frontiere[i,j]: 
        return V[i,j]
    else:
        return 0.25 * (V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1] + rhos[i,j])

# Q8
def itere_J(V,rhos,frontiere):
    N = len(V) - 1
    # On est obligé de passer par un intermédiaire, sinon les potentiels 
    # futurs se mélangent aux potentiels passés.
    Vnouveau = np.zeros(N+1,N+1)
    for i in range(N+1):
        for j in range(N+1):
            Vnouveau[i,j] = nouveau_potentiel(V,rhos,frontiere,i,j)
    # Calcul de l'erreur avec facilités Numpy (on aurait aussi pu sommer au fur et à mesure)
    erreur = np.sqrt(np.sum((Vnouvea - V) ** 2) / N**2)
    # On doit repasser tout en revue pour remplacer les valeurs
    # même si on aurait aussi pu écrire V[:,:] = Vnouveau
    for i in range(N+1):
        for j in range(N+1):
            V[i,j] = Vnouveau[i,j]
    # On renvoie l'erreur à la fin
    return erreur

# Q9
def poisson(f_iter,V,rhos,frontiere,eps):
    # On s'assure qu'on rentre initialement dans la boucle
    erreur = 2 * eps
    while erreur > eps: # et on n'en sort que lorsqu'on est satisfait
        erreur = f_iter(V,rhos,frontiere)


# Q11
def itere_GS(V,rhos,frontiere):
    # Plus besoin de potentiel intermédiaire puisque le futur est utilisé dans 
    # l'ordre où on le met à jour. En revanche, il faut faire attention pour 
    # le calcul de l'erreur.
    N = len(V) - 1
    erreur = 0
    for i in range(N+1):
        for j in range(N+1):
            Vold = V[i,j]
            V[i,j] = nouveau_potentiel(V,rhos,frontiere,i,j)
            erreur = erreur + (Vold - V[i,j])**2
    return np.sqrt(erreur / N**2)
    
# Q12
def nouveau_potentiel_SOR(V,rhos,frontiere,i,j,omega):
    # Toujours penser à utiliser les anciens.
    Vpredit = nouveau_potentiel(V,rhos,frontiere,i,j)
    return (1-omega) * V[i,j] + omega * Vpredit

# Q13
def itere_SOR(V,rhos,frontiere):
    # On reprend la structure de la Q11 avec le nouveau calcul du potentiel
    N = len(V) - 1
    omega_opt = 2 / (1 + np.pi/N)
    erreur = 0
    for i in range(N+1):
        for j in range(N+1):
            Vold = V[i,j]
            V[i,j] = nouveau_potentiel_SOR(V,rhos,frontiere,i,j,omega_opt)
            erreur = erreur + (Vold - V[i,j])**2
    return np.sqrt(erreur / N**2)
    
# Q14

# Il faut N itérations avec omega=omega_opt, sachant que chaque itération fait 
# intervenir un nombre proportionnel à N**2 d'opération, au final la 
# complexité est en N**3. Et efectivement, sur la courbe, quand N passe de 70 
# à 140 (soit une multiplication par 2), le temps d'exécution passe de 2s à 
# 16s (soit une multiplication par 8 = 2**3). Pour N = 1000 (soit 10 fois plus 
# que N=100), on devrait obtenir un temps d'environ 6000s (=2h !).

# Q15

# Il s'agit de calculer à chaque fois le gradient du potentiel en 
# s'intéressant soit aux différences sur une ligne (pour calculer Ex), soit à 
# celles sur une colonne (pour calculer Ey). On ne peut bien sûr calculer des 
# différences symétriques que si l'on n'est pas sur les bords 
# perpendiculairement à la direction de calcul.

def calc_ExEy(Ex,Ey,V,h):
    N = len(Ex)-1
    for i in range(N+1):
        for j in range(N+1):
            im,ip,jm,jp = i-1,i+1,j-1,j+1
            if i == 0: im = i
            if i == N: ip = i
            if j == 0: jm = j
            if j == N: jp = j
            Ex[i,j] = - (V[ip,j]-V[im,j])/(h*(ip-im))
            Ey[i,j] = - (V[i,jp]-V[i,jm])/(h*(jp-jm))

# Un premier exemple, celui du cylindre infini
if cylindre_infini:
    print('*'*70)
    print('Cylindre infini')
    # Initialisations
    eps0 = 8.85e-12   # epsilon_0
    R = 0.05           # Rayon du cylindre
    L = R*4       # Largeur L de la zone de visualisation
    xC,yC = L/2,L/2   # Centre du cylindre
    N = 100           # Nombre de points
    rho = 1e-5        # densité volumique de charge rho
    h = L/N               # Définition différente de celle de Q5
    rhos= rho/eps0 * h**2 # d'où une expression un peu différente pour rho''

    rhos_cyl = np.zeros((N+1,N+1))  # rho''
    V_cyl    = np.zeros((N+1,N+1))  # Potentiel
    Ex_cyl   = np.zeros((N+1,N+1))  # Composante du champ selon x
    Ey_cyl   = np.zeros((N+1,N+1))  # Composante du champ selon y
    
    frontiere_cyl = np.zeros((N+1,N+1),bool) # Initialement tout à False
    
    # Q19
    def dans_cylindre(x,y,xc,yc,R):
        return (x-xc)**2 + (y-yc)**2 <= R**2
        
    # Q20: 
    def initialise_rhos_cylindre(tab_rhos):
        # h, xC, yC, R et rhos sont des variables globales initialisées plus haut
        N = len(tab_rhos)-1
        for i in range(N+1):
            for j in range(N+1):
                x,y = h*i,h*j
                if dans_cylindre(x,y,xC,yC,R):
                    tab_rhos[i,j] = rhos

    # Q21: définition de la frontière
    def initialise_frontiere_cylindre(tab_f):
        N = len(tab_f)-1
        for i in range(N+1):
            for j in range(N+1):
                if i==0 or i==N or j==0 or j==N:
                    tab_f[i,j] = True
    
    # Faisons la résolution proprement dite
    print('Début des calculs')
    seuil = 1e-1
    initialise_rhos_cylindre(rhos_cyl)
    initialise_frontiere_cylindre(frontiere_cyl)
    poisson(itere_SOR,V_cyl,rhos_cyl,frontiere_cyl,seuil)
    calc_ExEy(Ex_cyl,Ey_cyl,V_cyl,h)
    print('Fin des calculs')
    
    # À présent, les représentations graphiques

    # D'abord les abscisses
    x = np.linspace(0,L,N+1)
    y = np.linspace(0,L,N+1)
    X,Y = np.meshgrid(x,y)
    cote = 6

    # Ensuite le premier graphique, en contour d'isopotentielles. On va faire 
    # des fonctions car on refait la même chose avec un cylindre légèrement 
    # différent juste après.    
    import matplotlib.cm as cm

    def graphe_rho(fichier,tab_rho):
        plt.figure(figsize=(cote,cote))
        plt.imshow(tab_rho, extent=(0, L, 0, L) )
        plt.title('Distribution de densite volumique')
        plt.xlabel('$x$ en m')
        plt.ylabel('$y$ en m')
        plt.savefig(fichier)
    graphe_rho('cyl_rhos.png',rhos_cyl)

    def graphe_isopot(fichier,V_cyl):
        plt.figure(figsize=(cote,cote))
        im = plt.imshow(V_cyl, interpolation='bilinear', origin='lower',
                    cmap=cm.jet, extent=(0, L, 0, L) )
        CS = plt.contour(X,Y,V_cyl)
        plt.clabel(CS,inline=1,fmt='%d')
        plt.title('Courbes isopotentielles (en volts)')
        plt.xlabel('$x$ en m')
        plt.ylabel('$y$ en m')
        plt.savefig(fichier)
    graphe_isopot('cyl_isoV.png',V_cyl)
    # On remarque que les courbes sont circulaires (comme on s'y attend) au 
    # centre, mais tende vers un carré vers les bords (imposé par les 
    # conditions au bord nulles sur la frontière)
    
    # Puis la coupe pour voir l'évolution du potentiel en passant par le centre
    def graphe_pot(fichier,V_cyl):
        plt.figure(figsize=(cote,cote))
        plt.plot(x,V_cyl[:,N//2+1])
        plt.title('Evolution du potentiel $V(x,y=L/2)$ en volts')
        plt.xlabel('$x$ en m')
        plt.ylabel('$V$ en volts')
        plt.xlim((min(x),max(x)))
        plt.ylim((0,1.05*np.max(V_cyl)))
        plt.savefig(fichier)
    graphe_pot('cyl_evoV.png',V_cyl)
    # On a une branche de parabole dans la zone interne au cylindre (dérivée 
    # seconde constante) et des droites dans les parties externes (dérivée 
    # seconde nulle)

    # Et finalement l'évolution du potentiel
    def graphe_Ex(fichier,Ex_cyl):
        plt.figure(figsize=(cote,cote))
        plt.plot(x,Ex_cyl[:,N//2+1])
        plt.title('Evolution du champ $E_x(x,y=L/2)$ en volts/metre')
        plt.xlabel('$x$ en m')
        plt.ylabel('$E_x$ en volts par metre')
        plt.xlim((min(x),max(x)))
        plt.ylim((1.05*np.min(Ex_cyl),1.05*np.max(Ex_cyl)))
        plt.savefig(fichier)
    graphe_Ex('cyl_evoEx.png',Ex_cyl)
    # Lui est linéaire là où le potentiel est parabolique à l'intérieur du 
    # cylindre et normalement nul en passant tout au centre (ce que confirme 
    # la valeur annoncée par l'énoncé, "petite" par rapport aux valeurs aux 
    # bords du cylindre)

    # Q23: on voit apparaître une troncature au centre où le potentiel est 
    # constant, phénomène bien connu qui avec ces symétries signifie que la 
    # densité de charge est retombée à 0. Partons des valeurs précédents et 
    # utilisons les facilités de Numpy pour annuler un cercle de charge 
    # volumique pour voir en quoi les calculs sont modifiés
    print('*'*70)
    print('Cylindre infini troué au centre')
    
    rhos_cyl[(X-xC)**2+(Y-yC)**2 < (R/2)**2] = 0

    graphe_rho('cyl2_rhos.png',rhos_cyl)

    print('Début des calculs')
    poisson(itere_SOR,V_cyl,rhos_cyl,frontiere_cyl,seuil)
    calc_ExEy(Ex_cyl,Ey_cyl,V_cyl,h)
    print('Fin des calculs')

    graphe_isopot('cyl2_isoV.png',V_cyl)
    graphe_pot('cyl2_evoV.png',V_cyl)
    graphe_Ex('cyl2_evoEx.png',Ex_cyl)
    # Cela correspond bien aux graphes de l'énoncé (à l'échelle verticale près)
