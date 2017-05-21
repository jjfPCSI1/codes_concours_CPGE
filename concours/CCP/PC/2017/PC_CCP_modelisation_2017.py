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
import matplotlib.cm as cm
import math

#plt.style.use('classic')

# Exemples que l'on veut exécuter:
cylindre_infini = True
deviation_electron = True

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


def graphe_isopot(fichier,V_cyl,levels=None):
    plt.figure(figsize=(cote,cote))
    im = plt.imshow(V_cyl, interpolation='bilinear', origin='lower',
                cmap=cm.jet, extent=(0, L, 0, L) )
    CS = plt.contour(X,Y,V_cyl,levels=levels)
    plt.clabel(CS,inline=1,fmt='%d')
    plt.title('Courbes isopotentielles (en volts)')
    plt.xlabel('$x$ en m')
    plt.ylabel('$y$ en m')
    plt.savefig(fichier)

def graphe_pot(fichier,V_cyl):
        plt.figure(figsize=(cote,cote))
        plt.plot(x,V_cyl[:,N//2+1])
        plt.title('Evolution du potentiel $V(x,y=L/2)$ en volts')
        plt.xlabel('$x$ en m')
        plt.ylabel('$V$ en volts')
        plt.xlim((min(x),max(x)))
        plt.ylim((0,1.05*np.max(V_cyl)))
        plt.savefig(fichier)

def graphe_Ex(fichier,Ex_cyl):
        plt.figure(figsize=(cote,cote))
        plt.plot(x,Ex_cyl[:,N//2+1])
        plt.title('Evolution du champ $E_x(x,y=L/2)$ en volts/metre')
        plt.xlabel('$x$ en m')
        plt.ylabel('$E_x$ en volts par metre')
        plt.xlim((min(x),max(x)))
        plt.ylim((1.05*np.min(Ex_cyl),1.05*np.max(Ex_cyl)))
        plt.savefig(fichier)


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
    seuil = 1e-5
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

    def graphe_rho(fichier,tab_rho):
        plt.figure(figsize=(cote,cote))
        plt.imshow(tab_rho, extent=(0, L, 0, L) )
        plt.title('Distribution de densite volumique')
        plt.xlabel('$x$ en m')
        plt.ylabel('$y$ en m')
        plt.savefig(fichier)
    graphe_rho('cyl_rhos.png',rhos_cyl)

    graphe_isopot('cyl_isoV.png',V_cyl)
    # On remarque que les courbes sont circulaires au centre (comme on s'y 
    # attendrait s'il n'y avait que le cylindre uniformément chargé), mais 
    # tendent vers un carré vers les bords (imposé par les conditions au bord 
    # nulles sur la frontière)
    
    # Puis la coupe pour voir l'évolution du potentiel en passant par le centre
    graphe_pot('cyl_evoV.png',V_cyl)
    # On a une branche de parabole dans la zone interne au cylindre (dérivée 
    # seconde constante) et des droites dans les parties externes (dérivée 
    # seconde nulle)

    # Et finalement l'évolution du champ
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
    # Cela correspond bien aux graphes de l'énoncé 


# Deuxième simulation, un électron entre deux plaques d'un condensateur, base 
# des oscilloscopes ancienne génération
if deviation_electron:
    print('*'*70)
    print('Deviation electron')
    # Q24: la vitesse est donnée par v0 = sqrt(2*e*V0/m)
    # Q25: le champ est donné par -2Vp/d
    # Q26: trajectoire parabolique
    # Q27: plus de forces => ligne droite
    # Q28: rhos est nul en tout point (on a des charges surfaciques mais pas volumiques)
    L = 0.1
    D = 7e-2
    d = 2e-2
    ell = 4e-2
    N = 100
    h = L/N
    V0= 950
    Vp= 180
    m = 9.11e-31
    e = 1.6e-19
    rhos_osc = np.zeros((N+1,N+1))
    V_osc    = np.zeros((N+1,N+1))
    Ex_osc   = np.zeros((N+1,N+1))
    Ey_osc   = np.zeros((N+1,N+1))
    frontiere_osc = np.zeros((N+1,N+1),dtype=bool)

    # Q29
    def initialise_frontiere_condensateur(V,frontiere):
        """ Les armatures font ici partie de la frontière et il faut donc 
        imposer la valeur du potentiel sur chacune d'elle. """
        debut = int((L - D - ell/2) / L * N)
        fin   = int((L - D + ell/2) / L * N)
        h1    = N//2 - int(d/2 / L * N) # Signe pour que l'affichage corresponde bien
        h2    = N//2 + int(d/2 / L * N) # en toute fin avec les histoires de transposée
        V[debut:fin+1,h2] = Vp
        V[debut:fin+1,h1] =-Vp
        frontiere[debut:fin+1,h1] = True
        frontiere[debut:fin+1,h2] = True
        frontiere[0,:] = True
        frontiere[N,:] = True
        frontiere[:,0] = True
        frontiere[:,N] = True
    
    # Le calcul proprement dit
    print('Début des calculs')
    seuil = 1e-5
    initialise_frontiere_condensateur(V_osc,frontiere_osc)
    poisson(itere_SOR,V_osc,rhos_osc,frontiere_osc,seuil)
    calc_ExEy(Ex_osc,Ey_osc,V_osc,h)
    print('Fin des calculs')


    # D'abord les abscisses
    x = np.linspace(0,L,N+1)
    y = np.linspace(0,L,N+1)
    X,Y = np.meshgrid(x,y)
    cote = 6
    
    levels = np.arange(-Vp+10,Vp,20)

    # La transposition est nécessaire car le premier indice est interprété 
    # comme une ligne donc est vertical alors que dans l'énoncé on travaille 
    # depuis le début avec le premier indice qui correspond aux x. Ce n'était 
    # pas embêtant dans l'exemple du cylindre du fait de la symétrie globale 
    # du système.
    graphe_isopot('osc_isoV.png',V_osc.transpose(),levels)
    
    # Q30: Calcul des composantes du champ en n'importe quel point x ou y
    def val_Ex(Ex,Ey,x,y,h):
        i,j = int(x/h),int(y/h)
        if i+1 >= len(Ex) or j+1 >= len(Ex): return 0 # Au cas où l'on dépasse (pour debug)
        rx,ry = x-i*h, y-j*h
        E = Ex[i,j] + ((Ex[i+1,j]-Ex[i,j])*rx + (Ex[i,j+1]-Ex[i,j])*ry)/h
        return E
        
    def val_Ey(Ex,Ey,x,y,h):
        i,j = int(x/h),int(y/h)
        if i+1 >= len(Ex) or j+1 >= len(Ex): return 0 # Au cas où l'on dépasse (pour debug)
        rx,ry = x-i*h, y-j*h
        return Ey[i,j] + ((Ey[i+1,j]-Ey[i,j])*rx + (Ey[i,j+1]-Ey[i,j])*ry)/h
    
    # On définit aussi celles pour le champ théorique sans effets de bords
    def Ey_theo(Ex,Ey,x,y,h):
        if x > L - D - ell/2 and x < L - D + ell/2:
            return -2*Vp / d
        else:
            return 0

    def Ex_theo(Ex,Ey,x,y,h):
        return 0
    
    # Q33
    # Initialisations diverses
    Npts = 200              # Nombre de points voulus
    v0 = np.sqrt(2*e*V0/m)  # Vitesse initiale de l'électron
    ttraversee = L/v0       # Temps pour la traversée horizontale théorique
    dt = ttraversee / (Npts+1)  # dt associé pour nombre de points voulu
    # Pour pouvoir proprement faire aussi l'intégration théorique, on va 
    # encapsuler tout cela dans une fonction annexe à qui on donne les deux 
    # fonctions qui renvoient le champ électrique attendu au point voulu avec 
    # les mêmes arguments que val_Ex et val_Ey précédents.
    def euler_electron(fEx,fEy):
        lx = np.zeros(Npts)     # Tableau des coordonnées en x        
        ly = np.zeros(Npts)     # Tableau des coordonnées en y
        lvx= np.zeros(Npts)     # Tableau des vitesses en x
        lvy= np.zeros(Npts)     # Tableau des vitesses en y
        lx[0] = 0
        ly[0] = L/2             # La hauteur est non nulle
        lvx[0]= v0              # mais la vitesse est suivant x
        lvy[0]= 0
        # Intégration proprement dite
        for i in range(1,Npts):
            lx[i] = lx[i-1] + lvx[i-1]*dt
            ly[i] = ly[i-1] + lvy[i-1]*dt
            lvx[i]= lvx[i-1]- e/m * fEx(Ex_osc,Ey_osc,lx[i-1],ly[i-1],h)*dt
            lvy[i]= lvy[i-1]- e/m * fEy(Ex_osc,Ey_osc,lx[i-1],ly[i-1],h)*dt
        return lx,ly
    
    lx_simu,ly_simu = euler_electron(val_Ex,val_Ey)
    lx_theo,ly_theo = euler_electron(Ex_theo,Ey_theo)
    
    # Reste à faire un zouli dessin pour contempler tout cela et comparer à la 
    # figure de l'énoncé. (toujours les histoires de transposée pour inverser 
    # x et y à l'affichage)
    plt.figure(figsize=(cote,cote))
    #im = plt.imshow(V_osc.transpose(), interpolation='bilinear', origin='lower',
    #            cmap=cm.jet, extent=(0, L, 0, L) )
    CS = plt.contour(X,Y,V_osc.transpose(),levels=levels,colors='k')
    plt.clabel(CS,inline=1,fmt='%d')
    plt.plot(lx_simu,ly_simu,label='Simulation',linewidth=3)
    plt.plot(lx_theo,ly_theo,label='Théorie',linewidth=3)
    plt.legend(loc='best')
    plt.xlabel('$x$ en m')
    plt.ylabel('$y$ en m')
    plt.title("Déviation d'un électron,\ncalcul théorique du cours vs effets de bords")
    plt.savefig('osc_simulation_complete.png')
    
    # On retrouve ce que l'énoncé suggérait à savoir que la version 
    # "théorique" va plus loin car les effets de bords vont ralentir 
    # l'électron dans la simulation, le pas de temps étant calculé pour que la 
    # version théorique puisse aller au bout mais pas l'autre. En tous cas, en 
    # pratique, l'approximation de l'absence d'effet de bords n'est pas 
    # franchement correcte et le calcul du cours s'en trouve mis à mal...
