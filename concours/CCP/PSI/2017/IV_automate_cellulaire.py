""" 
Auteur: JJ Fleck, PCSI1, Kléber
    
Proposition d'implémentation des diverses fonctions demandées dans la partie 
IV du sujet PSI CCP 2017 sur la modélisation par automate cellulaire
de la circulation automobile sur une autoroute.
"""


from numpy import *               # Imports présupposés par le sujet
from matplotlib.pyplot import *   # (mais mauvaise bonne idée en pratique)
from random import random as rand # Pour la fonction rand() de l'énoncé

dx = 7.5                 # Pas en espace (en m)
dt = 1.2                 # Pas en temps  (en s)
vm = 130                 # Vitesse maximale (en km/h)
La = 8500                # Longueur de l'autoroute (en m)
N  = int(La // dx + 1)   # Nombre de cellules dans la simulation
p  = 160                 # Nombre de pas de temps
vcell = int(ceil(vm/3.6 * dt / dx))  # Vitesse en nombre de cellules par seconde
lance_simulation = True  # Veut-on lancer la simulation ?


def maj(Route,Vitesses,p,v_max,i):
    Vitesses_suivantes = array(Vitesses[i]) # Copie de l'état actuel
    N = len(Vitesses_suivantes)
    for j in range(N):
        if Route[i][j] == 1: # Si il y a une voiture dans la case
            # Étape 1: Accélération
            if Vitesses_suivantes[j] < v_max:
                Vitesses_suivantes[j] = Vitesses_suivantes[j] + 1
            # Étape 2: Décélération
            dn = distance(Route,i,j)
            if Vitesses_suivantes[j] > dn - 1:
                Vitesses_suivantes[j] = dn - 1
            # Étape 3: Facteur aléatoire
            if rand() < p and Vitesses_suivantes[j] > 0:
                Vitesses_suivantes[j] = Vitesses_suivantes[j] - 1
    return Vitesses_suivantes
        

def distance(Route,i,j):
    N = len(Route[i])
    d = 1 # On est au moins à une distance 1 de la prochaine case
    while Route[i][(j+d)%N] != 1: # Tant qu'on ne rencontre pas d'autre voiture,
        d = d+1                   # on continue à regarder devant soi
    return d

def deplacement(Vitesses,Route,Vitesses_suivantes,i):
    N = len(Vitesses_suivantes)
    for j in range(N):
        if Route[i][j] == 1: # Si il y a une voiture,
            # on prédit sa nouvelle position
            prochain = (j + Vitesses_suivantes[j])%N 
            # et on met à jour
            Route[i+1][prochain] = 1
            Vitesses[i+1][prochain] = Vitesses_suivantes[j]
    return Route, Vitesses

if lance_simulation:
    Route = zeros((p,N),dtype=int)     # Initialisation du tableau des routes successives
    Vitesses = zeros((p,N),dtype=int)  # Pareil pour les vitesses
    proba_presence = 0.15              # "Proba" de trouver une voiture
    intervalle = int(1/proba_presence) # Intervalle entre deux voitures successives
    proba_attention= 0.3               # Proba que le conducteur freine
    for j in range(N): # Initialisation des voitures
        # On démarre avec les voitures bien espacées à la vitesse maximale
        if j%intervalle == 0:
            Route[0][j] = 1
            Vitesses[0][j] = vcell 
    for i in range(p-1): # On passe chaque temps en revue
        v_next = maj(Route,Vitesses,proba_attention,vcell,i)
        Route,Vitesses = deplacement(Vitesses,Route,v_next,i)
    imshow(1-Route)
    xlim(0,250)
    ylim(150,0)
    xlabel('Indice $j$ des positions')
    ylabel('Indice $i$ des temps')  
    title("Apparition spontanée de bouchons\nen partant d'une distribution uniforme de voitures\navec toutes la même vitesse initiale (zoom)")
    savefig('bouchons.png')
