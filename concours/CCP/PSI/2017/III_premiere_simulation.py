from numpy import *              # Imports présupposés par le sujet
from matplotlib.pyplot import *  # (mais mauvaise bonne idée en pratique)


diagramme_c_q = True

def debit(v_max,c_max,C_ligne):
    v = v_max * (1 - C_ligne/c_max)
    return C_ligne * v

def diagramme(v_max,c_max,C_ligne):
    Q = debit(v_max,c_max,C_ligne)
    plot(C_ligne,Q)

if diagramme_c_q:
    vmax = int(130 / 3.6)      # 130 km/h converti en m/s
    cmax = 1/3                 # Une voiture fait environ 3 m
    C = linspace(0,cmax,1000)
    Q = debit(vmax,cmax,C)
    diagramme(vmax,cmax,C)
    xlabel('Concentration $c$ en m$^{-1}$')
    ylabel('Debit $q$ en s$^{-1}$')
    xlim((0,round(2*cmax,1)/2))
    ylim((0,round(max(Q)+0.1,1)))
    title('Diagramme pour v_max$={}$ m/s et c_max={} m$^{{-1}}$'.format(vmax,round(cmax,2)))
    savefig('diagramme.png')
    clf()

def C_depart(dx,d1,d2,c1,c2,C):
    """ Renvoie le profil de concentration initiale sous forme d'un plateau
    rectangulaire de concentration c2 sur ]d1;d2] et de concentration c1 
    ailleurs. Les arguments sont
    * dx (float): pas d'espace
    * d1 (float): début du plateau
    * d2 (float): fin du plateau
    * c1 (float): concentration du plancher
    * c2 (float): concentration du plateau
    * C  (array): tableau bidimensionnel de largeur supérieure à d2//dx
    dont on veut initialiser la première ligne
    La fonction ne renvoie rien mais modifie au vol la tableau
    """
    n1,n2 = int(d1//dx), int(d2//dx)
    C[0,:] = c1         # On initialise tout à c1
    C[0,n1+1:n2+1] = c2 # Et on met à c2 les points de l'intervalle [[n1+1;n2]]

def resolution(C,dt,dx,c_max,v_max):
    p,n = C.shape # nombre de pas de temps (p) et d'espace (n)
    for i in range(0,p-1):    # On regarde instant après instant 
        # On commence par calculer le débit à l'instant i
        Q = debit(v_max,c_max,C[i])
        # On l'utilise alors pour calculer la concentration au temps suivant
        for j in range(n):    # en itérant sur tout l'espace
            Qjp1 = Q[(j+1)%n] # Conditions aux limites périodiques
            # Application de la formule
            C[i+1,j] = C[i,j] - (Qjp1 - Q[j]) / dx * dt
    return C # L'énoncé demande le renvoi (mais C est modifié anyway)

n = 6000  # Longueur
p = 6000 # Temps
dt= 0.001
dx= 0.1
La= n*dx
c1,c2 = cmax*0.7,cmax*0.9
d1,d2 = 2*La/5, 3*La/5

resolution_simple = True

if resolution_simple:
    C = zeros((p,n))
    C_depart(dx,d1,d2,c1,c2,C)
    resolution(C,dt,dx,cmax,vmax)
    x = np.linspace(0,La,n)
    plot(x,C[0],'--',label='Debut')
    nb = 4
    print(trapz(C[0]))
    for i in range(nb-1):
        print(trapz(C[(i+1)*p//nb-1]))
        plot(x,C[(i+1)*p//nb-1],'k')
    plot(x,C[-1],':',label='Fin')
    print(trapz(C[-1]))
    legend(loc='best')
    xlabel("Longueur sur l'autoroute en m")
    ylabel("Concentration en m$^{-1}$")
    savefig('resolution.png')
    clf()


def derive(f,dx,j):
    n = len(f)
    return (f[(j+1)%n] - f[(j-1)%n]) / (2*dx)

def resolution(C,dt,dx,c_max,v_max):
    p,n = C.shape # nombre de pas de temps (p) et d'espace (n)
    for i in range(0,p-1): # On regarde instant après instant
        # On commence par calculer le débit à l'instant i
        Q = debit(v_max,c_max,C[i])
        # On l'utilise alors pour calculer la concentration au temps suivant
        for j in range(n):    # en itérant sur tout l'espace
            # Evolution de la concentration
            Cij_moyen = (C[i,(j+1)%n]+C[i,(j-1)%n])/2
            C[i+1,j] = Cij_moyen - derive(Q,dx,j) * dt
    return C # L'énoncé demande le renvoi (mais C est modifié anyway)

c1 = 0.1*cmax
c2 = 0.3*cmax

resolution_raffinee = True

if resolution_raffinee:
    C = zeros((p,n))
    C_depart(dx,d1,d2,c1,c2,C)
    resolution(C,dt,dx,cmax,vmax)
    x = np.linspace(0,La,n)
    plot(x,C[0],'--',label='Debut')
    nb = 4
    print(trapz(C[0]))
    for i in range(nb-1):
        print(trapz(C[(i+1)*p//nb-1]))
        plot(x,C[(i+1)*p//nb-1],'k')
    plot(x,C[-1],':',label='Fin')
    print(trapz(C[-1]))
    legend(loc='best')
    xlabel("Longueur sur l'autoroute en m")
    ylabel("Concentration en m$^{-1}$")
    savefig('resolution2.png')
    clf()

