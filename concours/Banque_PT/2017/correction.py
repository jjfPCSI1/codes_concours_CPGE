# -*- coding: utf-8 -*-
"""
Created on Fri May  5 21:25:28 2017

@author: Kany

Proposition de correction pour l'épreuve de modélisation, banque PT, 2017
"""

#Q8

def indice(A,val):
    "méthode: recherche linéaire"
    "hypothèse: A[0]<=val and val<A[-1]"
    id_ = 0 #id_ car id est une built-in function 
    while A[id_] <= val:
        id_ += 1
    return id_ - 1

assert indice([0,1,2,3,4,5,6,7],3.5)==3
assert indice([0,1,2,3,4,5,6,7],3)==3
assert indice([0,1,2,3,4,5,6,7],0)==0
assert indice([0,1,2,3,4,5,6,7],6.99)==6

def indice(A,val):
    "méthode: recherche par dichotomie"
    "hypothèse: A[0]<=val and val<A[-1]"
    debut = 0
    fin = len(A)-1
    while fin-debut>1:
        milieu = (debut+fin)//2
        a = A[milieu]
        if a == val:
            return milieu
        elif a < val:
            debut = milieu
        else:
            fin = milieu
    return debut
    
assert indice([0,1,2,3,4,5,6,7],3.5)==3
assert indice([0,1,2,3,4,5,6,7],3)==3
assert indice([0,1,2,3,4,5,6,7],4.5)==4
assert indice([0,1,2,3,4,5,6,7],4)==4
assert indice([0,1,2,3,4,5,6,7],0)==0
assert indice([0,1,2,3,4,5,6,7],6.99)==6

#Q9
def extraire(T,P,Nm,i,j):
    ST=[[T[i][j],   T[i][j+1],   P[j],   Nm[i]],
        [T[i+1][j], T[i+1][j+1], P[j+1], Nm[i+1]]]
    return ST
    
#Q10
P = [0.295, 0.39, 0.48, 0.565, 0.645, 0.72, 0.79] #bar
Nm = [600,900,1300,1700,2200,2700,3200,3800,4400,5000,5600,6300,7000] #tr/min

T=[[2,420, 3.454, 4.408, 5.440, 6.484, 7.522, 8.548],
   [2.702, 3.776, 4.852, 5.900, 6.962, 8.004, 9.036],
   [3.064, 4.162, 5.248, 6.330, 7.418, 8.432, 9.434],
   [3.270, 4.412, 5.552, 6.644, 7.734, 8.774, 9.766],
   [3.432, 4.630, 5.804, 6.952, 8.062, 9.126, 10.154],
   [3.498, 4.710, 5.916, 7.090, 8.230, 9.334, 10.39],
   [3.560, 4.778, 6.022, 7.214, 3.388, 9.552, 10.614],
   [3.658, 4.878, 6.168, 7.358, 8.558, 9.742, 10.844],
   [3.740, 5.008, 6.324, 7.510, 8.738, 9.962, 11.066],
   [3.816, 5.134, 6.480, 7.708, 8.962, 10.204, 11.32],
   [3.908, 5.290, 6.622, 7.890, 9.174, 10.408, 11,582],
   [3.990, 5.390, 6.754, 8.076, 9.372, 10.612, 11,79],
   [4.030, 5.436, 6.808, 8.152, 9.452, 10.700, 11.898]]
   
#Un exemple au hasard
Pcol = 0.5 #bar
Nmot = 1600 #tr/min

j = indice(P,Pcol)
i = indice(Nm,Nmot)
ST= extraire(T,P,Nm,i,j)
print(ST)

#Q11
def Gauss(A,C):
    n = len(A)
    #on rend la matrice triangulaire
    for ligne in range(n):
        pivot = A[ligne][ligne]
        #il faudrait vérifier pivot != 0
        for ligne_dessous in range(ligne+1,n):
            coef = A[ligne_dessous][ligne]/pivot
            for colonne in range(ligne,n):
                A[ligne_dessous][colonne] -= A[ligne][colonne]*coef
            C[ligne_dessous] -= C[ligne]*coef
    #on resout à partir du bas
    solution = [0]*n
    for ligne in range(n-1,-1,-1):
        s = C[ligne]
        for colonne in range(ligne+1,n):
            s -= solution[colonne]*A[ligne][colonne]
        solution[ligne] = s/A[ligne][ligne]
    return solution
 
def test_Gauss(ST,Pcol,Nmot):
    A=[[1,0,0,0],[1,1,0,0],[1,0,1,0],[1,1,1,1]]
    C=[ST[0][0],ST[0][1],ST[1][0],ST[1][1]]
    b = Gauss(A,C)
    #print(b)
    
    Pj   = ST[0][2]
    Pjp1 = ST[1][2]
    Nmi  = ST[0][3]
    Nmip1 = ST[1][3]
    x = (Pcol-Pj)/(Pjp1-Pj)
    y = (Nmot-Nmi)/(Nmip1-Nmi)

    return b[0] + b[1]*x + b[2]*y + b[3]*x*y
 
print("interpolation par Gauss       :",test_Gauss(ST,Pcol,Nmot))
   
#Q14
def interpol(ST,Pcol,Nmot):
    Pj   = ST[0][2]
    Pjp1 = ST[1][2]
    Nmi  = ST[0][3]
    Nmip1 = ST[1][3]
    x = (Pcol-Pj)/(Pjp1-Pj)
    y = (Nmot-Nmi)/(Nmip1-Nmi)
    b = [0]*4
    b[0] = ST[0][0]
    b[1] = ST[0][1]-ST[0][0]
    b[2] = ST[1][0]-ST[0][0]
    b[3] = ST[1][1]-ST[0][1]-ST[1][0]+ST[0][0] #ERREUR DANS SUJET
    #print(b)
    
    return b[0] + b[1]*x + b[2]*y + b[3]*x*y

print("interpolation par substitution:",interpol(ST,Pcol,Nmot))  

#Q15
def sonde(r):
    "fonction seuil"
    lambda_ = 1/r #lambda_ car lambda est un mot reservé
    if lambda_ >= 1:
        return 0.1 #pauvre
    else:
        return 0.9 #riche
        
#Q16
def duree_injection(Usonde,Kpp,Kpn,Ki,tinjc0,integ_i,dt):
    """
    entrée : Usonde(t_i), integ_i
    sortie : Usonde(t_(i+1)), integ_(i+1)
    """
    
    integ_suivant = integ_i + tinjc0*dt    
    if Usonde <= 0.5:
        return tinjc0*Kpp + Ki*integ_suivant, integ_suivant
    else:
        return tinjc0*Kpn - Ki*integ_suivant, integ_suivant

#Q18
def Euler(tau,K,dt,s_i,e_i):
    """
    entrée : s(ti)
    sortie : s(ti+1)
    """
    return (K*e_i-s_i)*dt/tau + s_i
    
#Q19
def richesse(tau1,tau2,K,dt,r_i,w_i,tinjec_i):
    w_suivant = Euler(tau1,K,dt,w_i,tinjec_i)
    r_suivant = Euler(tau2,1,dt,r_i,w_i)
    return r_suivant,w_suivant
    
Tcycle = 60 / (2*4000) #s
tau1 = 10*Tcycle #s
tau2 = 20*Tcycle #s
Kpp = 1.03
Kpn = 0.95
Ki  = 0.8 #s**(-1)
tinjc0 = 12.3e-3 #s
t_0 = 0 #s
t_n = 100*Tcycle #s
dt  = 2e-6 #s
w_0 = 0.9 
K = 1/tinjc0 #s**(-1)
integ = 0

liste_richesse = [0.9]
liste_tinjec = [tinjc0]
liste_U = [sonde(liste_richesse[0])]

#Q21
liste_temps = [t_0+i*dt for i in range(int(t_n/dt)+1)]

r = liste_richesse[0]
t_injec = liste_tinjec[0]
U = liste_U[0]
w = w_0

cycle_mesure_U = int(Tcycle / dt)
print(cycle_mesure_U) #3750

for numero_t in range(1,len(liste_temps)):
    if numero_t % cycle_mesure_U == 0:
        U_nouveau = sonde(r)
        if U_nouveau != U: #franchissement de seuil
            integ = 0
        U = U_nouveau
    t_injec, integ = duree_injection(U,Kpp,Kpn,Ki,tinjc0,integ,dt)
    r, w = richesse(tau1,tau2,K,dt,r,w,t_injec)
    
    liste_U.append(U)
    liste_tinjec.append(t_injec)
    liste_richesse.append(r)

import matplotlib.pyplot as plt
plt.plot(liste_temps,liste_U,label="tension de sonde (V)")    
plt.xlabel("temps (s)")
plt.legend()
plt.grid()
plt.show()

plt.plot(liste_temps,liste_tinjec,label="durée d'injection (s)")
plt.xlabel("temps (s)")
plt.legend()
plt.grid()
plt.show()

plt.plot(liste_temps,liste_richesse,label="richesse")    
plt.xlabel("temps (s)")
plt.legend()
plt.grid()
plt.show()
    
   
