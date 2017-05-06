"""
Auteur: JJ Fleck, PCSI1, Kléber

Correction du sujet de l'X 2015 pour MP-PC sur les enveloppes convexes du plan.
"""

# On commence par l'exemple donné pour les tests:

tableau_exemple = [[0,1,1,4,4,5,5,7,7,8,11,13],
                   [0,4,8,1,4,9,6,-1,2,5,6,1]]

n = len(tableau_exemple[0])
assert n == len(tableau_exemple[1]) # On vérifie qu'il n'y ait pas d'erreur de frappe.

# Question 1 : Trouver le point de plus petite ordonnée
def plusBas(tab,n):
    jbas = 0
    for j in range(1,n): # On boucle sur les points suivants
        if   tab[1][j] <  tab[1][jbas]: # Si j'ai trouvé mieux
            jbas = j                    # je note
        # Si j'ai aussi bien avec un point plus à gauche
        elif tab[1][j] == tab[1][jbas] and tab[0][j] < tab[0][jbas]:
            jbas = j                    # je note aussi
    return jbas # Renvoi du bon point

# Vérification:
assert plusBas(tableau_exemple,n) == 7

# Question 3: Fonction d'orientation
# Pas vraiment de moyen simple de faire joli...
def orient(tab,i,j,k):
    a = tab[0][j] - tab[0][i] # Coordonnées de chaque vecteurs
    b = tab[1][j] - tab[1][i] # une
    c = tab[0][k] - tab[0][i] # à
    d = tab[1][k] - tab[1][i] # une
    determinant = a*d - b*c   # Le déterminant en question
    if   determinant > 0: return  1 # En fonction du signe du déterminant
    elif determinant < 0: return -1 # on détermine le signe
    else:                 return  0 # ou la nullité de l'aire

# Vérifications (question 2):
assert orient(tableau_exemple,0,3,4) ==  1
assert orient(tableau_exemple,8,9,10)== -1
assert orient(tableau_exemple,3,3,4) ==  0


# Question 5
def prochainPoint(tab,n,i,affiche=False):
    next = 0                         # On commence à 0
    if i==0: next = 1                # Si on regarde 0, il faut commencer à 1
    for j in range(n):               # On regarde chaque point
        # Si l'orientation est négative, c'est qu'on n'a pas le plus grand de la relation
        if orient(tab,i,next,j) < 0: # NB: la non égalité à 0 impose que i != j
            if affiche: print(next)  # Un peu de feedback pour les vérifications
            next = j                 # on note le meilleur candidat pour le moment
    return next # On renvoie le meilleur de tous les temps

# Vérification (question 6): 
#print('Devrait afficher 0, 1, 2 et 5')
#print(prochainPoint(tableau_exemple,n,10,affiche=True))

# Question 7: l'enveloppage proprement dit.
def convJarvis(tab,n):
    L = [plusBas(tab,n)]                  # On commence tout en bas
    candidat = prochainPoint(tab,n,L[-1]) # Le premier prochain point
    while candidat != L[0]:               # Tant qu'on n'a pas bouclé
        L.append(candidat)                # on rajoute le candidat
        candidat = prochainPoint(tab,n,candidat) # et on regarde le prochain candidat
    return L

# Vérification
assert convJarvis(tableau_exemple,n) == [7,11,10,5,2,0]

# Complexité: prochainPoint est linéaire en n et on doit l'exécuter autant de 
# fois qu'on trouve un point dans l'enveloppe finale (qui en contient m par 
# hypothèse), d'où O(n*m)

# On passe au second algorithme

# Question 9: Quicksort ou Fusionsort en n*log(n)

# On définit les piles à partir des listes Python de sorte à n'utiliser que 
# les primitives suivantes par la suite
def newStack(): return []
def isEmpty(s): return (len(s)==0)
def push(i,s): s.append(i)
def top(s): return s[-1]
def pop(s): return s.pop()

# Question 10: Mise à jour de l'enveloppe supérieure (ne renvoie rien mais 
# modifie la pile)
def majES(tab,es,i):
    if not(isEmpty(es)):  # Si on veut faire qc, il faut déjà que la pile soit non vide
        last = pop(es)    # On récupère le dernier
        # On ne s'arrête que si la pile est vide ou si l'orientation n'est plus la bonne
        while not(isEmpty(es)) and orient(tab,top(es),i,last) <= 0:
            last = pop(es)# On continue à dépiler
        push(last,es)     # On remet le dernier (enlevé en trop)
    push(i,es)            # On rajoute systématiquement le point d'étude

# Vérification:
es = [0,2,5,8]
majES(tableau_exemple,es,9)
assert es == [0,2,5,9]

# Question 11: Mise à jour de l'enveloppe inférieure (pareil qu'avant sauf 
# qu'on change de sens de rotation)
def majEI(tab,ei,i):
    if not(isEmpty(ei)):  # Si on veut faire qc, il faut déjà que la pile soit non vide
        last = pop(ei)    # On récupère le dernier
        # On ne s'arrête que si la pile est vide ou si l'orientation n'est plus la bonne
        while not(isEmpty(ei)) and orient(tab,top(ei),i,last) >= 0:
            last = pop(ei)# On continue à dépiler
        push(last,ei)     # On remet le dernier (enlevé en trop)
    push(i,ei)            # On rajoute systématiquement le point d'étude

# Vérification:
ei = [0,7,8]
majEI(tableau_exemple,ei,9)
assert ei == [0,7,9]

# Question 12:
def convGraham(tab,n):
    es = newStack()     # Initialisations des deux piles
    ei = newStack()
    for i in range(n):
        majES(tab,es,i) # Mise à jour au fur et à mesure
        majEI(tab,ei,i)
    pop(es) # On enlève le dernier point qui est en double
    while not(isEmpty(es)): # On prend l'enveloppe sup à rebrousse-poil
        val = pop(es)
        if val != 0:        # 0 est forcément en double aussi, on ne le rajoute donc pas
            push(val,ei)
    return ei # On a rajouté à l'enveloppe inférieure tous les points du dessus
    
# Vérification:
resultat = convGraham(tableau_exemple,n)
assert resultat == [0,7,11,10,5,2]

