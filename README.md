# Codes informatiques pour concours de CPGE

Banque de propositions d'implémentations pour les sujets proposés aux 
différentes épreuves d'informatique des concours en CPGE (a priori à la fois 
en info commune et en option info pour la filière MP).

## Notes pour les lecteurs

Les implémentations sont proposées "as is" sans garantie. Il s'agit de voir si 
les algorithmes proposés dans les sujets fonctionnent et parfois de les 
illustrer sur des exemples particuliers. Les contributions peuvent être plus 
ou moins commentées selon qu'il s'agissait en premier lieu de vérifier une 
intuition ou d'expliquer dans le détail à nos élèves, soyez donc compréhensifs 
et n'hésitez pas à proposer des améliorations/contributions/rapports de bugs 
si vous en détectez.

## Codes déjà disponibles

* [concours/Centrale/MP_PC_PSI_TSI/2015](concours/Centrale/MP_PC_PSI_TSI/2015/)

## Notes pour les contributeurs

L'idée est d'organiser le répertoire concours en sous-répertoires pour chaque concours:
* Banque_PT/
* CCP/
* Centrale/
* E3A/
* Mines/
* XENS/

Dans chaque sous-répertoire, on peut ajouter un sous-répertoire par filière 
(MP/, PC/, etc) ou par regroupement de filières (cas des sujets communs donc 
MP_PC/ pour l'X ou encore MP_PC_PSI_TSI/ pour Centrale)

On peut aussi imaginer un étage supplémentaire chez les MP pour prendre en 
compte les sujets d'option info qui sont distincts de ceux de l'informatique 
commune (mais les concours ayant option info et info commune ont souvent leur 
épreuve d'info commune en commun sur plusieurs filières donc cela ne devrait 
pas rajouter un étage).

Enfin, un dernier sous-répertoire avec l'année correspondante (2017/ par exemple)

Dans ce répertoire d'année, on peut mettre tous les fichiers pour le sujet en 
particulier tant qu'il n'y a qu'un seul contributeur. En revanche, si une 
autre personne veut proposer une autre implémentation, il faudra réorganiser 
les choses avec un sous-répertoire par contribution séparée.
