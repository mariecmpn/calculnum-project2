# recuperation des donnees pour le projet 2 - calcul numerique scientifique
# Partie 2 - determination de la frontiere libre

#%% importation packages

import numpy as np
import matplotlib.pyplot as plt
import os

#%% definition des variables

L = 1.
H = 1.

Alpha = np.zeros(22)

Alpha[0] = 10.
Alpha[1] = 1.
for i in range(2,22):
    Alpha[i] = 10**(-(i-1))
    
N = 51 #nombre de points de maillage pour x
    
X = np.linspace(0.,L,N)

#%% execution du programme partie1

#os.system('make clean')
#os.system('make partie2')
#os.system('./partie2') # on execute le fichier executable partie2 depuis notre console python
# penser a le compiler d'abord

#%% lecture des fichiers

# solution exacte sur Gamma_3
file = open('frontiere_ex.txt', 'r')

data = file.read()
Exact = data.split()
Ex = [float(i) for i in Exact]

#print(Ex)


# solution approchee sur Gamma_3
file = open('frontiere_app_noyaux.txt', 'r')

data = file.read()
Approche = data.split()
Approche = [float(i) for i in Approche]
App = np.array(Approche)

#print(App)

# On divise le tableau en 20 tableaux (1 tableau pour chaque alpha)
Gamma_alpha = np.split(App, 22)

# solution approchee sur Gamma_3
file = open('frontiere_app_adomain.txt', 'r')

data = file.read()
Approche = data.split()
Approche = [float(i) for i in Approche]
App_ado = np.array(Approche)

#print(App)

# On divise le tableau en 20 tableaux (1 tableau pour chaque alpha)
Gamma_alpha_ado = np.split(App_ado, 22)


 
#%% graphiques

alpha = 1
plt.plot(X, Ex, label = 'frontiere exacte')
plt.plot(X, Gamma_alpha[alpha], label = 'frontiere approchee avec noyaux')
plt.plot(X, Gamma_alpha_ado[alpha], label = 'frontiere approchee avec Adomain')
#plt.axis([0.,1.,0.,1.])
plt.legend()
plt.title('Frontieres Gamma exacte et approchee')
plt.show()