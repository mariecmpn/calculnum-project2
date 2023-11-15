#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fonctions.h"
#include "methodesnum.h"
#include "donnees.h"

double alpha = 1.;

double u0(double x){
	double res;
	res = (1./alpha)*h(x);
	return res;
}


/* Le fichier doit renvoyer la fonction f3 app par la méthode d'Adomain sous le format double f3_Ad(x) */

/* Créer une fonction pour intégrer la fonction calculée par la méthodes de Gauss */

