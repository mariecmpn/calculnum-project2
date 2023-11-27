#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "donnees.h"


double recup_L(float L) {
    /* fonction permettant de recuperer la valeur de la longueur (horizontale) du domaine */
    return 1.;
}

double recup_H(float H) {
    /* fonction permettant de recuperer la valeur de la longueur (verticale) du domaine */
    return 1.;
}

int recup_M(float M) {
    /* fonction permettant de recuperer le nombre de termes pour les sommes infinies que l'on a tronque */
    return 1;
}

int recup_n(float n) {
    /* fonction permettant de recuperer le nombre de points de quadrature utilises pour la methode de Gauss-Legendre */
    return 5;
}

int recup_nb_iter(int nb_iter) {
    /* fonction permettant de recuperer le nombre souhaite d'iterations pour les methodes d'approximation de f_3 */
    return 16;
}