#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "donnees.h" // besoin des donnees du probleme
#include "fonctions.h" // besoin des fonctions de la partie 1 pour calculer T_tilde
#include "methodesnum.h" // besoin de la methode de Newton
#include "fonctions2.h" // header de ce fichier
#include "noyaux.h"
#include "adomain.h"

double derivee_T_noyaux(double x, double y, double alpha) {
    /* fonction qui calcule la derivee de T_tilde(x,y) avec f_3 calculee par la methode des noyaux iteres
    x, y: coordonnees du point pour lequel on calcule la derivee
    alpha: parametre de Lavrentier utilise pour calculer f_3 */
    double S = B_0_noyaux(alpha);
    double L = recup_L(L);
    int i, M = recup_M(M);
    for (i = 1; i <= M; i++) {
        S = S + ((i*M_PI/L)*A_m_noyaux(i,alpha)*exp((i*M_PI/L)*y) - (i*M_PI/L)*B_m_noyaux(i,alpha)*exp((-i*M_PI/L)*y))*cos((i*M_PI/L)*x);
    }
    return S;
}

double derivee_T_adomain(double x, double y, double alpha) {
    /* fonction qui calcule la derivee de T_tilde(x,y) avec f_3 calculee par la methode d'Adomain
    x, y: coordonnees du point pour lequel on calcule la derivee
    alpha: parametre de Lavrentier utilise pour calculer f_3 */
    double S = B_0_adomain(alpha);
    double L = recup_L(L);
    int i, M = recup_M(M);
    for (i = 1; i <= M; i++) {
        S = S + ((i*M_PI/L)*A_m_adomain(i,alpha)*exp((i*M_PI/L)*y) - (i*M_PI/L)*B_m_adomain(i,alpha)*exp((-i*M_PI/L)*y))*cos((i*M_PI/L)*x);
    }
    return S;
}

double Gamma_ex(double x) {
    /* fonction de la frontiere libre exacte cherchee
    x: point de Gamma_0 */
    return 0.5*x + 0.3;
}

double T_0(double x) {
    /* fonction qui calcule T_0 pour un x donne
    on choisit une frontiere constante */
    double r = T_ex(x,Gamma_ex(x));
    return r;
}

double fonction_T_noyaux(double x, double y, double alpha) {
    /* fonction qui retourne la fonction utilisee pour Newton: T_tilde(x,y) - T_0(x) (pour un alpha et x donne)
    x,y: coordonnees dont on veut calculer l'image
    alpha: parametre de Lavrentier utilise dans l'approximation de T_tilde */
    double r = T_tilde_noyaux(x,y,alpha) - T_0(x);
    return r;
}

double fonction_T_adomain(double x, double y, double alpha) {
    /* fonction qui retourne la fonction utilisee pour Newton: T_tilde(x,y) - T_0(x) (pour un alpha et x donne)
    x,y: coordonnees dont on veut calculer l'image
    alpha: parametre de Lavrentier utilise dans l'approximation de T_tilde */
    double r = T_tilde_adomain(x,y,alpha) - T_0(x);
    return r;
}

