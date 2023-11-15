#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h" // header de ce fichier
#include "donnees.h" // besoin des donnees du probleme
#include "methodesnum.h" // besoin de la methode de quadrature de gauss-legendre
#include "noyaux.h"


/************************
 fonctions conditions 
 aux limites et exactes
*************************/

double T_ex(double x, double y) {
    /* fonction solution exacte 
    x: reel dont on veut calculer l'image */
    double r = cosh(M_PI*y)*cos(M_PI*x);
    return r;
}

double f_0(double x, int m, double alpha) {
    /* fonction f_0 condition limite sur Gamma_0 
    x: reel dont on veut calculer l'image */
    double r = cos(M_PI*x) + 0.*alpha;
    return r;
}

double f_3_ex(double x) {
    /* fonction f_3 exacte condition limite sur Gamma_3 
    x: reel dont on veut calculer l'image */

    //on recupere d'abord les donnees du probleme
    double H;
    H = recup_H(H);
    // puis on retourne la fonction souhaitee 
    double r = cosh(M_PI)*cos(M_PI*x);
    return r;
}

double q_0(double x) {
    /* fonction q_0 exacte condition limite sur Gamma_0 */
    return 0.;
}

/************************
 fonctions a approcher
*************************/

double inte_f3(double x, int m, double alpha) {
    /* fonction que l'on integre dans l'expression de f_3 */
    //on recupere les donnees du probleme
    float L;
    L = recup_L(L);
    double r = cos(m*M_PI*x/L)*h(x)+ 0.*alpha;
    return r;
}

double f_3(double x, double alpha){
    /* fonction qui approche f_3 
    x: reel dont on veut calculer l'image par f_3
    alpha: parametre de Lavrentier */
    int i;
    double S;
    //on recupere les donnees du probleme
    double L,H;
    double a;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);

    a = (1./alpha);
    S = a*h(x) - a*(H/(alpha*H + 1.))*gauss(n,inte_h,0,L,0,alpha);
    //printf("%f\n",S);
    for (i=1; i <= M; i++) {
        S = S - a * (L*sinh((i*M_PI*H)/L) / (i*M_PI+alpha*L*sinh((i*M_PI*H)/L))) * cos((i*M_PI*x)/L) * gauss(n,inte_f3,0,L,i,alpha); 
    }
    return S;
}

/*************************
 fonctions a integrer pour
 les coefficients
**************************/

double B0_f3_f0(double x, int m, double alpha){
    /* fonction a integer pour le coefficient B_0
    car en argument de gauss() il faut une fonction
    x: reel dont on veut calculer l'image */
    //double r = f_3(x,alpha)- f_0(x,m,alpha);
    double r = f_3(x,alpha) + m*0.;
    return r;
}

double Am_f3(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m
    x: reel dont on veut calculer l'image
    m: indice m pour lequel on calcule A_m */
    
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = f_3(x,alpha)*cos((m*M_PI*x)/L);
    return r;
}

double Am_f0(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m
    x: reel dont on veut calculer l'image
    m: indice m pour lequel on calcule A_m */
    
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (f_0(x,m,alpha) * exp((-m*M_PI*H)/L))*cos((m*M_PI*x)/L);
    return r;
}

double Bm_f0(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m ()
    x: reel dont on veut calculer l'image
    m: indice m pour lequel on calcule B_m */
    
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = exp((m*M_PI*H)/L)*f_0(x,m,alpha)*cos((m*M_PI*x)/L);
    return r;
}

double Bm_f3(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m
    x: reel dont on veut calculer l'image
    m: indice m pour lequel on calcule B_m */
    
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = f_3(x,alpha)*cos((m*M_PI*x)/L);
    return r;
}

/************************
 fonctions coefficients
*************************/

double A_0() {
    /* fonction qui calcule le coefficient A_0 pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    //double r = (1/L)*gauss(n,f_0,0,L,0,0);
    return 0.;
}

double B_0(double alpha) {
    /* fonction qui calcule le coefficient B_0 pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*H))*gauss(n,B0_f3_f0,0,L,0,alpha);
    return r;
}

double A_m(int m, double alpha) {
    /* fonction qui calcule le coefficient A_m pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    float L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*sinh((m*M_PI*H)/L)))*(gauss(n,Am_f3,0,L,m,alpha) - gauss(n,Am_f0,0,L,m,alpha));
    return r;
}

double B_m(int m, double alpha) {
    /* fonction qui calcule le coefficient B_m pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*sinh((m*M_PI*H)/L)))*(gauss(n,Bm_f0,0,L,m,alpha) - gauss(n,Bm_f3,0,L,m,alpha));
    return r;
}

/************************
 fonction T tilde approchee
*************************/

double T_tilde(double x, double y, double alpha) {
    /* fonction qui donne l'approximation de T_tilde */
    double T;
    int i;
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    T = A_0()+ B_0(alpha)*y;
    //printf("%f\n", T);
    for (i = 1; i<=M; i++) {
        T = T + (A_m(i,alpha)*exp(i*M_PI*y/L) + B_m(i,alpha)*exp(-i*M_PI*y/L))*cos(i*M_PI*x/L);
        //printf("%f\n", T);
    }
    return T;
}

