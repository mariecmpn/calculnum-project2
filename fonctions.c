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


/*************************
 fonctions a integrer pour
 les coefficients
**************************/

double B0_f3_f0_noyaux(double x, int m, double alpha){
    /* fonction a integer pour le coefficient B_0 pour la methode des noyaux iteres
    car en argument de gauss() il faut une fonction
    x: reel dont on veut calculer l'image */
    //double r = f_3(x,alpha)- f_0(x,m,alpha);
    double r = noyaux_iter(x,alpha) + m*0.;
    return r;
}

double Am_f3_noyaux(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m pour la methode des noyaux iteres
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
    double r = noyaux_iter(x,alpha)*cos((m*M_PI*x)/L);
    return r;
}

double Am_f0_noyaux(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m pour la methode des noyaux iteres
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

double Bm_f0_noyaux(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m pour la methode des noyaux iteres
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

double Bm_f3_noyaux(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m pour la methode des noyaux iteres
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
    double r = noyaux_iter(x,alpha)*cos((m*M_PI*x)/L);
    return r;
}

double B0_f3_f0_adomain(double x, int m, double alpha){
    /* fonction a integer pour le coefficient B_0 pour la methode de decomposition d'Adomain
    car en argument de gauss() il faut une fonction
    x: reel dont on veut calculer l'image */
    //double r = f_3(x,alpha)- f_0(x,m,alpha);
    double r = noyaux_iter(x,alpha) + m*0.;
    return r;
}

double Am_f3_adomain(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m pour la methode de decomposition d'Adomain
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
    double r = noyaux_iter(x,alpha)*cos((m*M_PI*x)/L);
    return r;
}

double Am_f0_adomain(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m pour la methode de decomposition d'Adomain
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

double Bm_f0_adomain(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m pour la methode de decomposition d'Adomain
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

double Bm_f3_adomain(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m pour la methode de decomposition d'Adomain
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
    double r = noyaux_iter(x,alpha)*cos((m*M_PI*x)/L);
    return r;
}

/************************
 fonctions coefficients
*************************/

double A_0_noyaux() {
    /* fonction qui calcule le coefficient A_0 pour l'approximation de T_tilde pour la methode des noyaux iteres */
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

double B_0_noyaux(double alpha) {
    /* fonction qui calcule le coefficient B_0 pour l'approximation de T_tilde pour la methode des noyaux iteres */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*H))*gauss(n,B0_f3_f0_noyaux,0,L,0,alpha);
    return r;
}

double A_m_noyaux(int m, double alpha) {
    /* fonction qui calcule le coefficient A_m pour l'approximation de T_tilde pour la methode des noyaux iteres */
    //on recupere d'abord les donnees du probleme
    float L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*sinh((m*M_PI*H)/L)))*(gauss(n,Am_f3_noyaux,0,L,m,alpha) - gauss(n,Am_f0_noyaux,0,L,m,alpha));
    return r;
}

double B_m_noyaux(int m, double alpha) {
    /* fonction qui calcule le coefficient B_m pour l'approximation de T_tilde pour la methode des noyaux iteres */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*sinh((m*M_PI*H)/L)))*(gauss(n,Bm_f0_noyaux,0,L,m,alpha) - gauss(n,Bm_f3_noyaux,0,L,m,alpha));
    return r;
}

double A_0_adomain() {
    /* fonction qui calcule le coefficient A_0 pour l'approximation de T_tilde pour la methode de decomposition d'Adomain */
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

double B_0_adomain(double alpha) {
    /* fonction qui calcule le coefficient B_0 pour l'approximation de T_tilde pour la methode de decomposition d'Adomain */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*H))*gauss(n,B0_f3_f0_adomain,0,L,0,alpha);
    return r;
}

double A_m_adomain(int m, double alpha) {
    /* fonction qui calcule le coefficient A_m pour l'approximation de T_tilde pour la methode de decomposition d'Adomain */
    //on recupere d'abord les donnees du probleme
    float L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*sinh((m*M_PI*H)/L)))*(gauss(n,Am_f3_adomain,0,L,m,alpha) - gauss(n,Am_f0_adomain,0,L,m,alpha));
    return r;
}

double B_m_adomain(int m, double alpha) {
    /* fonction qui calcule le coefficient B_m pour l'approximation de T_tilde pour la methode de decomposition d'Adomain */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    double r = (1/(L*sinh((m*M_PI*H)/L)))*(gauss(n,Bm_f0_adomain,0,L,m,alpha) - gauss(n,Bm_f3_adomain,0,L,m,alpha));
    return r;
}

/************************
 fonction T tilde approchee
*************************/

double T_tilde_noyaux(double x, double y, double alpha) {
    /* fonction qui donne l'approximation de T_tilde par la methode des noyaux iteres */
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
    T = A_0_noyaux()+ B_0_noyaux(alpha)*y;
    //printf("%f\n", T);
    for (i = 1; i<=M; i++) {
        T = T + (A_m_noyaux(i,alpha)*exp(i*M_PI*y/L) + B_m_noyaux(i,alpha)*exp(-i*M_PI*y/L))*cos(i*M_PI*x/L);
        //printf("%f\n", T);
    }
    return T;
}

double T_tilde_adomain(double x, double y, double alpha) {
    /* fonction qui donne l'approximation de T_tilde par la methode de decomposition d'Adomain */
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
    T = A_0_adomain()+ B_0_adomain(alpha)*y;
    //printf("%f\n", T);
    for (i = 1; i<=M; i++) {
        T = T + (A_m_adomain(i,alpha)*exp(i*M_PI*y/L) + B_m_adomain(i,alpha)*exp(-i*M_PI*y/L))*cos(i*M_PI*x/L);
        //printf("%f\n", T);
    }
    return T;
}

