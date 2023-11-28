#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"
#include "fonctions2.h"

int main() {

    /***************************
     definition des variables 
    ****************************/ 

    int i, j; // entiers pour les boucles for
    int N = 51; // nombre de points de maillage en x
    double pas = 1./(N-1); // pas de maillage en x
    //double *Points = malloc(sizeof(int[N])); // tableau qui contient les valeurs des x_i
    double eps = 1.E-8; // tolerance epsilon pour la methode de Newton

    // recuperation des donnees du probleme definies dans donnees.c
    double L = recup_L(L), H = recup_H(H);
    int M = recup_M(M), n = recup_n(n);

    double Alpha;
    //double Alpha[22]; // tableau qui contient les differentes valeurs de alpha
    double Gamma_app;
    //double *Gamma_app = malloc(sizeof(int[N])); // tableau qui contient les images des x_i par T_alpha calculees numeriquement
    double *Gamma_exact = malloc(sizeof(int[N])); // tableau qui contient les images des x_i par T_ex calculees numeriquement
    int choix; // entier qui designe avec quelle methode on calcule f_3

    /***************************
     Choix de la methode de
     resolution de l'equation
     de Fredholm
    ****************************/

    printf("%s\n", "Choix de la methode de resolution de l'equation de Fredholm: ");
    printf("%s\n", "1: Methode de resolution d'Adomain");
    printf("%s\n", "2: Methode des noyaux iteres");
    printf("%s\n", "(Par defaut: methode d'Adomain)");
    printf("%s", "Choix = ");
    scanf("%d", &choix);

    /***************************
     remplissage des tableaux
    ****************************/ 

    // remplissage de Alpha
    Alpha = 10.;


    /***************************
     Calcul de la frontiere libre
     pour plusieurs alpha
    ****************************/ 

   double x_i;
   if (choix == 2) { // avec f_3 trouve avec la methode des noyaux iteres
        FILE *approche_noyaux;
        approche_noyaux = fopen("frontiere_app_noyaux.txt", "w");
        x_i = 0.;
        for (i=0; i<22;i++) {
                x_i = 0.;
                for (j=0; j<N; j++) {
                    //Points[j] = x_i;
                    Gamma_app = newton(x_i,fonction_T_noyaux,derivee_T_noyaux,eps,Alpha);
                    fprintf(approche_noyaux, "%g", Gamma_app); // on l'enregistre dans le fichier frontiere_app_noyaux.txt
                    fputs(" ", approche_noyaux);
                    x_i  = x_i+pas;
                }
                Alpha = Alpha*0.1;
                fputs("\n", approche_noyaux); // on change de ligne quand on change de alpha
        }
         fclose(approche_noyaux); // on ferme le fichier qu'on a ouvert
   }
    else { // avec f_3 trouve avec la methode de decomposition d'Adomain
        FILE *approche;
        approche = fopen("frontiere_app_adomain.txt", "w");
        for (i=0; i<22;i++) {
                x_i = 0.;
                for (j=0; j<N; j++) {
                    //Points[j] = x_i;
                    Gamma_app = newton(x_i,fonction_T_adomain,derivee_T_adomain,eps,Alpha);
                    fprintf(approche, "%g", Gamma_app); // on l'enregistre dans le fichier frontiere_app_adomain.txt
                    fputs(" ", approche);
                    x_i  = x_i+pas;
                }
                Alpha = Alpha*0.1;
                fputs("\n", approche); // on change de ligne quand on change de alpha
        }
         fclose(approche); // on ferme le fichier qu'on a ouvert
    }

    /***************************
     Calcul de la solution exacte
     ****************************/ 

    FILE *exact;
    exact = fopen("frontiere_ex.txt", "w");
    x_i = 0.;
    for (j = 0; j < N; j++) {
        //Points[j] = x_i;
        Gamma_exact[j] = Gamma_ex(x_i);
        fprintf(exact, "%g", Gamma_ex(x_i)); // on l'enregistre dans le fichier frontiere_ex.txt
        fputs(" ", exact);
        x_i  = x_i+pas;
    }

    // on ferme les fichiers qu'on a ouvert
    fclose(exact);

    // on desalloue l'espace memoire des tableaux alloues dynamiquement
    //free(Points);
    //free(Gamma_app);
    free(Gamma_exact);

    // on retourne 0
    return 0;
}