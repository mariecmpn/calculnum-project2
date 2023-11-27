#include "adomain.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"
#include "noyaux.h"

// PROGRAMME PRINCIPAL POUR LA COMPARAISON ENTRE LA METHODE DIRECTE ET INDIRECTE DE LA METHODE D'ADOMAIN

int main() {
    int i, j; // entiers pour les boucles for
    int Nx = 51; // nombre de points de maillage en x
    double pasx = 1./(Nx-1); // pas de maillage en x
    double a_dir, a_indir, err;

    // on fixe alpha = 1 pour ce test
    double alpha = 1.;

    // on ouvre les fichiers
    FILE *direct;
    FILE *indirect;
    FILE *erreur;
    direct = fopen("adomain_direct.txt", "w");
    indirect = fopen("adomain_indirect.txt", "w");
    erreur = fopen("adomain_erreur.txt", "w");


    double x_i = 0.;
    for (i = 1; i<= Nx; i++) {
        a_dir = Adomain(x_i, alpha);
        a_indir = adomain_direct(x_i, alpha);
        err = fabs(a_dir - a_indir);
        fprintf(direct, "%g", a_dir);
        fprintf(indirect, "%g", a_indir);
        fprintf(erreur, "%g", err);
        fputs(" ", direct);
        fputs(" ", indirect);
        fputs(" ", erreur);
        x_i = x_i + pasx;
    }


    fclose(direct);
    fclose(indirect);
    fclose(erreur);
    return 0;
}