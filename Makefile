# Compiler used
# Changer le nom du compilateur si besoin
CC = gcc-11

# Regles pour la partie 1 du projet
partie1: main.o methodesnum.o donnees.o fonctions.o noyaux.o
	$(CC) -o partie1 main.o fonctions.o noyaux.o methodesnum.o donnees.o

methodenum.o: methodesnum.c donnees.h
	$(CC) -c methodesnum.c

donnees.o: donnees.c
	$(CC) -c donnees.c

fonctions.o: fonctions.c donnees.h methodesnum.h
	$(CC) -c fonctions.c

main.o: main.c noyaux.h adomain.h donnees.h methodesnum.h fonctions.h
	$(CC) -c main.c

noyaux.o: noyaux.c methodesnum.h donnees.h
	$(CC) -c noyaux.c

# Regles pour la partie 2 du projet
partie2: main2.o methodesnum.o donnees.o fonctions.o fonctions2.o
	$(CC) -o partie2 main2.o fonctions2.o fonctions.o methodesnum.o donnees.o

fonctions2.o: fonctions2.c donnees.h fonctions.h methodesnum.h
	$(CC) -c fonctions2.c

main2.o: main2.c fonctions2.h donnees.h fonctions.h methodesnum.h
	$(CC) -c main2.c 

# Regles pour la comparaison entre les methodes directe et indirecte d'Adomain

comp_adomain: main_adomain.o fonctions.o adomain.o methodesnum.o donnees.o
	$(CC) -o comp_adomain main_adomain.o fonctions.o adomain.o methodesnum.o donnees.o

main_adomain.o: main_adomain.c adomain.h donnees.h methodesnum.h fonctions.h
	$(CC) - c main_adomain.c

# Regle clean
clean: 
	rm *.o partie1 partie2