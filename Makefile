OBJ = mersenne_twister.o aes.o von_neumann.o

simul: main_debut.c $(OBJ)
	gcc -Wall -std=c99 -pedantic -o simul -lm main_debut.c $(OBJ) -lm #deux fois -lm car variable selon l'environnement

# Ces 2 lignes définissent la méthode de création d'un .o
.SUFFIXES: .o

.c.o:; gcc -Wall -std=c99 -pedantic -c -o $@ -lm $<

aes.o: aes.h aes.c

mersenne_twister.o: mersenne_twister.h mersenne_twister.c

von_neumann.o: von_neumann.h von_neumann.c

clean:
	rm *.o simul

