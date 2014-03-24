#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "von_neumann.h"
#include "aes.h"
#include "mersenne_twister.h"

#define ARRAY_MAX_SIZE 1000
#define OLDRAND_MAX 2147483647
#define CALL_NBR 1024

static int next;

int rdtsc()
{
	// cette fonction suivante cree un warning : c'est normal.
	__asm__ __volatile__("rdtsc");
}

void oldinit_rand(int seed)
{
	next = seed;
}

int oldrand()
{
	next = next * 1103515245 + 12345;
	return (unsigned int)(next % OLDRAND_MAX);
}

// === TEST DE FRÉQUENCE MONOBIT ===
//Principe du test : 
//Le but de ce test est de s’intéresser à la proportion de zéros et de uns dans les bits d’une séquence entière : 
//on regarde ici tous les bits des 1024 réalisations de la séquence. On teste donc si le nombre de uns et de zéros 
//d’une séquence sont approximativement les mêmes comme attendu dans une séquence vraiment aléatoire.
//n: nombre de mots à traiter
//nbBits: nombre de bits par mots à traiter
//values: tableau des mots
double Frequency(int n,int nbBits, long * values)
{
	int j;
	int sum = 0;
	for(j=0;j<n;j++)
	{
		int i = 0;
		for(i=0;i<nbBits;i++)
		{
			if((values[j] >> i) & 0x01)
			{
				sum ++;
			}
			else
			{
				sum --;
			} 
		}
	}

	double s_obs = abs((double)sum)/sqrt((double)n*nbBits);

	return erfc(s_obs/sqrt(2.0));
}

// === TEST DES RUNS ===
//Principe du test : Le but de ce test est de s’intéresser à la longueur des suites successives de zéros et de uns dans la séquence observée. 
//Il teste donc la longueur moyenne de ce qu’on appelle les “runs”, i.e. les suites consécutives de 0 ou de 1.
//n: nombre de mots à traiter
//nbBits: nombre de bits par mots à traiter
//values: tableau des mots
double Runs(int n,int nbBits, long * values)
{
	int j = 0;
	int i = 0;
	//double nb_un =0;
	int sum = 0;
	int * valuesRand = (int *) calloc(n*nbBits,sizeof(int));
	for(j=0;j<n;j++)							//génération du tableau de bits et Pré-Test
	{
		for(i=0;i<nbBits;i++)
		{
			valuesRand[j*nbBits+i] = (values[j] >> i) & 0x01; 	// un bit par case
			sum += valuesRand[j*nbBits+i];				// TODO: se passer de valuesRand
		}
	}

	// Pré-Test
	double pi = ((double)sum)/(n*nbBits);
	double tau = 2.0/sqrt(n*nbBits);
	if(abs(pi-0.5)>tau)
	{
		return 0; 							// Fin du test
	}

	int vn_obs = 1; 							// init à 1 car Vn(obs) = Sum(...) +1
	
	int k = 0;
	for(k = 0; k <n*nbBits; k++)
	{
		if(valuesRand[k]!=valuesRand[k+1])
		{
			vn_obs++;
		}
	}

	return erfc( abs((double)vn_obs - 2*n*nbBits*pi*(1.0-pi)) / (2.0*sqrt(2.0*n*nbBits)*pi*(1.0-pi)));

	//for(j=0;j<n;j++)
	//{
	//	for(i=0;i<nbBits;i++)
	//	{
	//		if((values[j] >> i) & 0x01)
	//		{
	//			sum ++;
	//		}
	//		else
	//		{
	//			sum --;
	//		} 
	//	}
	//}

	//double s_obs = abs((double)sum)/sqrt((double)n*nbBits);

	//return erfc(s_obs/sqrt(2.0));
}
int main()
{
	word16 x=1111; 			// nombre entre 1000 et 9999 pour Von Neumann
	struct mt19937p mt; 		// Pour Mersenne-Twister
	int tmp = rand(), seed; 	// Pour Mersenne-Twister
	u32 Kx[NK], Kex[NB*NR], Px[NB]; // pour l'AES

	int output_rand; 		// sortie du rand du C
	int output_randFaible; 		// sortie du rand du C 
	int output_randFort; 		// sortie du rand du C
	word32 output_AES; 		// sortie pour l'AES
	word16 output_VN; 		// sortie pour pour Von Neumann
	word32 output_MT; 		// sortie pour Mersenne-Twister

               
	// initialisation des graines des generateurs

	srand(rdtsc());  		// rand du C 
	seed = rand();
	oldinit_rand(seed);
	sgenrand(time(NULL)+(tmp), &mt);// Mersenne-Twister

	// Initialisation de la clé et du plaintext pour l'AES 
	// 45 est un paramètre qui doit changer à chaque initialisation
	init_rand(Kx, Px, NK, NB, 45);
	KeyExpansion(Kex,Kx); 		// AES : sous-clefs

 	FILE* fichierOldRand = NULL;
    	fichierOldRand = fopen("oldRand.txt", "w+");

 	FILE* fichierOldRandFaible = NULL;
    	fichierOldRandFaible = fopen("oldRandFaible.txt", "w+");

 	FILE* fichierOldRandFort = NULL;
    	fichierOldRandFort = fopen("oldRandFort.txt", "w+");

 	FILE* fichierVN = NULL;
    	fichierVN = fopen("VN.txt", "w+");

 	FILE* fichierMT = NULL;
    	fichierMT = fopen("MT.txt", "w+");

 	FILE* fichierAES = NULL;
    	fichierAES = fopen("AES.txt", "w+");

	int i;
	for(i=0;i<1024;i++)
	{
		// sorties des generateurs	
		output_rand = oldrand(); 	// rand du C
		output_VN = Von_Neumann(&x); 	// Von Neumann
		output_MT = genrand(&mt); 	// Mersenne-Twister
		output_AES = AES(Px, Kex); 	// AES

		output_randFaible = output_rand & 0x0F;
		output_randFort = (output_rand >> 27) & 0x0F;
		
		//printf("- Generation de nombres aleatoires -\n");
		//printf("rand du C : %u \n",output_rand); 
		//printf("Von Neumann : %u\n",output_VN);
		//printf("Mersenne Twister : %u\n",output_MT);
		//printf("AES : %u\n",output_AES);
		
		// affichage
		fprintf(fichierOldRand,"%u\n",output_rand);
 		fprintf(fichierOldRandFaible,"%u\n",output_randFaible);
		fprintf(fichierOldRandFort,"%u\n",output_randFort);
		fprintf(fichierVN,"%u\n",output_VN);
		fprintf(fichierMT,"%u\n",output_MT);
		fprintf(fichierAES,"%u\n",output_AES);
	}
	fclose(fichierOldRand);
	fclose(fichierOldRandFaible);
	fclose(fichierOldRandFort);
	fclose(fichierVN);
	fclose(fichierMT);
	fclose(fichierAES);

	int size = 1024;
	long * valuesRand = (long *) calloc(size,sizeof(long));
	int k = 0;
	for(k=0;k<size;k++)
	{
		valuesRand[k] = AES(Px,Kex);
	}

	printf("Test Frequency: %f\n",Frequency(size,31,valuesRand));
	printf("Test Runs: %f\n",Runs(size,31,valuesRand));

	free(valuesRand);

	return 0;
}
