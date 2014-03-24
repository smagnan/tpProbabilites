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
#define ARRAY_MAX_SIZE 1000

static int next;

struct file_attente
{
	double *arrivee;
	int taille_arrivee;
	double *depart;
	int taille_depart;
};
typedef struct file_attente file_attente;

struct evolution
{
	double *temps;
	unsigned int *nombre;
};
typedef struct evolution evolution;

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

double Alea(u32 Kex[NB*NR], u32 Px[NB])
{
	double divide = 4294967296;
	word32 output_AES;
	output_AES = AES(Px, Kex); 	// AES
	return (double)output_AES/divide;

}

double Exponentielle(double lambda,u32 Kex[NB*NR], u32 Px[NB]) 
{
	return -log(1-Alea(Kex,Px))/lambda;
}

file_attente FileMM1(double lambda, double mu, double D,u32 Kex[NB*NR], u32 Px[NB])
{
	file_attente resultats;
	resultats.arrivee = (double *) calloc(ARRAY_MAX_SIZE,sizeof(double));
	resultats.taille_arrivee = 0;
	resultats.depart = (double *) calloc(ARRAY_MAX_SIZE,sizeof(double));
	resultats.taille_depart = 0;
	
	double time = Exponentielle(lambda,Kex,Px);
	resultats.arrivee[0] = time;
	resultats.taille_arrivee++;
	
	while(time<=D)
	{	
		time = Exponentielle(lambda,Kex,Px) + resultats.arrivee[resultats.taille_arrivee-1];
		if(time>D)
		{
			break;
		}
		resultats.arrivee[resultats.taille_arrivee] = Exponentielle(lambda,Kex,Px) + resultats.arrivee[resultats.taille_arrivee-1];
		resultats.taille_arrivee++;	
	}

	resultats.depart[0] = resultats.arrivee[0]+ Exponentielle(mu,Kex,Px);
	resultats.taille_depart++;
	
	int i;
	for(i=1 ; i< resultats.taille_arrivee; i++)
	{
		if(resultats.arrivee[i] > resultats.depart[i-1])
		{
			resultats.depart[i] = resultats.arrivee[i]+ Exponentielle(mu,Kex,Px);
		} 
		else
		{
			resultats.depart[i] = resultats.depart[i-1] +Exponentielle(mu,Kex,Px);
		}
	}

	resultats.taille_depart = resultats.taille_arrivee;

	return resultats;
}

evolution evolClient(file_attente file)
{
	evolution resultats;
	resultats.temps =  (double *) calloc(ARRAY_MAX_SIZE,sizeof(double));
	resultats.nombre =  (unsigned int *) calloc(ARRAY_MAX_SIZE,sizeof(unsigned int));
	
	int a = 0; // arrivees
	int d = 0; // departs
	resultats.temps[0] = file.arrivee[0];
	resultats.nombre[0]  = 1;
	a++;
	
	while((a+d)<=2*file.taille_arrivee)
	{

		if(file.arrivee[a]<file.depart[d] && a<file.taille_arrivee)// si le client suivant arrive avant que le client précédent soit parti
		{
			resultats.temps[a+d] = file.arrivee[a];
                        resultats.nombre[a+d] = resultats.nombre[a+d-1]+1;
			a++;
		}
		else
		{
			resultats.temps[a+d] = file.depart[d];
                        resultats.nombre[a+d] = resultats.nombre[a+d-1]-1;
			d++;
		}

	}


	//for(i = 0; i <file.taille_arrivee; i++)
	//{
	//	if(file.arrivee[i+1]<file.depart[i])  // si le client suivant arrive avant que le client précédent soit parti
	//	{
	//		resultats.temps[i+1] = file.arrivee[i+1];
	//		resultats.nombre[i+1] = resultats.nombre[i]+1;
	//	}
	//	else
	//	{
	//		resultats.temps[i+1] = file.depart[i];
        //                resultats.nombre[i+1] = resultats.nombre[i]-1;

	//	}
	//}

	return resultats;
}

double nbClientMoyen(evolution evol,int nbValeurs)
{
	int i = 0;
	double numerateur;
	double denomin;
	for(i = 1; i< nbValeurs;i++)
	{
		numerateur += (evol.temps[i]-evol.temps[i-1])*evol.nombre[i];
		denomin += evol.temps[i]-evol.temps[i-1];
	}

	return numerateur/denomin;
}

double attenteMoyClient(file_attente file)
{
	int i = 0;
	double sum;
	for(i = 1; i< file.taille_arrivee;i++)
        {
		sum += file.depart[i]-file.arrivee[i];
	}

	return sum/file.taille_arrivee;
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
	printf("Alea: %f\n",Alea(Px,Kex));
	printf("Exponentielle: %f\n",Exponentielle(8.0,Px,Kex));
	
	file_attente f_a = FileMM1(0.2,0.33,180,Px,Kex);
	evolution evol  = evolClient(f_a);
	
	for(i=0;i<f_a.taille_arrivee;i++)
	{
		printf("%f-%f\n",f_a.arrivee[i],f_a.depart[i]);
	} 

	for(i=0;i<2*f_a.taille_arrivee;i++)
        {
                printf("time: %f- val: %d\n",evol.temps[i],evol.nombre[i]);
        }

	printf("N_Moy: %f\n",nbClientMoyen(evol,2*f_a.taille_arrivee));
	printf("Attente moy: %f\n",attenteMoyClient(f_a));
	printf("lambda * attente moy =? N_Moy: %f\n",attenteMoyClient(f_a)*0.2);//TODO: mettre 0.2 dans une variable lambda
	free(valuesRand);

	return 0;
}
