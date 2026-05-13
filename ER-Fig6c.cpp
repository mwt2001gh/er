//"Exonuclease ribozyme encoded in a circular genome: On the initiation of the RNA world"
//by Minglun Ling, Chunwu Yu, Hongyu Jiang, Wentao Ma*
//C source codes for the simulation program -- the version corresponds to the case shown in Fig. 6c of the paper

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <string.h>

/******* RANDOM NUMBER GENERATOR BY ZIFF **********/
#define A1 471
#define B1 1586
#define C1 6988
#define D1 9689
#define M 16383
#define RIMAX 2147483648.0        // = 2^31 
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-A1) & M] ^ ra[(nd-B1) & M] ^ ra[(nd-C1) & M] ^ ra[(nd-D1) & M])
void seed(long seed);  // Random number initialization 
static long ra[M + 1], nd;
/**************************************************/

#define SD 999  		   // Random seed
#define STEPNUM 10000000   // Total time steps of Monte Carlo simulation     
#define STAREC 0           // The step to start record
#define RECINT 10000       // The interval steps of recording

#define MAX_RNA_LENGTH 100   // The maximum RNA length allowed in the simulation
#define LONG_CHAIN_LEN 30    // The maximum RNA length for monitoring
#define MAX_CHAR_LENGTH 20   // The maximum RNA length of characteristic sequence of a ribozyme
#define MINI_CIRCLE_LENGTH 6 // The minimum RNA length for a circular RNA
#define LEN sizeof(struct rna)

#define C 2
#define G 3
#define A 1
#define U 4

//***********sequence****************
#define ERSEQ  A,C,G,A,A,C,U,G   // The characteristic sequence of exonuclease ribozyme (ER)
#define REPSEQ G,A,G,U,C,U,C,U   // The characteristic sequence of RNA replicase ribozyme (REP)
#define CTSEQ  G,U,G,A,U,C,C,A   // The characteristic sequence of an RNA control species 
#define NCSEQ  G,A,G,U,C,A,C,U   // The characteristic sequence of an assumed noncoding sequence (only 1 nt different from that of REP)
#define ERNCSEQ ERSEQ,NCSEQ      // The sequence of ER + the noncoding sequence in a circular genome inoculated in the beginning
#define CTNCSEQ CTSEQ,NCSEQ      // The sequence of a control for the circular genome inoculated in the beginning
#define ERREPSEQ ERSEQ,REPSEQ    // The sequence of a two-gene circular genome (ER+REP)
//***********sequence***************
#define INOCUSEQ  ERNCSEQ   // The inoculated species  
#define INOCUSEQ1 CTNCSEQ   // Another inoculated species

#define INOCUNUM 10      // The molecular number of the inoculated species

#define N 30			 // The side length of the two-dimensional grid
#define ROOMNUM (N*N)    // The room number in the grid
#define TNPB 100000      // Total nucleotide precursors introduced in the beginning

#define PSP 0.5          // Probability of separation of a base-pair
#define PBB  0.000001    // Probability of breaking of a phosphodiester bond
#define PRL 0.0000001    // Probability of random ligation between RNA molecules 
#define PEL  0.000005    // Probability of end-to-end ligation of an RNA molecule (circularization)  
#define PAT 0.5          // Probability of attracting a substrate by a template when the substrate could base-pair with the template
#define PTL 0.05         // Probability of a template-directed ligation (non-enzymatic)
#define PTLR 0.9		 // Probability of a template-directed ligation catalyzed by REP
#define PMNP 0.01		 // Probability of movement of nucleotide precursor
#define PMN 0.005        // Probability of movement of nucleotides
#define PNF 0.05         // Probability of a nucleotide forming from a nucleotide precursor (non-enzymatic)
#define PND 0.05         // Probability of a nucleotide decay to a nucleotide precursor
#define PNDE 0.001       // Probability of a nucleotide residue decaying at RNA's chain end
#define PFP 0.001        // probability of false base-pairing
#define PBER 0.9         // Probability of a phosphodiester bond breaking catalyzed by ER at the end of an RNA chain  
#define PRP 0.5          // Probability of a ribozyme being protected (from the cleaving of ER) due to its folding  
#define FDA 5            // The factor concerning the de novo attraction of a substrate onto an RNA template
#define FLT 0.5          // The factor concerning a linear RNA acting as a template
#define TER 10           // The times for an ER to function in a time step  
#define TREP 10			 // The times for a REP to function in a time step

#define RMRW (pow(p->length1+p->length2,1/2.0))    
						 // The relationship between an RNA's movement and its mass (the Zimm model)

long randl(long);         // Generating a random number (long) 
double randd(void);       // Generating a Random number (double)         
void avail_xy_init(void); // Initialization for xy_choose
void xy_choose(void);     // Picking a room at random 
void fresh_unit(void);    // Updating a unit for the next time step
int findseq(char seq[], int seqlength, struct rna* p); // Find a specific subsequence in an RNA sequence 
void inits(void);         // Initialization of the model system
void inoculate(char,int); // Inoculating specific RNA species
void unit_case(void);     // Action of units (molecules) in the system
void record(void);        // Data recording
void freepool(void);      // Memory releasing

struct rna                // Representation of an RNA (including nucleotide) molecule
{
	char information[2][MAX_RNA_LENGTH];  // Sequences of the RNA and its complementary chain
	int length1;
	int length2;
	struct c2_frag* chain2;  // A substrate (nucleotide or oligonucleotide) on the RNA (as a template) 
	char type1;              // Topological type of the RNA: circular (1) or linear (0)
	char type2;              // Topological type of the RNA's complementary chain: circular (1) or linear (0)
	struct rna* next;
	struct rna* prior;
};
struct rna* room_head[2][N][N];
struct rna* p, * p1, * p2, * p3, * p4, * ps, * ps1, * ps2;

struct c2_frag          // Representation of a substrate (nucleotide or oligonucleotide) on an RNA template 
{
	int start;
	int length;
	struct c2_frag* next;
	struct c2_frag* prior;
};
struct c2_frag* chain2;
struct c2_frag* c2f, * c2f1, * c2f2, * c2f3, * c2f4;

static char erseq[50] = { ERSEQ };
static char repseq[50] = { REPSEQ };
static char ctseq[50] = { CTSEQ };

static char erncseq[50] = { ERNCSEQ };
static char ctncseq[50] = { CTNCSEQ };
static char errepseq[50] = { ERREPSEQ };

static char inocuseq[50] = { INOCUSEQ };
static char inocuseq1[50] = { INOCUSEQ1 };

static int np_arr[N][N];
char temp_information[2][MAX_RNA_LENGTH];

int over_max_len = 0;
int x, y;                 // The coordinate of rooms in the grid 
int erlength, replength, ctlength, ernclength, ctnclength, erreplength;
int inoculength, inoculength1;
int randcase, randcase1, randcaser, randcaser1, length3, g = 0, h = 0, g_end = 0, gi, g_stop = 0;
int flag, flag1, flag2, flag3, flag4, flag5;
long i;                  // Cycle variable for Monte Carlo steps
long available;
long availabl[ROOMNUM];

float er[(STEPNUM - STAREC) / RECINT + 1];  // For data recording in time steps
float er_rib[(STEPNUM - STAREC) / RECINT + 1];
float cir_er[(STEPNUM - STAREC) / RECINT + 1];

float rep[(STEPNUM - STAREC) / RECINT + 1];  
float rep_rib[(STEPNUM - STAREC) / RECINT + 1];
float cir_rep[(STEPNUM - STAREC) / RECINT + 1];

float ct[(STEPNUM - STAREC) / RECINT + 1];
float cir_ct[(STEPNUM - STAREC) / RECINT + 1];

float ernc[(STEPNUM - STAREC) / RECINT + 1];  
float cir_ernc[(STEPNUM - STAREC) / RECINT + 1];

float ctnc[(STEPNUM - STAREC) / RECINT + 1];  
float cir_ctnc[(STEPNUM - STAREC) / RECINT + 1];

float errep[(STEPNUM - STAREC) / RECINT + 1];  
float cir_errep[(STEPNUM - STAREC) / RECINT + 1];

float total_mat_num[(STEPNUM - STAREC) / RECINT + 1];
float unit[(STEPNUM - STAREC) / RECINT + 1];
float cir_unit[(STEPNUM - STAREC) / RECINT + 1];
float np_num[(STEPNUM - STAREC) / RECINT + 1];


void seed(long seed)  // Random number initialization
{
	int a;

	if (seed < 0) { puts("SEED error."); exit(1); }
	ra[0] = (long)fmod(16807.0 * (double)seed, 2147483647.0);
	for (a = 1; a <= M; a++)
	{
		ra[a] = (long)fmod(16807.0 * (double)ra[a - 1], 2147483647.0);
	}
}

//------------------------------------------------------------------------------
long randl(long num)  // Generating a random number (long) 
{
	return(RandomInteger % num);
}

//------------------------------------------------------------------------------
double randd(void)        // Generating a random number (double) 
{
	return((double)RandomInteger / RIMAX);
}

//------------------------------------------------------------------------------
void avail_xy_init(void)   // Initialization for xy_choose
{
	int j;
	for (j = 0; j < ROOMNUM; j++)
	{
		availabl[j] = j + 1;
	}
	available = ROOMNUM;
}

//------------------------------------------------------------------------------
void xy_choose(void)       // Picking a grid room at random
{
	long rl, s;
	rl = randl(available);
	s = availabl[rl];
	x = (s - 1) % N;
	y = (s - 1) / N;
	availabl[rl] = availabl[available - 1];
	available--;
}

//------------------------------------------------------------------------------
void fresh_unit(void)     // Updating a unit for the next time step
{
	p1 = p->prior;
	p2 = p->next;
	p3 = room_head[!h][y][x]->next;
	room_head[!h][y][x]->next = p;
	p->next = p3;
	p->prior = room_head[!h][y][x];
	if (p3 != room_head[!h][y][x])p3->prior = p;
	p1->next = p2;
	if (p2 != room_head[h][y][x])p2->prior = p1;
	p = p1;
}

//------------------------------------------------------------------------------
int findseq(char seq[], int seqlength, struct rna* p)  // Find a specific subsequence in an RNA sequence
{
	int flag2, a, b;

	char inf[MAX_RNA_LENGTH + MAX_CHAR_LENGTH];
	for (a = 0; a < MAX_RNA_LENGTH + MAX_CHAR_LENGTH; a++)inf[a] = 0;

	for (a = 0; a < p->length1 + seqlength; a++)
	{
		if (a < p->length1) inf[a] = p->information[0][a];
		else inf[a] = p->information[0][a - p->length1];
	}

	flag2 = 0;

	if (p->length1 >= seqlength)
	{
		if (p->type1 == 0)
		{
			for (b = 0; p->length1 - seqlength - b >= 0; b++)
			{
				flag2 = 0;
				for (a = 0; a < seqlength; a++)
				{
					if (inf[b + a] == seq[a])continue;
					else { flag2 = 1; break; }  // No this subsequence at this location
				}
				if (flag2 == 0)break; // Containing this subsequence at this location
			}
		}
		else if (p->type1 == 1)
		{

			for (b = 0; b <= p->length1 - 1; b++)
			{
				flag2 = 0;
				for (a = 0; a < seqlength; a++)
				{
					if (inf[b + a] == seq[a])continue;
					else { flag2 = 1; break; }
				}
				if (flag2 == 0)break;
			}
		}
	}
	else flag2 = 1;

	if (flag2 == 0)return(0);   // Yes, the sequence contains the subsequence
	else return(1);   // No, the sequence does not contain the subsequence
}

//------------------------------------------------------------------------------
void inits(void)         // Initialization of the model system
{
	int j, m, k;
	seed(SD);

	erlength = 0;
	for (j = 0; erseq[j] != 0; j++)
		erlength++;
	replength = 0;
	for (j = 0; repseq[j] != 0; j++)
		replength++;
	ctlength = 0;
	for (j = 0; ctseq[j] != 0; j++)
		ctlength++;

	ernclength = 0;
	for (j = 0; erncseq[j] != 0; j++)
		ernclength++;
	ctnclength = 0;
	for (j = 0; ctncseq[j] != 0; j++)
		ctnclength++;
	erreplength = 0;
	for (j = 0; errepseq[j] != 0; j++)
		erreplength++;

	inoculength = 0;
	for (j = 0; inocuseq[j] != 0; j++)
		inoculength++;

	inoculength1 = 0;
	for (j = 0; inocuseq1[j] != 0; j++)
		inoculength1++;

	for (m = 0; m < 2; m++)
	{
		for (y = 0; y < N; y++)
		{
			for (x = 0; x < N; x++)
			{
				p1 = (struct rna*)malloc(LEN);
				if (!p1) { printf("\tinit1--memeout\n"); exit(0); }
				room_head[m][y][x] = p1;
				p1->next = room_head[m][y][x];
			}
		}
	}

	for (k = 0; k < TNPB; k++)  // Initial introduction of nucleotide precursors
	{
		x = randl(N);
		y = randl(N);
		np_arr[y][x]++;
	}
}

////////////////////////////////////
void inoculate(char ch_typ, int num)   // Inoculating specific RNA species
{
	int k, k1;

	for (k = 0; k < num; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		memset(p2->information, 0, sizeof(p2->information));
		for (k1 = 0; k1 < inoculength; k1++) p2->information[0][k1] = inocuseq[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength;
		p2->length2 = 0;
		p2->type1 = ch_typ;
		p2->type2 = 0;
		p2->next = room_head[h][y][x]->next;
		if (p2->next != room_head[h][y][x])(p2->next)->prior = p2;
		room_head[h][y][x]->next = p2;
		p2->prior = room_head[h][y][x];

		c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
		if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
		p2->chain2 = c2f;
		c2f->next = p2->chain2;
		c2f->prior = p2->chain2;
	}

	for (k = 0; k < num; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		memset(p2->information, 0, sizeof(p2->information));
		for (k1 = 0; k1 < inoculength1; k1++) p2->information[0][k1] = inocuseq1[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength1;
		p2->length2 = 0;
		p2->type1 = ch_typ;
		p2->type2 = 0;
		p2->next = room_head[h][y][x]->next;
		if (p2->next != room_head[h][y][x])(p2->next)->prior = p2;
		room_head[h][y][x]->next = p2;
		p2->prior = room_head[h][y][x];

		c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
		if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
		p2->chain2 = c2f;
		c2f->next = p2->chain2;
		c2f->prior = p2->chain2;
	}
}

//------------------------------------------------------------------------------
void unit_case(void)      // Action of units (molecules) in the system
{
	int a, b, c, d, j, k, randnt, np_bef, nt_turn, rep_turn, er_turn, length, randseq;
	double f, f1, rtdaddlig, rtdaddphili;

	avail_xy_init();      // Initialization for xy_choose
	for (d = 0; d < ROOMNUM; d++)
	{
		xy_choose();     // Picking a room at random 
		np_bef = np_arr[y][x];     //Events of nucleotide precursors
		for (k = 0; k < np_bef; k++)
		{
			randcaser = randl(2);
			switch (randcaser)
			{
			case 0:  // Forming nt
				if (randd() < PNF)
				{
					np_arr[y][x]--;
					p3 = (struct rna*)malloc(LEN);
					if (!p3) { printf("\t%d form_monomer--memeout\n", k + 1); exit(0); }
					memset(p3->information, 0, sizeof(p3->information));
					randnt = randl(4) + 1;
					switch (randnt)
					{
					case 1:  p3->information[0][0] = A; break;
					case 2:  p3->information[0][0] = C; break;
					case 3:  p3->information[0][0] = G; break;
					case 4:  p3->information[0][0] = U; break;
					default: printf("form randnt error");
					}
					p3->length1 = 1;
					p3->length2 = 0;
					p3->type1 = 0;
					p3->type2 = 0;

					p3->prior = room_head[!h][y][x];
					p3->next = room_head[!h][y][x]->next;
					if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
					room_head[!h][y][x]->next = p3;

					c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
					if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
					p3->chain2 = c2f;
					c2f->next = p3->chain2;
					c2f->prior = p3->chain2;
				}
				break;

			case 1:   // Nucleotide precursors moving
				if (randd() < PMNP)
				{
					np_arr[y][x]--;
					randcaser1 = randl(4);   // Four possible directions
					switch (randcaser1)
					{
					case 0:
						np_arr[y][(N + x - 1) % N]++;  // Toroidal topology
						break;
					case 1:
						np_arr[y][(x + 1) % N]++;
						break;
					case 2:
						np_arr[(N + y - 1) % N][x]++;
						break;
					case 3:
						np_arr[(y + 1) % N][x]++;
						break;
					default: printf("np moving error");
					}
				}
				break;

			default:printf("np case error");
			}
		}

		for (p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)  //Events of RNA (including nucleotides)
		{
			randcase = randl(6);
			switch (randcase)
			{
			case 0:       // RNA ligation
				if (p->type1 == 0)    // Linear chain
				{
					if (p->length1 >= MINI_CIRCLE_LENGTH)   // Intramolecular ligation (circularization) 
					{
						if (randd() < PEL)
						{
							p->type1 = 1;
							fresh_unit();
							break;
						}
					}

					for (p3 = p->next; p3 != p; p3 = p3->next)  // Intermolecular ligation
					{
						if (p3 == room_head[h][y][x]) { p3 = room_head[h][y][x]->next; if (p3 == p)break; }
						if (p3->type1 == 0)
						{
							if (randd() < PRL / (p->length1 * p3->length1))   // Longer chains should be more difficult to ligate
							{
								if (p->length1 + p3->length1 > MAX_RNA_LENGTH - 1)
								{
									over_max_len++; continue;
								}

								for (a = 0; a < p3->length1; a++)
								{
									p->information[0][a + p->length1] = p3->information[0][a];
									p->information[1][a + p->length1] = p3->information[1][a];
								}

								for (c2f3 = p3->chain2->next; c2f3 != p3->chain2; c2f3 = c2f3->next)
								{
									c2f3->start += p->length1;
								}
								p->chain2->prior->next = p3->chain2->next;
								p3->chain2->next->prior = p->chain2->prior;
								p->chain2->prior = p3->chain2->prior;
								p3->chain2->prior->next = p->chain2;

								p->length1 = p->length1 + p3->length1;
								p->length2 = p->length2 + p3->length2;

								(p3->prior)->next = p3->next;
								if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
								free(p3->chain2);
								free(p3);
								break;
							}
						}
					}
				}
				fresh_unit();
				break;

			case 1:            // Decay and degradation 
				if (p->length1 == 1)  // Decay of nucleotides
				{
					if (p->length2 == 0)
					{
						if (randd() < PND)
						{
							np_arr[y][x]++;
							(p->prior)->next = p->next;
							if (p->next != room_head[h][y][x])(p->next)->prior = p->prior;
							p3 = p;
							p = p->prior;
							free(p3->chain2);
							free(p3); break;
						}
					}
					else if (p->length2 == 1)    // The decay of the paired nucleotides at the same time
					{
						if (randd() < PND * sqrt(PND))
						{
							np_arr[y][x] += 2;
							(p->prior)->next = p->next;
							if (p->next != room_head[h][y][x])(p->next)->prior = p->prior;
							p3 = p;
							p = p->prior;
							free(p3->chain2);
							free(p3); break;
						}
					}
					else { printf("unexpected error on chain-length"); exit(0); }
				}
				else                  //Degradation of RNA
				{
					if (p->type1 == 0)  // Linear chain
					{
						c2f1 = p->chain2->prior;

						// Nucleotide residue decaying at the end of RNA
						if (p->information[1][p->length1 - 1] == 0)         // Single chain at the end 
						{
							if (randd() < PNDE)
							{
								p->information[0][p->length1 - 1] = 0;
								p->length1--;
								np_arr[y][x]++;
							}
						}
						else if (p->information[1][p->length1 - 1] != 0)  // Double chain at the end
						{
							if (randd() < PNDE * sqrt(PNDE))                 // The decay of the paired residues at the same time
							{
								p->information[0][p->length1 - 1] = 0;
								np_arr[y][x]++;

								p->information[1][p->length1 - 1] = 0;
								np_arr[y][x]++;

								p->length1--;
								p->length2--;
								c2f1->length--;

								if (c2f1->length == 0)
								{
									(c2f1->prior)->next = c2f1->next;
									(c2f1->next)->prior = c2f1->prior;
									free(c2f1);
								}
							}
						}

						if (p->length1 == 1)
						{
							fresh_unit();
							break;
						}

						// Nucleotide residue decaying at the start of RNA
						c2f2 = p->chain2->next;
						if (p->information[1][0] == 0)      // Single chain at the start
						{
							if (randd() < PNDE)
							{
								for (b = 1; b < p->length1; b++)
								{
									p->information[0][b - 1] = p->information[0][b];
									p->information[1][b - 1] = p->information[1][b];
								}
								p->information[0][p->length1 - 1] = 0;
								p->information[1][p->length1 - 1] = 0;
								p->length1--;
								np_arr[y][x]++;

								for (c2f3 = c2f2; c2f3 != p->chain2; c2f3 = c2f3->next)
								{
									c2f3->start--;
								}
							}
						}
						else if (p->information[1][0] != 0) // Double chain at the start
						{
							if (randd() < PNDE * sqrt(PNDE)) { // The decay of paired residues at the same time

								for (b = 1; b < p->length1; b++)
								{
									p->information[0][b - 1] = p->information[0][b];
									p->information[1][b - 1] = p->information[1][b];
								}
								p->information[0][p->length1 - 1] = 0;
								p->information[1][p->length1 - 1] = 0;
								np_arr[y][x]++;
								np_arr[y][x]++;

								p->length1--;
								p->length2--;
								c2f2->length--;

								for (c2f3 = c2f2->next; c2f3 != p->chain2; c2f3 = c2f3->next)
								{
									c2f3->start--;
								}

								if (c2f2->length == 0)
								{
									(c2f2->prior)->next = c2f2->next;
									(c2f2->next)->prior = c2f2->prior;
									free(c2f2);
								}
							}
						}

						if (p->length1 == 1)
						{
							fresh_unit();
							break;
						}

						while (1)  // Breaking of phosphodiester bonds within an RNA chain
						{
							for (j = p->length1; j > 1; j--)
							{
								f = PBB;
								for (c2f1 = p->chain2->prior; c2f1 != p->chain2; c2f1 = c2f1->prior)
								{
									if (j > c2f1->start + 1)
									{
										if (j <= c2f1->length + c2f1->start)    // Falling into double chain region
										{
											f = f * sqrt(f);
										}
										break;
									}
								}

								c2f1 = p->chain2->prior;
								c2f2 = p->chain2->next;
								if (randd() < f)
								{
									p3 = (struct rna*)malloc(LEN);
									if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }
									memset(p3->information, 0, sizeof(p3->information));
									c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
									if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
									p3->chain2 = c2f;
									c2f->next = p3->chain2;
									c2f->prior = p3->chain2;

									for (b = 0; b < p->length1 - j + 1; b++)
									{
										p3->information[0][b] = p->information[0][b + j - 1];
										p->information[0][b + j - 1] = 0;
									}
									p3->length1 = p->length1 - j + 1;
									p->length1 = j - 1;
									p3->length2 = 0;
									p3->type1 = 0;
									p3->type2 = 0;

									if (c2f1 != p->chain2 && c2f1->start + c2f1->length > j - 1)
									{
										for (b = 0; b < c2f1->start + c2f1->length - j + 1; b++)
										{
											p3->information[1][b] = p->information[1][b + j - 1];
											p->information[1][b + j - 1] = 0;
										}

										while (c2f1 != p->chain2)
										{
											c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
											if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
											c2f3->prior = p3->chain2;
											c2f3->next = p3->chain2->next;
											p3->chain2->next->prior = c2f3;
											p3->chain2->next = c2f3;

											if (j <= c2f1->start + 1)
											{
												c2f3->length = c2f1->length;
												c2f3->start = c2f1->start - j + 1;
												p3->length2 = p3->length2 + c2f3->length;
												p->length2 = p->length2 - c2f3->length;

												(c2f1->prior)->next = c2f1->next;
												(c2f1->next)->prior = c2f1->prior;
												c2f = c2f1;
												c2f1 = c2f1->prior;
												free(c2f);

												if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
												if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p3-chain-length"); exit(0); }

												if (j >= c2f1->start + c2f1->length + 1)
												{
													break;
												}
											}
											else
											{
												c2f3->length = c2f1->start + c2f1->length - j + 1;
												c2f3->start = 0;
												c2f1->length = c2f1->length - c2f3->length;

												p3->length2 = p3->length2 + c2f3->length;
												p->length2 = p->length2 - c2f3->length;

												if (p->length1 < p->length2) { printf("unexpected error on p-chain-length"); exit(0); }
												if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
												if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p-chain-length"); exit(0); }
												break;
											}
										}
									}
									p3->prior = room_head[!h][y][x];
									p3->next = room_head[!h][y][x]->next;
									if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
									room_head[!h][y][x]->next = p3;
									break;
								}
							}
							if (j == 1) break;
						}
					}
					else if (p->type1 == 1)  // Circular chain
					{
						flag5 = 0;     // A flag indicating whether the circle chain breaking has happened
						for (j = p->length1; j > 0; j--)
						{
							f = PBB; flag1 = 0; 
							for (c2f1 = p->chain2->prior; c2f1 != p->chain2; c2f1 = c2f1->prior)
							{
								if (j > c2f1->start + 1)
								{
									if (j <= c2f1->length + c2f1->start)         // Falling into double chain region
									{
										f = f * sqrt(f); flag1 = 1;
									}
									break;
								}
							}
							if (j == 1 && p->type2 == 1)f = f * sqrt(f);

							if (randd() < f)
							{
								flag5 = j;  // The location of the circle chain breaking 
								p->type1 = 0;
								p->type2 = 0;
								if (j == 1)break;

								for (b = 0; b < j - 1; b++)
								{
									temp_information[0][b] = p->information[0][b];
									temp_information[1][b] = p->information[1][b];
								}
								for (b = 0; b < p->length1; b++)
								{
									if (b < p->length1 - j + 1)
									{
										p->information[0][b] = p->information[0][b + j - 1];
										p->information[1][b] = p->information[1][b + j - 1];
									}
									else
									{
										p->information[0][b] = temp_information[0][b - p->length1 + j - 1];
										p->information[1][b] = temp_information[1][b - p->length1 + j - 1];
									}
								}
								if (p->chain2->next == p->chain2)break;

								if (flag1 == 1)
								{
									c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
									if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
									c2f3->prior = c2f1;
									c2f3->next = c2f1->next;
									c2f1->next->prior = c2f3;
									c2f1->next = c2f3;

									c2f3->start = j - 1;
									c2f3->length = c2f1->start + c2f1->length - j + 1;
									c2f1->length = c2f1->length - c2f3->length;
								}
								else
								{
									c2f3 = c2f1->next;
								}

								for (c2f2 = c2f3; c2f2 != p->chain2; c2f2 = c2f2->next) c2f2->start = c2f2->start - j + 1;
								for (c2f2 = p->chain2->next; c2f2 != c2f1->next; c2f2 = c2f2->next) c2f2->start = c2f2->start + p->length1 - j + 1;

								if (c2f3 != p->chain2 && c2f1 != p->chain2)   
								{
									p->chain2->prior->next = p->chain2->next;
									p->chain2->next->prior = p->chain2->prior;
									p->chain2->next = c2f3;
									p->chain2->prior = c2f1;
									c2f1->next = p->chain2;
									c2f3->prior = p->chain2;
								}

								if (p->length1 < p->length2) { printf("unexpected error on p-chain-length5"); exit(0); }

								break;
							}
						}
						if (flag5 != 0)  // Circlar chain breaking has happened
						{
							flag4 = p->length1 - flag5 + 1;
							while (1)
							{
								for (j = p->length1; j > flag4; j--)  // Only consider those bonds that have not been considered concerning the circular chain breaking
								{
									f = PBB;
									for (c2f1 = p->chain2->prior; c2f1 != p->chain2; c2f1 = c2f1->prior)
									{
										if (j > c2f1->start + 1)
										{
											if (j <= c2f1->length + c2f1->start)         // Falling into double chain region
											{
												f = f * sqrt(f);
											}
											break;
										}
									}

									c2f1 = p->chain2->prior;
									c2f2 = p->chain2->next;
									if (randd() < f)
									{
										p3 = (struct rna*)malloc(LEN);
										if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }
										memset(p3->information, 0, sizeof(p3->information));
										c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
										p3->chain2 = c2f;
										c2f->next = p3->chain2;
										c2f->prior = p3->chain2;

										for (b = 0; b < p->length1 - j + 1; b++)
										{
											p3->information[0][b] = p->information[0][b + j - 1];
											p->information[0][b + j - 1] = 0;
										}
										p3->length1 = p->length1 - j + 1;
										p->length1 = j - 1;
										p3->length2 = 0;
										p3->type1 = 0;
										p3->type2 = 0;

										if (c2f1 != p->chain2 && c2f1->start + c2f1->length > j - 1)   
										{
											for (b = 0; b < c2f1->start + c2f1->length - j + 1; b++)
											{
												p3->information[1][b] = p->information[1][b + j - 1];
												p->information[1][b + j - 1] = 0;
											}

											while (c2f1 != p->chain2)
											{
												c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
												if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
												c2f3->prior = p3->chain2;
												c2f3->next = p3->chain2->next;
												p3->chain2->next->prior = c2f3;
												p3->chain2->next = c2f3;

												if (j <= c2f1->start + 1)
												{
													c2f3->length = c2f1->length;
													c2f3->start = c2f1->start - j + 1;
													p3->length2 = p3->length2 + c2f3->length;
													p->length2 = p->length2 - c2f3->length;

													(c2f1->prior)->next = c2f1->next;
													(c2f1->next)->prior = c2f1->prior;
													c2f = c2f1;
													c2f1 = c2f1->prior;
													free(c2f);

													if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
													if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p3-chain-length"); exit(0); }

													if (j >= c2f1->start + c2f1->length + 1)
													{
														break;
													}
												}
												else
												{
													c2f3->length = c2f1->start + c2f1->length - j + 1;
													c2f3->start = 0;
													c2f1->length = c2f1->length - c2f3->length;

													p3->length2 = p3->length2 + c2f3->length;
													p->length2 = p->length2 - c2f3->length;

													if (p->length1 < p->length2) { printf("unexpected error on p-chain-length"); exit(0); }
													if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
													if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p-chain-length"); exit(0); }
													break;
												}
											}
										}
										p3->prior = room_head[!h][y][x];
										p3->next = room_head[!h][y][x]->next;
										if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
										room_head[!h][y][x]->next = p3;
										break;
									}
								}
								if (j == flag4) break;
							}
						}
					}
				}
				fresh_unit();
				break;

			case 2:                         //Template-directed synthesis
				for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next)         // Template-directed ligation
				{
					rtdaddlig = randd();
					c2f3 = c2f2->next;
					if (c2f3 != p->chain2)
					{
						if (c2f2->length + c2f2->start == c2f3->start && rtdaddlig < PTL)
						{
							c2f2->length = c2f2->length + c2f3->length;
							c2f2->next = c2f3->next;
							(c2f3->next)->prior = c2f2;
							free(c2f3);

							c2f2 = c2f2->prior;  continue; // Equal chance for the ligation of every nick
						}
					}
					else
					{
						if (p->type1 == 1 && c2f2->length + c2f2->start == p->length1 && p->chain2->next->start == 0 && rtdaddlig < PTL)
						{
							c2f1 = p->chain2->prior;
							c2f4 = p->chain2->next;
							if (c2f1 == c2f4)
							{
								p->type2 = 1;
								break;
							}

							for (b = 0; b < c2f1->length; b++)
							{
								temp_information[0][b] = p->information[0][b + p->length1 - c2f1->length];
								temp_information[1][b] = p->information[1][b + p->length1 - c2f1->length];
							}
							for (b = p->length1 - 1; b >= 0; b--)
							{
								if (b >= c2f1->length)
								{
									p->information[0][b] = p->information[0][b - c2f1->length];
									p->information[1][b] = p->information[1][b - c2f1->length];
								}
								else
								{
									p->information[0][b] = temp_information[0][b];
									p->information[1][b] = temp_information[1][b];
								}
							}

							c2f4->length = c2f4->length + c2f1->length;   
							for (c2f4 = c2f4->next; c2f4 != p->chain2; c2f4 = c2f4->next)  
							{
								c2f4->start = c2f4->start + c2f1->length;
							}
							c2f1->prior->next = p->chain2;
							p->chain2->prior = c2f1->prior;
							free(c2f1);
							break;
						}
					}
				}

				if (p->length1 < 4) { fresh_unit(); break; }  // Short RNA cannot act as template

				if (p->type1 == 1) f1 = PAT;
				else f1 = PAT * FLT;   // Linear chain as a worse template
				for (p3 = p->next; p3 != p; p3 = p3->next)  // Substrates attraction by an RNA template
				{
					if (p3 == room_head[h][y][x])
					{
						p3 = room_head[h][y][x]->next;
						if (p3 == p)break;
					}
					if (p3->length2 == 0 && p3->length1 <= p->length1)
					{
						c2f2 = p->chain2->next;
						if (p->type1 == 0 && (p->length2 == 0 || c2f2->start != 0))
						{
							if (p->length2 == 0)
							{
								length = p->length1;
							}
							else if (c2f2->start != 0)
							{
								length = c2f2->start;
							}

							if (p3->length1 <= length)
							{
								for (c = 0; c <= length - p3->length1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + b]) == 5)continue; // Base-pairing
										else if (randd() < PFP)continue;   // False base-pairing is allowed to an extent
										else { flag = 1; break; }
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd() * FDA;  // De novo attraction would be more difficult (the primer effect)
									if (rtdaddphili < f1)
									{
										for (a = 0; a < p3->length1; a++)
										{
											p->information[1][c + a] = p3->information[0][p3->length1 - 1 - a];
											if (p->information[1][c + a] == 0) { printf("unexpected error on add-nucleotide3"); exit(0); };
										}

										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->prior = p->chain2;
										c2f3->next = p->chain2->next;
										p->chain2->next->prior = c2f3;
										p->chain2->next = c2f3;

										c2f3->start = c;
										c2f3->length = p3->length1;

										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										break;
									}
								}
							}
						}
						else if (p->type1 == 1 && p->length2 == 0)
						{
							if (p->length1 >= p3->length1)
							{
								for (c = 0; c <= p->length1 - 1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if (c + b < p->length1)
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + b]) == 5)continue;
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
										else
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + b - p->length1]) == 5)continue;
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd() * FDA;
									if (rtdaddphili < f1)
									{
										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->prior = p->chain2;
										c2f3->next = p->chain2->next;
										p->chain2->next->prior = c2f3;
										p->chain2->next = c2f3;
										c2f3->length = p3->length1;

										if (p3->length1 + c <= p->length1)
										{
											for (a = 0; a < p3->length1; a++) {
												p->information[1][c + a] = p3->information[0][p3->length1 - 1 - a];
												if (p->information[1][c + a] == 0) { printf("unexpected error on add-nucleotide41"); exit(0); };
											}

											c2f3->start = c;
										}
										else
										{
											for (b = 0; b < p->length1 - c; b++)
												temp_information[0][b] = p->information[0][b + c];

											for (b = p->length1 - 1; b >= 0; b--)
											{
												if (b >= p->length1 - c) {
													p->information[0][b] = p->information[0][b - p->length1 + c];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide42"); exit(0); };
												}
												else {
													p->information[0][b] = temp_information[0][b];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide43"); exit(0); };
												}
											}
											for (a = 0; a < p3->length1; a++) { p->information[1][a] = p3->information[0][p3->length1 - 1 - a]; }

											c2f3->start = 0;
										}
										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										break;
									}
								}
							}
						}
						else if (p->type1 == 1 && c2f2->start != 0)
						{
							c2f1 = p->chain2->prior;
							length = c2f2->start + p->length1 - c2f1->start - c2f1->length;

							if (length >= p3->length1)
							{
								for (c = 0; c <= length - p3->length1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if (c + c2f1->start + c2f1->length + b < p->length1)
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + c2f1->start + c2f1->length + b]) == 5)continue;
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
										else
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + c2f1->start + c2f1->length + b - p->length1]) == 5)continue;
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd();
									if (c != 0)rtdaddphili = rtdaddphili * FDA;
									if (rtdaddphili < f1)
									{
										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->length = p3->length1;

										if (p3->length1 + c <= p->length1 - c2f1->start - c2f1->length)
										{
											for (a = 0; a < p3->length1; a++)
												p->information[1][c2f1->start + c2f1->length + c + a] = p3->information[0][p3->length1 - 1 - a];

											c2f3->next = p->chain2;
											c2f3->prior = p->chain2->prior;
											p->chain2->prior->next = c2f3;
											p->chain2->prior = c2f3;
											c2f3->start = c + c2f1->start + c2f1->length;
										}
										else if (c >= p->length1 - c2f1->start - c2f1->length)
										{
											for (a = 0; a < p3->length1; a++)
												p->information[1][c - p->length1 + c2f1->start + c2f1->length + a] = p3->information[0][p3->length1 - 1 - a];

											c2f3->prior = p->chain2;
											c2f3->next = p->chain2->next;
											p->chain2->next->prior = c2f3;
											p->chain2->next = c2f3;
											c2f3->start = c - p->length1 + c2f1->start + c2f1->length;
										}
										else
										{
											for (b = 0; b < p->length1 - c2f1->start - c2f1->length - c; b++)
											{
												temp_information[0][b] = p->information[0][b + c2f1->start + c2f1->length + c];
												temp_information[1][b] = p->information[1][b + c2f1->start + c2f1->length + c];
											}
											for (b = p->length1 - 1; b >= 0; b--)
											{
												if (b >= p->length1 - c2f1->start - c2f1->length - c)
												{
													p->information[0][b] = p->information[0][b - p->length1 + c2f1->start + c2f1->length + c];
													p->information[1][b] = p->information[1][b - p->length1 + c2f1->start + c2f1->length + c];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide51"); exit(0); }
												}
												else
												{
													p->information[0][b] = temp_information[0][b];
													p->information[1][b] = temp_information[1][b];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide52"); exit(0); };
												}
											}
											for (a = 0; a < p3->length1; a++) { p->information[1][a] = p3->information[0][p3->length1 - 1 - a]; }

											for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next)
												c2f2->start = c2f2->start + p->length1 - c2f1->start - c2f1->length - c;

											c2f3->prior = p->chain2;
											c2f3->next = p->chain2->next;
											p->chain2->next->prior = c2f3;
											p->chain2->next = c2f3;
											c2f3->start = 0;
										}
										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										break;
									}
								}
							}
						}
						flag3 = 0;  // A flag indicating whether the attraction has happened
						for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next)
						{
							if (c2f2->next == p->chain2)
							{
								if (p->type1 == 1 && p->chain2->next->start != 0) break;
								length = p->length1 - c2f2->start - c2f2->length;
							}
							else
							{
								length = c2f2->next->start - c2f2->start - c2f2->length;
							}

							if (p3->length1 <= length)
							{
								for (c = 0; c <= length - p3->length1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c2f2->start + c2f2->length + c + b]) == 5)continue;
										else if (randd() < PFP)continue;
										else { flag = 1; break; }
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd();
									if (c != 0)rtdaddphili = rtdaddphili * FDA;
									if (rtdaddphili < f1)
									{
										for (a = 0; a < p3->length1; a++) {
											p->information[1][c2f2->start + c2f2->length + c + a] = p3->information[0][p3->length1 - 1 - a];
											if (p->information[1][c2f2->start + c2f2->length + c + a] == 0) { printf("unexpected error on add-nucleotide6"); exit(0); };
										}

										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->prior = c2f2;
										c2f3->next = c2f2->next;
										c2f2->next->prior = c2f3;
										c2f2->next = c2f3;

										c2f3->start = c2f2->start + c2f2->length + c;
										c2f3->length = p3->length1;

										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										flag3 = 1;   // The attraction has happened.
										break;
									}
								}
							}
						}
						if (flag3 == 1) break;
					}
				}
				fresh_unit();
				break;

			case 3:                     // Separation of RNA duplex
				if (p->length2 != 0)
				{
					j = 0;
					for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next) j++;
					randseq = randl(j);
					c2f2 = p->chain2->next;
					for (j = 1; j <= randseq; j++)
					{
						c2f2 = c2f2->next;
					}

					if (randd() < pow(PSP, sqrt(c2f2->length * 1.0)))
					{
						p3 = (struct rna*)malloc(LEN);
						if (!p3) { printf("\t%dsep--memeout\n", i); exit(0); }
						memset(p3->information, 0, sizeof(p3->information));
						c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
						if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
						p3->chain2 = c2f;
						c2f->next = p3->chain2;
						c2f->prior = p3->chain2;

						for (b = 0; b < c2f2->length; b++)
						{
							p3->information[0][b] = p->information[1][c2f2->start + c2f2->length - 1 - b];
							if (p3->information[0][b] == 0) { printf("unexpected error on seperation"); exit(0); }
							p->information[1][c2f2->start + c2f2->length - 1 - b] = 0;
						}
						p3->length1 = c2f2->length;
						p->length2 = p->length2 - c2f2->length;
						p3->length2 = 0;
						if (p->type2 == 1)
						{
							p3->type1 = 1;
							p->type2 = 0;
						}
						else p3->type1 = 0;
						p3->type2 = 0;

						(c2f2->prior)->next = c2f2->next;
						(c2f2->next)->prior = c2f2->prior;
						free(c2f2);

						p3->prior = room_head[!h][y][x];
						p3->next = room_head[!h][y][x]->next;
						if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
						room_head[!h][y][x]->next = p3;
					}
				}
				fresh_unit();
				break;

			case 4:    // Functioning of ER
				if (p->type1 == 0 && p->length2 == 0 && p->length1 >= erlength && p->length1 < 1.5 * erlength)  
				{																// A ribozyme cannot be much longer than its characteristic domain
					flag = findseq(erseq, erlength, p);
					if (flag == 0)
					{
						er_turn = TER;    // The times of functioning for an ER

						for (p3 = p->next; p3 != p; p3 = p3->next)
						{
							if (p3 == room_head[h][y][x])
							{
								p3 = room_head[h][y][x]->next;
								if (p3 == p)break;
							}
							if (p3->length1 == 1 || p3->type1 == 1) continue;
							if (findseq(erseq, erlength, p3) == 0 && p3->length2 == 0 && p3->length1 < 1.5 * erlength)
								if (randd() < PRP)continue;  // A ribozyme will be protected to an extent (from the cleaving of ER) because of its folding
							if (findseq(repseq, replength, p3) == 0 && p3->type1 == 0 && p3->length2 == 0 && p3->length1 < 1.5 * replength)
								if (randd() < PRP)continue;  // A ribozyme will be protected because of its folding

							if (p3->information[1][p3->length1 - 1] == 0)
							{
								if (er_turn <= 0) break;  // Because the ER is used out
								er_turn--;
								if (randd() < PBER)
								{
									p4 = (struct rna*)malloc(LEN);
									if (!p4) { printf("\ter form_monomer--memeout\n"); exit(0); }
									memset(p4->information, 0, sizeof(p4->information));
									c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
									if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
									p4->chain2 = c2f;
									c2f->next = p4->chain2;
									c2f->prior = p4->chain2;

									p4->information[0][0] = p3->information[0][p3->length1 - 1];
									p4->length1 = 1;
									p4->length2 = 0;
									p4->type1 = 0;
									p4->type2 = 0;
									p4->prior = room_head[!h][y][x];
									p4->next = room_head[!h][y][x]->next;
									if (p4->next != room_head[!h][y][x])(p4->next)->prior = p4;
									room_head[!h][y][x]->next = p4;

									p3->information[0][p3->length1 - 1] = 0;
									p3->length1--;
								}
							}
						}
					}
				}
				if (p->type1 == 0 && p->length2 == 0 && p->length1 >= replength && p->length1 < 1.5 * replength)  //:A ribozyme cannot be much longer than its characteristic domain
				{
					flag = findseq(repseq, replength, p);    // REP functioning
					if (flag == 0)
					{
						rep_turn = TREP;  // The times of functioning
						flag3 = 0;
						for (p3 = p->next; p3 != p; p3 = p3->next)
						{
							if (flag3 == 1) break;   // Because REP is used out
							if (p3 == room_head[h][y][x])
							{
								p3 = room_head[h][y][x]->next;
								if (p3 == p)break;
							}
							if (p3->length2 == 0 || (p3->chain2->next->next == p3->chain2 && p3->type1 != 1))continue;

							for (c2f2 = p3->chain2->next; c2f2 != p3->chain2; c2f2 = c2f2->next)
							{
								c2f3 = c2f2->next;
								if (c2f3 != p3->chain2)
								{
									if (c2f2->start + c2f2->length == c2f3->start)
									{
										if (rep_turn <= 0) { flag3 = 1; break; }  // REP is used out
										rep_turn--;
										if (randd() < PTLR)
										{
											c2f2->length = c2f2->length + c2f3->length;
											c2f2->next = c2f3->next;
											(c2f3->next)->prior = c2f2;
											free(c2f3);

											c2f2 = c2f2->prior;  continue; // Equal chance for the ligation of every nick
										}
									}
								}
								else
								{
									if (p3->type1 == 1 && c2f2->length + c2f2->start == p3->length1 && p3->chain2->next->start == 0)
									{
										if (rep_turn <= 0) { flag3 = 1; break; }   // Rep is used out
										rep_turn--;
										if (randd() < PTLR)
										{
											c2f1 = p3->chain2->prior;
											c2f4 = p3->chain2->next;
											if (c2f1 == c2f4)
											{
												p3->type2 = 1;
												break;     // Double circular chain without nick
											}

											for (b = 0; b < c2f1->length; b++)
											{
												temp_information[0][b] = p3->information[0][b + p3->length1 - c2f1->length];
												temp_information[1][b] = p3->information[1][b + p3->length1 - c2f1->length];
											}
											for (b = p3->length1 - 1; b >= 0; b--)
											{
												if (b >= c2f1->length)
												{
													p3->information[0][b] = p3->information[0][b - c2f1->length];
													p3->information[1][b] = p3->information[1][b - c2f1->length];
												}
												else
												{
													p3->information[0][b] = temp_information[0][b];
													p3->information[1][b] = temp_information[1][b];
												}
											}

											c2f4->length = c2f4->length + c2f1->length;
											for (c2f4 = c2f4->next; c2f4 != p3->chain2; c2f4 = c2f4->next)
											{
												c2f4->start = c2f4->start + c2f1->length;
											}
											c2f1->prior->next = p3->chain2;
											p3->chain2->prior = c2f1->prior;
											free(c2f1);
											break;
										}
									}
								}
							}
						}
					}
				}
				fresh_unit();
				break;

			case 5:            // RNA moving
				if (randd() * RMRW < PMN)    // The movement of an RNA molecule is related to its mass (the Zimm model)
				{
					randcase1 = randl(4);   // Four possible directions
					switch (randcase1)
					{
					case 0:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][y][(N + x - 1) % N]->next;   // Toroidal topology
						room_head[!h][y][(N + x - 1) % N]->next = p;
						p->next = p3;
						p->prior = room_head[!h][y][(N + x - 1) % N];
						if (p3 != room_head[!h][y][(N + x - 1) % N])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					case 1:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][y][(x + 1) % N]->next;
						room_head[!h][y][(x + 1) % N]->next = p;
						p->next = p3;
						p->prior = room_head[!h][y][(x + 1) % N];
						if (p3 != room_head[!h][y][(x + 1) % N])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					case 2:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][(N + y - 1) % N][x]->next;
						room_head[!h][(N + y - 1) % N][x]->next = p;
						p->next = p3;
						p->prior = room_head[!h][(N + y - 1) % N][x];
						if (p3 != room_head[!h][(N + y - 1) % N][x])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					case 3:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][(y + 1) % N][x]->next;
						room_head[!h][(y + 1) % N][x]->next = p;
						p->next = p3;
						p->prior = room_head[!h][(y + 1) % N][x];
						if (p3 != room_head[!h][(y + 1) % N][x])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					default:printf("rna moving error");
					}
				}
				else fresh_unit();
				break;

			default: printf("rna case error");
			}
		}
	}
}

//------------------------------------------------------------------------------
void record(void)               // Data recording
{
	int ch_num[MAX_RNA_LENGTH], ch_c_num[MAX_RNA_LENGTH], ch_er_num[MAX_RNA_LENGTH], ch_rep_num[MAX_RNA_LENGTH], ch_errep_num[MAX_RNA_LENGTH], ch_0_num[MAX_RNA_LENGTH], long_chain_num, si;

	FILE* fptxt, * fptxt1;
	errno_t err, err1;
	err = fopen_s(&fptxt, "monitor.txt", "at");
	if (err != 0) { printf("cannot open file");  exit(-1); }
	err1 = fopen_s(&fptxt1, "picture.txt", "at");
	if (err1 != 0) { printf("cannot open file1");  exit(-1); }

	er[g] = 0;
	er_rib[g] = 0;
	cir_er[g] = 0;

	rep[g] = 0;
	rep_rib[g] = 0;
	cir_rep[g] = 0;

	ct[g] = 0;
	cir_ct[g] = 0;

	ernc[g] = 0;
	cir_ernc[g] = 0;

	ctnc[g] = 0;
	cir_ctnc[g] = 0;

	errep[g] = 0;
	cir_errep[g] = 0;

	unit[g] = 0;
	cir_unit[g] = 0;
	total_mat_num[g] = 0;
	np_num[g] = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			np_num[g] += np_arr[y][x];
			for (p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)
			{
				total_mat_num[g] += p->length1 + p->length2;
				if (p->type1 == 0) unit[g]++;
				else cir_unit[g]++;

				flag1 = findseq(erseq, erlength, p);
				flag2 = findseq(repseq, replength, p);
				if (flag1 == 0 && flag2 != 0)
				{
					if (p->type1 == 0)
					{
						er[g]++;
						if (p->length1 < 1.5 * erlength)er_rib[g]++;
					}
					else cir_er[g]++;
				}
				else if (flag1 != 0 && flag2 == 0)
				{
					if (p->type1 == 0)
					{
						rep[g]++;
						if (p->length1 < 1.5 * replength)rep_rib[g]++;
					}
					else cir_rep[g]++;
				}
				flag = findseq(ctseq, ctlength, p);
				if (flag == 0)
				{
					if (p->type1 == 0) ct[g]++;
					else cir_ct[g]++;
				}

				flag = findseq(erncseq, ernclength, p);
				if (flag == 0 )
				{
					if (p->type1 == 0)ernc[g]++;
					else cir_ernc[g]++;
				}
				flag = findseq(ctncseq, ctnclength, p);
				if (flag == 0 )
				{
					if (p->type1 == 0)ctnc[g]++;
					else cir_ctnc[g]++;
				}
				flag = findseq(errepseq, erreplength, p);
				if (flag == 0 )
				{
					if (p->type1 == 0)errep[g]++;
					else cir_errep[g]++;
				}
			}
		}
	}
	total_mat_num[g] += np_num[g];

	printf("----- step=%d: er=%d-%d(%d), rep=%d-%d(%d), ct=%d(%d), ernc=%d(%d), ctnc=%d(%d), errep=%d(%d), unit=%d (%d)  (tm=%d, np=%d)\n", i,
		(int)er[g], (int)er_rib[g], (int)cir_er[g],
		(int)rep[g], (int)rep_rib[g], (int)cir_rep[g],
		(int)ct[g], (int)cir_ct[g],
		(int)ernc[g], (int)cir_ernc[g],
		(int)ctnc[g], (int)cir_ctnc[g],
		(int)errep[g], (int)cir_errep[g],
		(int)unit[g], (int)cir_unit[g], (int)total_mat_num[g], (int)np_num[g]);
	fprintf(fptxt, "----- step=%d: er=%d-%d(%d), rep=%d-%d(%d), ct=%d(%d), ernc=%d(%d), ctnc=%d(%d), errep=%d(%d), unit=%d (%d)  (tm=%d, np=%d)\n", i,
		(int)er[g], (int)er_rib[g], (int)cir_er[g],
		(int)rep[g], (int)rep_rib[g], (int)cir_rep[g],
		(int)ct[g], (int)cir_ct[g],
		(int)ernc[g], (int)cir_ernc[g],
		(int)ctnc[g], (int)cir_ctnc[g],
		(int)errep[g], (int)cir_errep[g],
		(int)unit[g], (int)cir_unit[g], (int)total_mat_num[g], (int)np_num[g]);
	fprintf(fptxt1, "step=%d: er=%d-%d(%d), rep=%d-%d(%d), ct=%d(%d), ernc=%d(%d), ctnc=%d(%d), errep=%d(%d), unit=%d (%d)  (tm=%d, np=%d)\n", i,
		(int)er[g], (int)er_rib[g], (int)cir_er[g],
		(int)rep[g], (int)rep_rib[g], (int)cir_rep[g],
		(int)ct[g], (int)cir_ct[g],
		(int)ernc[g], (int)cir_ernc[g],
		(int)ctnc[g], (int)cir_ctnc[g],
		(int)errep[g], (int)cir_errep[g],
		(int)unit[g], (int)cir_unit[g], (int)total_mat_num[g], (int)np_num[g]);

	for (si = 0; si < LONG_CHAIN_LEN; si++)    // For monitoring the evolution of RNAs' chain lengths
	{
		ch_num[si] = 0;  ch_c_num[si] = 0;

		ch_er_num[si] = 0;

		ch_0_num[si] = 0;

		ch_rep_num[si] = 0;

		ch_errep_num[si] = 0;
	}
	long_chain_num = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				ch_num[p->length1 - 1]++;
				if (p->type1 == 1)ch_c_num[p->length1 - 1]++;
				if (p->length1 > LONG_CHAIN_LEN)
				{
					long_chain_num++;
					for (int t = 0; t < p->length1; t++)    // The sequences of long RNA chains are shown
					{
						switch (p->information[0][t])
						{
						case 1: printf("A"); fprintf(fptxt, "A"); break;
						case 2: printf("C"); fprintf(fptxt, "C"); break;
						case 3: printf("G"); fprintf(fptxt, "G"); break;
						case 4: printf("U"); fprintf(fptxt, "U"); break;
						default: printf("error");
						}
					}
					printf("\t%d\t%d\n", p->length1, p->type1);
					fprintf(fptxt, "\t%d\t%d\n", p->length1, p->type1);
					for (int t = 0; t < p->length1; t++)
					{
						switch (p->information[1][t])
						{
						case 0: printf("-"); fprintf(fptxt, "-"); break;
						case 1: printf("A"); fprintf(fptxt, "A"); break;
						case 2: printf("C"); fprintf(fptxt, "C"); break;
						case 3: printf("G"); fprintf(fptxt, "G"); break;
						case 4: printf("U"); fprintf(fptxt, "U"); break;
						default: printf("error");
						}
					}
					printf("\t\t%d\n", p->type2);
					fprintf(fptxt, "\t\t%d\n", p->type2);
				}

				flag1 = findseq(erseq, erlength, p);
				flag2 = findseq(repseq, replength, p);

				if (flag1 == 0 && flag2 != 0)ch_er_num[p->length1 - 1]++;
				else if (flag1 != 0 && flag2 == 0)ch_rep_num[p->length1 - 1]++;
				else if (flag1 == 0 && flag2 == 0)ch_errep_num[p->length1 - 1]++;
				else ch_0_num[p->length1 - 1]++;
			}
		}
	}

	for (si = 0; si < LONG_CHAIN_LEN; si++)
	{
		printf("%dnt-%d(%d)(%d|%d|%d|%d), ", si + 1, ch_num[si], ch_c_num[si], ch_er_num[si], ch_rep_num[si], ch_errep_num[si], ch_0_num[si]);
		fprintf(fptxt, "%dnt-%d(%d)(%d|%d|%d|%d), ", si + 1, ch_num[si], ch_c_num[si], ch_er_num[si], ch_rep_num[si], ch_errep_num[si], ch_0_num[si]);
	}
	printf("\nChains over %d = %d   step=%d\n\n", LONG_CHAIN_LEN, long_chain_num, i);
	fprintf(fptxt, "\nChains over %d = %d   step=%d\n\n", LONG_CHAIN_LEN, long_chain_num, i);

	g++;

	fclose(fptxt);
	fclose(fptxt1);
}

//------------------------------------------------------------------------------
void freepool(void)        // Memory releasing  
{
	int m;
	for (m = 0; m < 2; m++)
	{
		for (y = 0; y < N; y++)
		{
			for (x = 0; x < N; x++)
			{
				while (1)
				{
					if (room_head[m][y][x]->next != room_head[m][y][x])
					{
						p = room_head[m][y][x]->next;
						room_head[m][y][x]->next = p->next;
						while (1)
						{
							if (p->chain2->next != p->chain2)
							{
								c2f = p->chain2->next;
								p->chain2->next = c2f->next;
								free(c2f);
							}
							else break;
						}
						free(p->chain2);
						free(p);
					}
					else break;
				}
				free(room_head[m][y][x]);
			}
		}
	}
}

//------------------------------------------------------------------------------ 
int main()
{
	inits();        // Initialization of the system
	for (i = 0; i <= STEPNUM; i++)      // Monte-Carlo cycle
	{
		if (i == 10000)inoculate(1,INOCUNUM); //Inoculating RNA species
		if (i >= STAREC && i % RECINT == 0)
		{
			record();  // Data recording
		}
		unit_case();  // Action of units (molecules) in a time step
		h = !h;       // Information transfering for the next step
	}
	freepool();
	return (0);
}
//========================================================  End of the program

