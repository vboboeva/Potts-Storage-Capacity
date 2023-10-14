#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "functions.cpp"
#include "generate_pattern_v4.cpp"
//#include "generate_pattern_rand.cpp"
#include <mpi.h>
//#include <igraph/igraph.h>
//#include "/usr/local/include/igraph/igraph.h"

/*choose range and step of pattern storage and retrieval*/
int Pmin = 10;
int Pstep = 10;
int p;

/*choose how many runs with different sets of patterns for quenched avg*/
int number_different_set_patterns = 1;

/*NET PARAMETERS*/
int N = 2000;
int Cm = 200;
int S = 5;
double U = 0.5;
double nu = 0.0;
double a = 0.1;

/*pattern parameters*/
int tot_num_trials = 1;
int NumSet = 10;

/*latching parameters*/
double b1 = 0.5;
double b2 = 0.0;
double b3 = 0.0;
double w = 0.0;
double g = 10.0;

/*network update*/
int Tnet = 20; 
int time_toprint = 500;

/*type of connectivity*/
int full_connectivity = 0;
int random_dilution = 1;
int variable_connectivity = 0;
int symmetric_random_dilution = 0;
int do_avg_Cij = 0;

/*cue parameters and type*/
double fraction_units_cued = 1.0;
int exp_cue = 0;
int tau = 2*N;
int full_cue = 1;
int partial_cue = 0;

/*graded or discrete update*/
int update_h = 1;
int update_r = 0;
int discrete = 0;
int graded = 1;
int beta = 200;

/*pattern correlations*/
int p_fact;
int Num_fact = 150;
double a_pf = 0.4;
double fact_eigen_slope;
double eps = 0.000001;
double f = 0.05;

/* GLOBAL VARIABLES */
int		**xi; //pattern p*N
int		*number_connections;
double	**s; // state of network N*(S+1)
double	*sold; // previous state of network N*(S+1)
double	****J; // connections with covariance learning rule N*Cm*S*S
double	*units_receive_cue; //  to construct a cue correlated with a pattern
double	**sparsity_by_state; 


double	*m; // overlap of network with each of p patterns
double	*m_sorted; // overlap of network with each of p patterns SORTED
int 	*indices;

double  m_retr; // overlap of pattern retrieved
int     retr; // pattern retrieved
int		retr_f;
int	*fraction_fact_each;

double	m_cue; // overlap of pattern cued
int		cue; // pattern cued

double	*mf; // overlap of network with each of p PARENTS
double	*mf_sorted; // overlap of network with each of p PARENTS
int 	*indices_f;


double	**h; // field that each active state recieves due to connections during learning
double	**r; // field that each active state recieves due to connections during learning
double	**theta; // variable threshold
int 	**Permut; // sequences of random numbers without repetition used for updating
int 	**C; // a nonsymmetric dilution matrix
double 	**C_sym; // a symmetric dilution matrix
double	***fin_state; // the final state of the network needed for computing C0, C1, C2 and q
double	**C0; // fraction of units co-inactive
double	**C1; // fraction of units co-same-state-active
double	**C2; // fraction of units co-active
double	H_interaction, H_quad; // each of the energy terms

int **Factors;
int **PC;
int **CP;

double fraction_retr[10];
double fraction_corr[10];
double fraction_fact[10];

double info;
double info2;
double sparsity_net;
/* FILES */
FILE *storage_cap;
FILE *overlap_intime;
FILE *fin_sparsity;
FILE *fin_overlaps;
FILE *cij;
FILE *q_file;
FILE *state_file;
FILE *C0_file;
FILE *C1_file;
FILE *C2_file;
FILE *fin_state_file;
FILE *H_quad_file;
FILE *H_interaction_file;
FILE *field_file;
FILE *fileconfigs;


/*from functions*/

//extern void read_pattern(int); // reads pattern from file generated by genero_pattern
extern void initializing(); // initializes the network state
extern void construct_CJ(); // constructs J and C
extern void update_state(); // takes trial, unit index to update and counter n as input

extern void getmemory(); 
extern void deletememory(); 
extern void SetUpTables(); // sets up update tables

extern int compare();
extern double compute_info();
extern void compute_m(); // computes overlap of network state with each pattern
extern void compute_m_sorted(); // computes overlap of network state with each pattern and then sorts
extern void compute_m_cue(); // computes overlap of network state with pattern cued
extern void compute_C0_C1_C2();
extern void compute_sparsity_by_state();

extern void compute_m_parents();
extern void compute_m_parents_sorted();
extern void test_m_mf();
extern int compare_mf();


/*functions generate pattern v2*/

extern void GetMemory();
extern void SetFactors();
extern void AssignChildren();
extern void SetPatterns();
extern void SavePatterns();

/*functions generate pattern rand*/

extern void GenerateRandomPatterns();
extern void SaveRandomPatterns();

/*---------------------------RUNS FOR mu_f PATTERNS-----------------------*/

int run_net(int mu_f)
{
	int i, n, trial, x, mu, iii, ttt;
	double I, sp;
	
	for(i=0; i<10; i++)
	{
		fraction_retr[i] = 0;
		fraction_corr[i] = 0;
		fraction_fact[i] = 0;
	}
	info = 0.0;
	info2 = 0.0;
	sparsity_net = 0.0;

	for(i=0; i<Num_fact; i++)
	{
	  fraction_fact_each[i]=0;
	}
	
	/* cueing each pattern*/
	/* for controlling the size of the data set*/
	for(mu=0; mu<mu_f; mu++)
	{
	  /* each trial corresponds to retrieving the same pattern, retr, corrupted differently */
	  for(trial=0;trial<tot_num_trials;trial++)
	  {
        cue=mu; // retrieve this pattern
	    retr=-1;
	    retr_f=-1;
	    m_retr=-100;
	    m_cue=-100;
		n=0; // unit update counter
		
		//printf("retrieval of pattern %d\n", retr);
		/* initialize network with cue */
		initializing();
        //compute_m_cue(a);
        //printf("m_cue=%.2f\n", m_cue);
		//printf("after initializing\n");
		/* update network */
		for(ttt=0;ttt<Tnet;ttt++)
		{
		  /*x is an updating sequence*/
		  x=(int)(NumSet*drand48()); 
		  for(iii=0;iii<N;iii++)
		  {
			i=Permut[iii][x]; // for asynchronous updating without repetition!
			update_state(trial, i, n); 
		 	//if((n%time_toprint)==0)
		 	//{
		 	  //t=(double)n/N; // effective time
		 	  //compute_m_cue(a);
		 	  //printf("m_cue=%.2f\n", m_cue);
		 	//}
			n++;
		  }
		}
		/* end of update after cue */

		//compute_m_cue(); // computes overlap with cued pattern
		compute_m(); // computes overlap with all patterns
		compute_m_sorted();

		compute_m_parents();
		compute_m_parents_sorted();

		test_m_mf();
		sp = compute_sparsity_net();
		I = compute_info();

		printf( "S=%d p=%d U=%.2f a=%.2f cue=%d retr=%d retr_f=%d m_cue=%.2f m_retr=%.2f I=%.3f entropy=%.3f sp=%.3f\n ", S, p, U, a, cue, retr, retr_f, m_cue, m_retr, I, ((1-a)*log(1/(1-a))+a*log(S/a))/(log(2)), sp);
		printf( "----------------------------------------------------\n");
		
		print_final_state();
		
		if (cue == retr and retr_f < 0) {fraction_retr[0]=fraction_retr[0]+m_retr;}
		else if (cue != retr and retr_f < 0) {fraction_corr[0]=fraction_corr[0]+m_retr;}
		else if (retr_f >= 0) {fraction_fact[0]=fraction_fact[0]+m_retr;}
		
		double thresh;
		for(i=1; i<10; i++)
		{
			thresh = i*0.1;
			//printf("%f\n",thresh);
			if (cue == retr and retr_f < 0 and m_retr >= thresh) {fraction_retr[i]++;}
		
			else if (cue != retr and retr_f < 0 and m_retr >= thresh) {fraction_corr[i]++;}
		
			else if (retr_f >= 0 and m_retr >= thresh) {fraction_fact[i]++;}
		}

		info+=I;
		info2+=I*I;
		
		sparsity_net+=sp;

	  }

	}
return 0;
}

/*---------------------------- START SIMULATION----------------------*/

int main(int argc, char **argv)
{

MPI_Init(NULL, NULL);

fact_eigen_slope = atof(argv[1]);
char *outdir;
outdir = new char [100*sizeof(char)];
sprintf(outdir, "mkdir -p zeta%.6f", fact_eigen_slope); // the -p flag means it will create dir only if it doesn't already exist
system(outdir);

//srand48(time(0));
srand48(1987);

/*Initialize the MPI environment. The two arguments to MPI Init are not
currently used by MPI implementations, but are there in case future
implementations might need the arguments. */
int file_free = 0;

/*Get the number of processes*/
int np;
MPI_Comm_size(MPI_COMM_WORLD, &np);
//printf("np=%d\n", np);

/*Get the rank of the process*/
int proc_id;
MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
MPI_Status status;
if ( proc_id == 0 ) file_free = 1;
else MPI_Recv(&file_free, 1, MPI_INT, proc_id-1, 1, MPI_COMM_WORLD, &status);
if (file_free == 1) // this process calls the routine to read and process input file
/*give read file permission to the next process*/
if (proc_id != np-1) MPI_Send(&file_free, 1, MPI_INT, proc_id+1, 1, MPI_COMM_WORLD);


/*p related things*/
p=Pstep*proc_id+Pmin;

int mu_f;
mu_f = p;
p_fact = f*p;

getmemory();

/*--------------------- Generate correlated patterns ------------------*/
GetMemory();
SetFactors();
AssignChildren();
SetPatterns();
SavePatterns();

/*--------------------- Generate random patterns ----------------------*/
//GenerateRandomPatterns();
//SaveRandomPatterns();

//compute_sparsity_by_state();
SetUpTables();

construct_CJ();

FILE *filefconfigs;
char *buf; 
buf = new char [1000*sizeof(char)];
sprintf(buf, "zeta%.6f/fconfigs%d", fact_eigen_slope, p);
fileconfigs = fopen(buf, "w");

run_net(mu_f);

FILE *fileretr;
char *buffer; 
buffer = new char [1000*sizeof(char)];
sprintf(buffer, "zeta%.6f/p%d", fact_eigen_slope, p);
fileretr = fopen(buffer, "w");

fprintf(fileretr, "%d\t", p);

int i;
for(i=0; i<10; i++)
{
	fprintf(fileretr, "%f\t", fraction_retr[i]/(float)(mu_f*tot_num_trials));
}
for(i=0; i<10; i++)
{
	fprintf(fileretr, "%f\t", fraction_corr[i]/(float)(mu_f*tot_num_trials));
}
for(i=0; i<10; i++)
{
	fprintf(fileretr, "%f\t", fraction_fact[i]/(float)(mu_f*tot_num_trials));
}

fprintf(fileretr, "%f\t%f\t%f\n", info/(float)(mu_f*tot_num_trials), info2/(float)(mu_f*tot_num_trials), sparsity_net/(float)(mu_f*tot_num_trials));	


FILE *filepars;
char *buffer1; 
buffer1 = new char [1000*sizeof(char)];
sprintf(buffer1, "zeta%.6f/factors%d", fact_eigen_slope, p);
filepars = fopen(buffer1, "w");

fprintf(filepars, "%d\t", p);

for(i=0; i<Num_fact; i++)
{
	fprintf(filepars, "%d\t", fraction_fact_each[i]);
}

fclose(fileretr);
fclose(filepars);
fclose(fileconfigs);
// Finalize the MPI environment. No more MPI calls can be made after this

MPI_Finalize();
return 0;

}
