#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
/*----------------------------------------------------------------------*
*		New pattern generator written by Vizhe, version 4,				*
* 		that defines factors as random patterns acting on random		*
* 		subsets of patterns	and not units, as does Ale (version 0).		*
* 		The objective is to reverse the									*
* 		correlation statistics of Ale's algorithm, that is,				*
* 		skewed correlations distribution between patterns, 			 	*
* 		and symmetric between units. Difference with version 1: 		*
* 		parents are random patterns that send influence in a 			*
* 		different Potts state for each unit of a pattern.				*
* 		Difference with version 3: random input included to solve		*
* 		small sparsity for small values of a_pf.						*	
*		Limit cases:													*
* 		1) a_pf ---> 0 gives us randomly correlated patterns			*
* 		2) f*Num_fact = 1 gives us ultrametric patterns	(f=p_fact/p)	*
* ----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
*																		*
*							Enter Parameters							*
* 																		*		
* ----------------------------------------------------------------------*/

/*network parameters*/
int N = 1000;
int S = 5;
int p = 1000;

/*pattern correlations*/
int p_fact = 10; 
int Num_fact = 150;

double eps = 0.000001;

/*correlation parameters*/
double a = 0.3;
double fact_eigen_slope;
double a_pf = 0.4;

int **Factors;
int **Patt;
int **PC;
double **h;
double *hmax;
int *smax;

FILE *fpattern;
FILE *fparents;


/*----------------------------------------------------------------------*
*																		*
*							FUNCTIONS									*
* 																		*		
* ----------------------------------------------------------------------*/

void GetMemory(double a)
{
int i, mu;
	
Factors = new int*[Num_fact];
for(i=0; i<Num_fact; i++)
  Factors[i] = new int[N];
  
Patt = new int*[p];
for(mu=0; mu<p; mu++)
  Patt[mu] = new int[N]; 
  
PC = new int*[Num_fact];
for(i=0; i<Num_fact; i++)
  PC[i] = new int[p_fact];
 
h = new double*[N];
for(i=0; i<N; i++)
  h[i] = new double[S+1];

hmax = new double[N];
smax = new int[N];
  
}

typedef struct {
int index;
double field;
int state;
} mystruct;
  

int cmpf (const void *x, const void *y)
{
  double xx = (*(mystruct*)x).field;
  double yy = (*(mystruct*)y).field;
  /*sort by decreasing order*/
  if (xx > yy) return -1;
  if (xx < yy) return  1;
  return 0;
}

/*Factors are just random patterns*/
void SetFactors(double a)
{
  //char *buffer; 
  //buffer = new char [1000*sizeof(char)];
  //sprintf(buffer, "parents_S%d_a%.1f_apf%.2f_pfact%d_p%d_Numfact%d", S, a, a_pf, p_fact, p, Num_fact);
  //fparents = fopen(buffer, "w");

  mystruct F[N];
  
  double rand_num;
  int rand_state;
	
  int n, i;
  for(n=0;n<Num_fact;n++)
  {
    for(i=0;i<N;i++)
    {
      //rand_num = drand48();
      //rand_state = (int)((double)(S)*drand48());
      //if (rand_num < a) Factors[n][i] = rand_state;
      //else Factors[n][i] = S;
      /*first initialize each pattern to null state*/
      Factors[n][i] = S;		
      
      F[i].index = i;
      F[i].field = drand48();
      F[i].state = (int)((double)(S)*drand48());
    }
    /*sort the units in order of decreasing field keeping the indices: this is done to make patterns with fixed sparsity*/
    qsort(F, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
    for(i=0;i<N;i++)
    {
	  Factors[n][F[i].index] = F[i].state; 
    }  
    //for(i=0;i<N;i++)
    //{
	  //fprintf(fparents,"%d \t", Factors[n][i]); 
    //}
    //fprintf(fparents, "\n"); 
  }
  //fclose(fparents);
}

void AssignChildren(double a)
{
  int n, m, i, patt_picked;
  double rand_num;

  for(n=0; n<Num_fact; n++)
  {
	m=0;
	while(m<p_fact)
	{
      rand_num = (int)((double)p*drand48());
      patt_picked = 0;		
      /*check that it hasn't been assigned before*/
      for(i=0; i<m; i++)
      {
		if(PC[n][i] == rand_num) patt_picked = 1;
	  }
	  if (patt_picked == 0) 
	  {
		PC[n][m] = rand_num;
		m++;
	  }
    }
  }
  
  //for(n=0; n<Num_fact; n++)
  //{
	//for(mu=0;mu<p_fact; mu++)
	//{
		//printf("%d \t", PC[n][mu]);
	//}
	//printf("\n");
  //}  
}

void SumFields(double a, int mu)
{

int i, n, m, k, rand_state;
double y, eigen_fact, expon;

/*set fields to 0*/
for(i=0;i<N;i++)
{
	for(k=0;k<S+1;k++)
	{
		h[i][k] = 0.0;
	}
}

/*each pattern sums field coming from all its parents*/
for(n=0; n<Num_fact; n++)
{
    expon = -fact_eigen_slope*n;
	for(m=0;m<p_fact; m++)
	{
		if(PC[n][m] == mu)
		{
			for(i=0;i<N;i++)
			{
			  // component coming from parents	
		      y = (double)drand48();
		      if(y<=a_pf)
		      {	
		        eigen_fact = exp(expon)*y/a_pf;
				h[i][Factors[n][i]] += eigen_fact;
				//printf("%.2f \t", h[i][Factors[n][i]]);
		      }
		    }
		}
	}
}
/*for children that have no parents, or when a_p is too small (this is epsilon in the paper + thesis)*/
for(i=0;i<N;i++)
{
  //small input to a random state
  rand_state = (int)((double)(S)*drand48());
  h[i][rand_state] += eps*drand48();
}

/*find state of maximal field*/
double hM;
int SM;
for(i=0;i<N;i++)
{
	hM = 0.0;
	SM = S;
	for(k=0;k<S+1;k++)
	{
		//printf("%d \t %d \t %.2f \n", i, k, h[i][k]);
		if(h[i][k] > hM) 
		{
			hM = h[i][k];
			SM = k;
		}
	}
	hmax[i] = hM;
	smax[i] = SM;
	//printf("%.2f\t", hmax[i]);
	//printf("%d\n", smax[i]);
}

}

void SetPatterns(double a)
{
  int i, mu, count;
  for(mu=0;mu<p;mu++)
  {
	/*for each pattern sum the fields*/  
    SumFields(a, mu);
    
    mystruct X[N];
 
    for(i=0;i<N;i++)
    {
      /*first initialize each pattern to null state*/
      Patt[mu][i] = S;		
      
      X[i].index = i;
      X[i].field = hmax[i];
      X[i].state = smax[i];
    }
    /*sort the patterns in order of decreasing field keeping the indices*/
    qsort(X, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
   
    i=0;
    count=0;
    while (count<N*a && i<N)
    {
	  if (X[i].state != S) 
	  {
		Patt[mu][X[i].index] = X[i].state;
		count++;
      }
	i++;
    }
  }
}

int SavePatterns(double a)
{
int mu, i;
	
char *buffer; 
buffer = new char [1000*sizeof(char)];
sprintf(buffer, "pattern_S%d_a%.1f_apf%.2f_pfact%d_p%d_Numfact%d_zeta%.3f", S, a, a_pf, p_fact, p, Num_fact, fact_eigen_slope);
fpattern = fopen(buffer, "w");
int sp;
int num_children_generated = 0;

for(mu=0;mu<p;mu++) 
{
  /*check sparsity of patterns before saving them: those that have 0 parents and don't receive any field and are globally null are discarded*/
  sp=0;
  for(i=0;i<N;i++)
  {
	if(Patt[mu][i] !=S) sp++;
  }  
  //printf("mu%d %d--------------------\n", mu, sp);
  
  if (sp == N*a)
  {
	num_children_generated++;
    for(i=0;i<N;i++)
    {
      fprintf(fpattern, "%d ", Patt[mu][i]);
    }
  fprintf(fpattern, "\n");
  }
}
fclose(fpattern);
return num_children_generated;	
}

/*----------------------------------------------------------------------*
*																		*
*							Main										*
* 																		*		
* ----------------------------------------------------------------------*/

int main(int argc, char **argv)
{

MPI_Init(&argc, &argv);


// Get the number of processes
int np;
MPI_Comm_size(MPI_COMM_WORLD, &np);
printf("np=%d\n", np);

// Get the rank of the process
int proc_id;
MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

printf("proc_id=%d\n", proc_id);

fact_eigen_slope = atof(argv[1]);

srand48(2017);

GetMemory(a);

/*Generate the factors*/
SetFactors(a);

/*Assign children to parents*/
AssignChildren(a);
//printf("here\n");

/*Generate the patterns*/
SetPatterns(a);

int num_children_generated;
/* to save patterns with different sparsities in different files */
num_children_generated = SavePatterns(a);
printf("%d\n", num_children_generated);

MPI_Finalize();

return 0;

}
