#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*----------------------------------------------------------------------*
*		New pattern generator written by Vizhe, version 5,				*
* 		that defines factors as random patterns acting on random		*
* 		subsets of patterns	and not units (as does Ale, version 0).		*
* 		The objective is to reverse the									*
* 		correlation statistics of Ale's algorithm, that is,				*
* 		skewed correlations distribution between patterns, 			 	*
* 		and symmetric between units. Difference with version 1: 		*
* 		parents are random patterns that send influence in a 			*
* 		different Potts state for each unit of a pattern.				*
* 		Difference with version 3: random input (eps) included to solve		*
* 		small sparsity for small values of a_pf.						*
* 		Difference with version 4: parents are defined also on 			*
* 		subsets of units as well as subsets of patterns					*
* 																		*
*Limit cases:															*
*1) a_pf ---> 0 gives us randomly correlated patterns					*
*2) f*Num_fact = 1 yields (effectively) ultrametric patterns(f=p_fact/p)*
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
int p_fact = 30; 
int N_fact = 1000; /*Setting N_fact=N, we recover version 4 of algorithm*/
int Num_fact = 150;

double eps = 0.000001;

/*correlation parameters*/
double a = 0.25;
double fact_eigen_slope = 0.005;
double a_pf = 0.4;

int **Factors;
int **Patt;
int **PC;
int **PU;
double **hh;
double *hmax;
int *smax;

FILE *fpattern;
FILE *fparents;


/*----------------------------------------------------------------------*
*																		*
*							FUNCTIONS									*
* 																		*		
* ----------------------------------------------------------------------*/

void GetMemory()
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
  
PU = new int*[Num_fact];
for(i=0; i<Num_fact; i++)
  PU[i] = new int[N_fact];
 
hh = new double*[N];
for(i=0; i<N; i++)
  hh[i] = new double[S+1];

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
void SetFactors()
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

void AssignChildren()
{
  int n, m, i, patt_picked, unit_picked;
  int rand_num;

  for(n=0; n<Num_fact; n++)
  {
	/*assign children to parents*/
	m=0;
	while(m<p_fact)
	{
      rand_num = (int)((double)p*drand48());
      patt_picked = 0;		
      /*check that it hasn't been assigned before*/
      for(i=0; i<m; i++)
      {
		if(PC[n][i] == rand_num) {patt_picked = 1;}/*{printf("%d \t %d\t %f\t %d \n", n, PC[n][i], rand_num, patt_picked);}*/
	  }
	  if (patt_picked == 0) 
	  {
		PC[n][m] = rand_num;
		m++;
	  }
    }
    /*assign units to parents*/
	m=0;
	while(m<N_fact)
	{
      rand_num = (int)((double)N*drand48());
	  /*printf("%i\t%d\n", m, rand_num);*/ 
      unit_picked = 0;		
      /*check that it hasn't been assigned before*/
      for(i=0; i<m; i++)
      {
		if(PU[n][i] == rand_num) unit_picked = 1;
	  }
	  if (unit_picked == 0) 
	  {
		/*printf("%d\t%d\n", m, rand_num);*/
		PU[n][m] = rand_num;
		m++;
	  }
    }    
  }
}

void SumFields(int mu)
{

int i, n, m, k, rand_state;
double y, eigen_fact, expon;

/*set fields to 0*/
for(i=0;i<N;i++)
{
	for(k=0;k<S+1;k++)
	{
		hh[i][k] = 0.0;
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
			/*printf("%d \t %d \n", PC[n][m], mu);*/
			for(i=0;i<N_fact;i++)
			{
			  // component coming from parents	
		      y = (double)drand48();
		      if(y<=a_pf)
		      {	
		        eigen_fact = exp(expon)*y/a_pf;
				hh[PU[n][i]][Factors[n][i]] += eigen_fact;
				//printf("%.2f \t", hh[i][Factors[n][i]]);
		      }
		    }
		}
	}
}
/*for children that have no parents, or when a_p is too small (this is epsilon in the paper + thesis)*/
for(i=0;i<N;i++)
{
  /*small input to a random state*/
  rand_state = (int)((double)(S)*drand48());
  hh[i][rand_state] += eps*drand48();
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
		//printf("%d \t %d \t %.2f \n", i, k, hh[i][k]);
		if(hh[i][k] > hM) 
		{
			hM = hh[i][k];
			SM = k;
		}
	}
	hmax[i] = hM;
	smax[i] = SM;
	//printf("%.2f\t", hmax[i]);
	//printf("%d\n", smax[i]);
}

}

void SetPatterns()
{
  int i, mu, count;
  for(mu=0;mu<p;mu++)
  {
	/*for each pattern sum the fields*/  
    SumFields(mu);
    
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

void SavePatterns()
{
int mu, i;
	
char *buffer; 
buffer = new char [1000*sizeof(char)];
sprintf(buffer, "pattern_S%d_a%.2f_apf%.2f_pfact%d_Nfact%d_Numfact%d_zeta%.3f", S, a, a_pf, p_fact, N_fact, Num_fact, fact_eigen_slope);
fpattern = fopen(buffer, "w");

for(mu=0;mu<p;mu++) 
{
    for(i=0;i<N;i++)
    {
      fprintf(fpattern, "%d ", Patt[mu][i]);
    }
  fprintf(fpattern, "\n");
}
fclose(fpattern);
}

/*----------------------------------------------------------------------*
*																		*
*							Main										*
* 																		*		
* ----------------------------------------------------------------------*/

int main() /*int argc, char **argv*/
{
/*fact_eigen_slope = atof(argv[1]);*/

/*set seed to make the same dataset*/
srand48(2017);

GetMemory();

/*Generate the factors*/
SetFactors();

/*Assign children to parents*/
AssignChildren();

/*Generate the children*/
SetPatterns();

/*Save the children!*/
SavePatterns();

return 0;

}
