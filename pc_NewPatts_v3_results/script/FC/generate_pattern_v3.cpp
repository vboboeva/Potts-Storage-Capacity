/*----------------------------------------------------------------------*
*		New pattern generator written by Vizhe (v3),		*	
* ----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
*																		*
*		No need to enter parameters, done in main!						*
* 																		*		
* ----------------------------------------------------------------------*/

/*network parameters*/
extern int N;
extern int S;
extern int p;
extern double a;

/*pattern correlations*/
extern int p_fact; 
extern int Num_fact;
extern double eps;
extern double f;

/*correlation parameters*/
extern double a_pf;
extern double fact_eigen_slope;

extern int **xi;

int **Factors;
int **PC;
double **hh;
double *hmax;
int *smax;

FILE *corrpattern;

/*----------------------------------------------------------------------*
*																		*
*							FUNCTIONS									*
* 																		*		
* ----------------------------------------------------------------------*/

void GetMemory()
{
int i;
	
Factors = new int*[Num_fact];
for(i=0; i<Num_fact; i++)
  Factors[i] = new int[N];
  
PC = new int*[Num_fact];
for(i=0; i<Num_fact; i++)
  PC[i] = new int[p_fact];
 
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

void AssignChildren()
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

void SumFields(int mu)
{

int i, n, m, k;
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
			for(i=0;i<N;i++)
			{
		      y = (double)drand48();
		      if(y<=a_pf && Factors[n][i] !=S)
		      {	
		        eigen_fact = exp(expon)*y/a_pf;
				hh[i][Factors[n][i]] += eigen_fact + eps*drand48();
				//printf("%.2f \t", hh[i][Factors[n][i]]);
		      }
		    }
		}
	}
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
  int sp, i, mu, count;
  for(mu=0;mu<p;mu++)
  {
	/*for each pattern sum the fields*/  
    SumFields(mu);
    
    mystruct X[N];
 
    for(i=0;i<N;i++)
    {
      /*first initialize each pattern to null state*/
      xi[mu][i] = S;		
      
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
		xi[mu][X[i].index] = X[i].state;
		count++;
      }
	i++;
    }
  
  /*check sparsity of patterns before saving them: those that have 0 parents and don't receive any field and are globally null are discarded*/
  //sp=0;
  //if (mu == 280)
  //{
  //for(i=0;i<N;i++)
  //{
	////if(xi[mu][i] !=S) sp++;
	//printf("%d  ", xi[mu][i]);
  //} 
  //}
}
}

int SavePatterns()
{
int mu, i;
	
char *buffer; 
buffer = new char [1000*sizeof(char)];
sprintf(buffer, "corrpattern_S%d_a%.1f_apf%.2f_pfact%d_p%d_Numfact%d", S, a, a_pf, p_fact, p, Num_fact);
corrpattern = fopen(buffer, "w");
int sp;
int num_children_generated = 0;

for(mu=0;mu<p;mu++) 
{
  /*check sparsity of patterns before saving them: those that have 0 parents and don't receive any field and are globally null are discarded*/
  sp=0;
  for(i=0;i<N;i++)
  {
	if(xi[mu][i] !=S) sp++;
  }  
  //printf("mu%d %d--------------------\n", mu, sp);
  
  if (sp != 0)
  {
	num_children_generated++;
    for(i=0;i<N;i++)
    {
      fprintf(corrpattern, "%d ", xi[mu][i]);
    }
  fprintf(corrpattern, "\n");
  }
}
fclose(corrpattern);
return num_children_generated;	
}
