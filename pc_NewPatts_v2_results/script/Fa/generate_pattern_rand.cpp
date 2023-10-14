/*----------------------------------------------------------------------*
*		Generate Random Patterns written by Vizhe						*		
* ----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
*																		*
*		No need to enter parameters, done in main!						*
* 																		*		
* ----------------------------------------------------------------------*/

/*network parameters*/
extern int N;
extern int S;

FILE *randpattern;

extern int **xi;

typedef struct {
int index;
double field;
int state;
} randstruct;
  

int cmpfrand (const void *x, const void *y)
{
  double xx = (*(randstruct*)x).field;
  double yy = (*(randstruct*)y).field;
  /*sort by decreasing order*/
  if (xx > yy) return -1;
  if (xx < yy) return  1;
  return 0;
}

void GenerateRandomPatterns(double a, int p)
{
  int i, mu, nu;
  
  xi= new int*[p];
  for(i=0; i<p; i++)
    xi[i]=new int[N];
    
  char *buffer; 
  buffer = new char [1000*sizeof(char)];
  sprintf(buffer, "pattern_S%d_a%.1f", S, a);
  randpattern = fopen(buffer, "w");    
  
  for(mu=0;mu<p;mu++)
  {
    randstruct X[N];
 
    for(i=0;i<N;i++)
    {
      /*first initialize each pattern to null state*/
      xi[mu][i] = S;		
      
      X[i].index = i;
      X[i].field = drand48();
      X[i].state = (int)((double)S*drand48());
    }
    /*sort the patterns in order of decreasing field keeping the indices*/
    qsort(X, N, sizeof(randstruct), cmpfrand);  

    /*set N*a largest of them to become active*/
    for(i=0;i<N*a;i++)
    {
	  xi[mu][X[i].index] = X[i].state; 
    }  
  }
 
}


void SaveRandomPatterns(double a, int p)
{
int mu, i;
	
char *buffer; 
buffer = new char [1000*sizeof(char)];
sprintf(buffer, "randpattern_S%d_a%.1f_apf%.2f", S, a, a_pf);
randpattern = fopen(buffer, "w");

for(mu=0;mu<p;mu++) 
{
  for(i=0;i<N;i++)
    {
      fprintf(randpattern, "%d ", xi[mu][i]);
    }
  fprintf(randpattern, "\n");
}

fclose(randpattern);
}
