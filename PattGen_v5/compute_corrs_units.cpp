#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "const_potts.h"

FILE *C0_file;
FILE *C1_file;
FILE *C1C1_file;
FILE *C2_file;
FILE *C2C2_file;
FILE *C1C2_file;
FILE *C3_file;
FILE *C4_file;
FILE *C_file;
FILE *pat;

int **patterns;
double **C0;
double **C1;
double **C1C1;
double **C2;
double **C2C2;
double **C1C2;
double **C3;
double **C4;
double **C;

/*----------------------------------------------------------------------*
*																		*
*							Main										*
* 																		*		
* ----------------------------------------------------------------------*/

int main()
{
	
int i, j, k, l, mu;

patterns = new int*[p];
for(mu=0; mu<p; mu++)
  patterns[mu] = new int[N];
  
C0 = new double*[N];
for(i=0; i<N; i++)
  C0[i] = new double[N];

C1 = new double*[N];
for(i=0; i<N; i++)
  C1[i] = new double[N];

C1C1 = new double*[N];
for(i=0; i<N; i++)
  C1C1[i] = new double[N];
  
C2 = new double*[N];
for(i=0; i<N; i++)
  C2[i] = new double[N];

C2C2 = new double*[N];
for(i=0; i<N; i++)
  C2C2[i] = new double[N];

C1C2 = new double*[N];
for(i=0; i<N; i++)
  C1C2[i] = new double[N];

C3 = new double*[N];
for(i=0; i<N; i++)
  C3[i] = new double[N];
  
C4 = new double*[N];
for(i=0; i<N; i++)
  C4[i] = new double[N];

C = new double*[N];
for(i=0; i<N; i++)
  C[i] = new double[N];

char *buffer1;
buffer1 = new char [200*sizeof(char)];

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij0_file_a%.1f", S, a_pf, fact_eigen_slope, a);
//C0_file = fopen(buffer1, "w");

sprintf(buffer1, "Cij1_file_a%.2f_apf%.2f_pfact%d_p%d_Nfact%d_Numfact%d_zeta%.3f", a, a_pf, p_fact, p, N_fact, Num_fact, fact_eigen_slope);
C1_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij1Cij1_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C1C1_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij2_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C2_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij2Cij2_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C2C2_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij1Cij2_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C1C2_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij3_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C3_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij4_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C4_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cij_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C_file = fopen(buffer1,"w");

char *buffer;

buffer = new char [200*sizeof(char)];

sprintf(buffer, "pattern_S%i_a%.2f_apf%.2f_pfact%d_Nfact%d_Numfact%d_zeta%.3f", S, a, a_pf, p_fact, N_fact, Num_fact, fact_eigen_slope);

pat = fopen(buffer,"r");

/*read all patterns*/
for(mu=0;mu<p;mu++)
{
  for(i=0;i<N;i++)
  {
    fscanf(pat, "%d", &patterns[mu][i]);
    
  }
}

fclose(pat);

for(i=0;i<N;i++)
{
  for(j=0; j<N; j++)
  {
      C0[i][j]=0;
      C1[i][j]=0;
      C1C1[i][j]=0;      
      C2[i][j]=0;
      C2C2[i][j]=0;
      C1C2[i][j]=0;
      C3[i][j]=0;
      C4[i][j]=0;
     
   for(mu=0;mu<p;mu++)
   {
	/*co-inactive*/
	C0[i][j]+=(patterns[mu][i]==S)*(patterns[mu][j]==S);

	for(k=0;k<S;k++)
	{
	  /*co-active and in SAME states*/
	  C1[i][j]+=(patterns[mu][i]==k)*(patterns[mu][j]==k);

	  /*inactive and active*/
	  C3[i][j]+=(patterns[mu][i]==S)*(patterns[mu][j]==k);
	  
	  /*active and inactive*/
	  C4[i][j]+=(patterns[mu][i]==k)*(patterns[mu][j]==S);
	  
	  /*co-active but DIFFERENT states*/
	  
	  for(l=0;l<S;l++)
	  {
	    if (k != l)
	    {
	      C2[i][j]+=(patterns[mu][i]==k)*(patterns[mu][j]==l);
	    }
	  }
	}
   }

      C0[i][j]=C0[i][j]/double(p);//(N*(1.-a));
      C1[i][j]=C1[i][j]/double(p);//(N*a);
      C2[i][j]=C2[i][j]/double(p);//(N*a);
      C3[i][j]=C3[i][j]/double(p);//(N*a);
      C4[i][j]=C4[i][j]/double(p);//(N*a);
      
      C1C2[i][j]=C1[i][j]*C2[i][j];  
      C1C1[i][j]=C1[i][j]*C1[i][j];
      C2C2[i][j]=C2[i][j]*C2[i][j];

      
      //fprintf(C0_file, "%f	", C0[i][j]);
      fprintf(C1_file, "%f	", C1[i][j]);
      //fprintf(C1C1_file, "%f	", C1C1[i][j]);
      //fprintf(C2_file, "%f	", C2[i][j]);
      //fprintf(C2C2_file, "%f	", C2C2[i][j]);
      //fprintf(C1C2_file, "%f	", C1C2[i][j]);
      //fprintf(C3_file, "%f	", C3[i][j]);
      //fprintf(C4_file, "%f	", C4[i][j]);
      //fprintf(C_file, "%f	", C0[i][j]+C1[i][j]+C2[i][j]+C3[i][j]+C4[i][j]);
     
      //fflush(C0_file);
      fflush(C1_file);
      //fflush(C1C1_file);
      //fflush(C2_file);
      //fflush(C2C2_file);
      //fflush(C1C2_file);
      //fflush(C3_file);
      //fflush(C4_file);
      //fflush(C_file);
  }
  //fprintf(C0_file, "\n");
  fprintf(C1_file, "\n");
  //fprintf(C1C1_file, "\n");
  //fprintf(C2_file, "\n");
  //fprintf(C2C2_file, "\n");
  //fprintf(C1C2_file, "\n");
  //fprintf(C3_file, "\n");
  //fprintf(C4_file, "\n");
  //fprintf(C_file, "\n");
}

//fclose(C0_file);
fclose(C1_file);
//fclose(C1C1_file);
//fclose(C2_file);
//fclose(C2C2_file);
//fclose(C1C2_file);
//fclose(C3_file);
//fclose(C4_file);
//fclose(C_file);

}
