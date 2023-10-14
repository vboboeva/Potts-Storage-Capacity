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
	
int i, j, k, l, mu, rho;

patterns = new int*[p];
for(mu=0; mu<p; mu++)
  patterns[mu] = new int[N];
  
C0 = new double*[p];
for(mu=0; mu<p; mu++)
  C0[mu] = new double[p];

C1 = new double*[p];
for(mu=0; mu<p; mu++)
  C1[mu] = new double[p];

C1C1 = new double*[p];
for(mu=0; mu<p; mu++)
  C1C1[mu] = new double[p];
  
C2 = new double*[p];
for(mu=0; mu<p; mu++)
  C2[mu] = new double[p];

C2C2 = new double*[p];
for(mu=0; mu<p; mu++)
  C2C2[mu] = new double[p];

C1C2 = new double*[p];
for(mu=0; mu<p; mu++)
  C1C2[mu] = new double[p];

C3 = new double*[p];
for(mu=0; mu<p; mu++)
  C3[mu] = new double[p];
  
C4 = new double*[p];
for(mu=0; mu<p; mu++)
  C4[mu] = new double[p];

C = new double*[p];
for(mu=0; mu<p; mu++)
  C[mu] = new double[p];

char *buffer1;
buffer1 = new char [200*sizeof(char)];

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn0_file_a%.1f", S, a_pf, fact_eigen_slope, a);
//C0_file = fopen(buffer1, "w");

sprintf(buffer1, "Cmn1_file_a%.2f_apf%.2f_pfact%d_p%d_Nfact%d_Numfact%d_zeta%.3f", a, a_pf, p_fact, p, N_fact, Num_fact, fact_eigen_slope);
C1_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn1Cmn1_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C1C1_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn2_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C2_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn2Cmn2_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C2C2_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn1Cmn2_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C1C2_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn3_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C3_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn4_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
//C4_file = fopen(buffer1, "w");

//sprintf(buffer1, "S%d_apf%.2f_zeta%.1f/Cmn_file_a%.1f", S,  a_pf, fact_eigen_slope, a);
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

for(mu=0;mu<p;mu++)
{
  for(rho=0; rho<p; rho++)
  {
      C0[mu][rho]=0;
      C1[mu][rho]=0;
      C1C1[mu][rho]=0;      
      C2[mu][rho]=0;
      C2C2[mu][rho]=0;
      C1C2[mu][rho]=0;
      C3[mu][rho]=0;
      C4[mu][rho]=0;
     
   for(i=0;i<N;i++)
   {
	/*co-inactive*/
	C0[mu][rho]+=(patterns[mu][i]==S)*(patterns[rho][i]==S);

	for(k=0;k<S;k++)
	{
	  /*co-active and in SAME states*/
	  C1[mu][rho]+=(patterns[mu][i]==k)*(patterns[rho][i]==k);

	  /*inactive and active*/
	  C3[mu][rho]+=(patterns[mu][i]==S)*(patterns[rho][i]==k);
	  
	  /*active and inactive*/
	  C4[mu][rho]+=(patterns[mu][i]==k)*(patterns[rho][i]==S);
	  
	  /*co-active but DIFFERENT states*/
	  
	  for(l=0;l<S;l++)
	  {
	    if (k != l)
	    {
	      C2[mu][rho]+=(patterns[mu][i]==k)*(patterns[rho][i]==l);
	    }
	  }
	}
   }

      C0[mu][rho]=C0[mu][rho]/double(N);//(N*(1.-a));
      C1[mu][rho]=C1[mu][rho]/double(N);//(N*a);
      C2[mu][rho]=C2[mu][rho]/double(N);//(N*a);
      C3[mu][rho]=C3[mu][rho]/double(N);//(N*a);
      C4[mu][rho]=C4[mu][rho]/double(N);//(N*a);
      
      C1C2[mu][rho]=C1[mu][rho]*C2[mu][rho];  
      C1C1[mu][rho]=C1[mu][rho]*C1[mu][rho];
      C2C2[mu][rho]=C2[mu][rho]*C2[mu][rho];

      
      //fprintf(C0_file, "%f	", C0[mu][rho]);
      fprintf(C1_file, "%f	", C1[mu][rho]);
      //fprintf(C1C1_file, "%f	", C1C1[mu][rho]);
      //fprintf(C2_file, "%f	", C2[mu][rho]);
      //fprintf(C2C2_file, "%f	", C2C2[mu][rho]);
      //fprintf(C1C2_file, "%f	", C1C2[mu][rho]);
      //fprintf(C3_file, "%f	", C3[mu][rho]);
      //fprintf(C4_file, "%f	", C4[mu][rho]);
      //fprintf(C_file, "%f	", C0[mu][rho]+C1[mu][rho]+C2[mu][rho]+C3[mu][rho]+C4[mu][rho]);
     
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
