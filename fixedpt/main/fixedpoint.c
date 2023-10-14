#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

#define cdivN 1.0
#define Nxi 100
#define Nz 1000
#define S 5
#define a 0.5
#define alpha 2.4
#define U 0.5
#define eps 0.00001
#define T 200
#define	myexp 0.5
#define Runs 1

/* order parameters are global */
double m;
double q;
double o;

/* additional variables */
double psi;
double rho;

/* vectors containing the Nxi samples of the xi variable, for m, q and omega */
int xim[Nxi], xiq[Nxi], xio[Nxi];

/* vectors containing the Nz samples of S gaussian variables, for m, q and omega */
double zm[(S+1)*Nz], zq[(S+1)*Nz], zo[(S+1)*Nz];

/* mean-field vectors, for m, q and omega */
double Hm[S+1], Hq[S+1], Ho[S+1];

/*
 *  CAN BE PUT IN A MATRIX (TABLE)
 */
double v (int x, int y){
    if(y!=0)
    {
        if (x!=y) return -a/(double)S;
        else return 1-a/(double)S;
    }
    else return 0.;
}

double max (double x[], int dim){
    int i;
    double r = x[0];
    for (i=1; i<dim; i++){
        if (x[i]>r) r = x[i];
    }
    return r;
}

int max_index (double x[], int dim){
    int i;
    int max=0;
    double aux=x[0];
    for (i=1; i<dim; i++){
        if (x[i]>=aux){
			aux=x[i];
			max=i;
		}
    }
    return max;
}

double mean_field(int xi, int k, double z[], int index){
    int n;
    double s;
    if (k!=0){
        s=0.;
        for (n=1; n<=S; n++)  s += v(n,k)*z[(S+1)*index+n];
    
        return v(xi,k)*m + alpha*cdivN*psi/(2*S) + s*rho - U;
    }
    else return 0.;
}

void find_fixed_pt()
{
	FILE *file;
	char buf[0x100];
	snprintf(buf, sizeof(buf), "alpha%.2f.txt", alpha);
	file = fopen(buf, "w");	
	
	int run, check, ts;
    int k, l, n;
    int lxi, lz;
    double m_old, q_old, o_old;
    double Fm, Fq, Fo;
    double Ym, Yq, Yo;
    double sum;
    double eta;
    
    /* values and probabilities for the random variable xi */
    int xi[S+1];
    double p[S+1];
    xi[0] = 0;
    p[0] = 1-a;
    for (k=1; k<=S; k++){
        xi[k] = k;
        p[k] = a/(double)S;
    }
		
		for(run=0; run<Runs; run++)
		{
			/* Initialization of the random number generator */
			int seed = 1987;/*rlxd_seed();*/
			rlxd_init(1,seed);

			/* Initial values of the order parameters */
			ts=0;
			m = 10.0;
			q = 0.7;
			o = 0.1;
			
			psi = (o/(double)S)/(1.-o/(double)S);
			rho = sqrt((alpha*a/(double)S)*q*(1. + 2.*cdivN*psi + cdivN*psi*psi)/((double)S*(1.-a/(double)S)));
			
			printf("%f\t%f\t%f\t%f\t%d\n",a,U,myexp,alpha,S);
		
			while(ts<T)
			{
				eta = 1./pow((double)ts+1.,myexp);
				
				/* Extract discrete random numbers with prob mass function p */
				Fm=0.; discr_int(xim, xi, p, S+1, Nxi);
				Fq=0.; discr_int(xiq, xi, p, S+1, Nxi);
				Fo=0.; discr_int(xio, xi, p, S+1, Nxi);
				
				for(lxi=0; lxi<Nxi; lxi++)
				{
					/*printf("%d\t%d\t%f\t%f\t%f\n", ts, lxi, Fm, Fq, Fo);*/

					/* Extract standard gaussian random numbers */
					Ym=0.; gauss_dble(zm, (S+1)*Nz);
					Yq=0.; gauss_dble(zq, (S+1)*Nz); 
					Yo=0.; gauss_dble(zo, (S+1)*Nz);
					
					for (lz=0; lz<Nz; lz++)
					{
						/* function for m */
						for (n=0; n<=S; n++) Hm[n] = mean_field(xim[lxi],n,zm,lz); /*printf("%f\n",mean_field(xim[lxi],n,zm,lz));*/
						l=max_index(Hm,S+1);
						Ym += v(xim[lxi],l)/((double)Nz);
			
						/* function for q */
						for (n=0; n<=S; n++) Hq[n] = mean_field(xiq[lxi],n,zq,lz);
						l=max_index(Hq,S+1);
						if (l!=0) Yq += 1./((double)Nz);
			
						/* function for o */
						for (n=0; n<=S; n++) Ho[n] = mean_field(xio[lxi],n,zo,lz);
						l=max_index(Ho,S+1);
						sum=0.;
						for (k=1; k<=S; k++) sum += v(k,l)*zo[(S+1)*lz+k]; /* k also from 0???? */
						Yo += sum/((double)Nz);
					}

					/* function for m */
					Fm += Ym/((double)Nxi);
			
					/* function for q */
					Fq += Yq/((double)Nxi);
			
					/* function for o */
					Fo += Yo/((double)Nxi);
				}

				/* backup current order parameters */
				m_old = m;
				q_old = q;
				o_old = o;

				/* update order parameters */
				m = (1.-eta)*m_old + eta*Fm/(a*(1.-a/(double)S));
				q = (1.-eta)*q_old + eta*Fq/a;
				o = (1.-eta)*o_old + eta*Fo*(double)S/sqrt(alpha*q_old*(1. + 2.*cdivN*psi + cdivN*psi*psi)*a*(1.-a/(double)S));

				psi = (o/(double)S)/(1.-o/(double)S);
				rho = sqrt((alpha*a/(double)S)*q*(1. + 2.*cdivN*psi + cdivN*psi*psi)/((double)S*(1.-a/(double)S)));

				/* check convergence 
				if(sqrt((pow(m-m_old,2)+pow(q-q_old,2)+pow(o-o_old,2)))<eps) check=1; */
				printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", ts, eta, psi, rho, m, q, o);

				ts++;
			}
			fprintf(file, "%d\t%lf\t%lf\t%lf\n ", ts, m, q, o);        
			fflush(file);  
		}
		fclose(file);
}

int main ()
{
find_fixed_pt();
return 0;
}
