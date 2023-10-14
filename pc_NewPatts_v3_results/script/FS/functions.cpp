/*NET PARAMETERS*/
extern int N;
extern int Cm;
extern int S;
extern double U;
extern double nu;
extern int p;
extern double a;

/*pattern parameters*/
extern int tot_num_trials;
extern int PattSet;
extern int NumSet;

/*latching parameters*/
extern double b1;
extern double b2;
extern double b3;
extern double w;
extern double g;

/*network update*/
extern int Tnet; 
extern int time_toprint;

/*type of connectivity*/
extern int full_connectivity;
//if (full_connectivity == 1) Cm = N-1;
extern int random_dilution;
extern int symmetric_random_dilution;
extern int do_avg_Cij;
extern int variable_connectivity;

/*cue parameters and type*/
extern double fraction_units_cued;
extern int exp_cue;
extern int tau;
extern int full_cue;
extern int partial_cue;

/*graded or discrete update*/
extern int update_h;
extern int update_r;
extern int discrete;
extern int graded;
extern int beta;

/*pattern correlations*/
extern int p_fact;
extern int Num_fact;
extern double a_pf;
extern double fact_eigen_slope;
extern double eps;
extern double f;

extern int 			*number_connections;
extern double		*units_receive_cue; 
extern double    	*m;
extern double    	*m_sorted;
extern double		**sparsity_by_state;
extern double		m_retr;
extern double		m_cue;
extern int    	 	**xi;
extern double   	**s;
extern double    	****J; 
extern int       	retr;
extern int       	cue;
extern double   	**theta;
extern double    	**r;
extern double    	**h; 
extern int 			**Permut;
extern int 			**C;
extern double		***fin_state;
extern double		**C0;
extern double		**C1;
extern double		**C2;
extern double		H_interaction, H_quad;

extern FILE *storage_cap;
extern FILE *overlap_intime;
extern FILE *fin_sparsity;
extern FILE *fin_overlaps;
extern FILE *cij;
extern FILE *q_file;
extern FILE *state_file;
extern FILE *C0_file;
extern FILE *C1_file;
extern FILE *C2_file;
extern FILE *fin_state_file;
extern FILE *H_interaction_file;
extern FILE *H_quad_file;
extern FILE *field_file;

/*compare function for sorting algorithm*/
int compare(const void *x, const void *y)
{
  int ix = *(int *)x;
  int iy = *(int *)y;
  return m_sorted[ix] > m_sorted[iy] ? -1 : m_sorted[ix] < m_sorted[iy]; // decreasing order
}
/*----------------------------------------------------------*/
/*prints overlap of network state with each pattern in time and takes time as input */
void compute_sparsity_by_state()
{
int mu, k, i;	
for(mu=0;mu<p;mu++)
{
	for(k=0;k<S;k++)
	{
		sparsity_by_state[mu][k]=0.;
		for(i=0; i<N; i++)
		{
			if (xi[mu][i] == k) 
			{
				sparsity_by_state[mu][k]+=1./double(N);
			}
		}
	//printf("sparsity_by_state=%.2f\n", sparsity_by_state[mu][k]);
	}
}
}
/*----------------------------------------------------------*/
/*prints overlap of network state with each pattern in time and takes time as input */
double compute_sparsity_net()
{
int kmax, k, i;	
double smax, sp;
sp=0.0;
for(i=0; i<N; i++)
{
	kmax = -1;
	smax = -1.0;
	/*find state of maximal field*/
	for(k=0; k<S+1; k++)
	{
		if(s[i][k] > smax) 
		{
			smax = s[i][k];
			kmax = k;
		}
	}
	if(kmax != S) sp++;
}
sp=sp/(float)N;
return sp;
}

/*-----------------------------------------------------------*/
/*computes overlap of network state with each pattern*/

void compute_m()
{
int  k, i, mu;
double ma, maa, invdenN;

invdenN=1./(a*(1.-a/S)*(double)N); // normalizing factor of overlap

for(mu=0;mu<p;mu++)
{
  maa=0.;
  for(i=0;i<N;i++)
  {
    ma=0.;
    for(k=0;k<S;k++)
    {
      ma+=((double)(xi[mu][i]==k)-a/(float)S)*s[i][k]; //to calculate m  
    }
    maa+=ma;
  }
m[mu]=maa*invdenN; //value of m[mu] for each mu
}
}
/*-----------------------------------------------------------*/
/*sorts overlap*/

void compute_m_sorted()
{
int indices[p]; // contains indices of overlaps
int mu;
for(mu=0;mu<p;mu++)
{
  indices[mu] = mu;
}
memcpy(m_sorted, m, sizeof(double)*p); // copies overlaps inside global array m_sorted
qsort(indices, p, sizeof (*indices), compare); // sorts m_sorted, keeping track of indices

if (indices[0] == cue) 
{
	retr = cue;
	m_retr = m[indices[0]];
	m_cue = m[indices[0]];
}
else 
{
	retr = indices[0];
	m_retr = m[indices[0]];
}
}

/*-----------------------------------------------------------*/
/*computes overlap of network state with pattern CUED*/

void compute_m_cue()
{
int  k, i;
double ma, maa, invdenN;

invdenN=1./(a*(1.-a/S)*(double)N); // normalizing factor of overlap

maa=0.;
for(i=0;i<N;i++)
{
  ma=0.;
  for(k=0;k<S;k++)
  {
    ma+=((double)(xi[cue][i]==k)-a/(float)S)*s[i][k]; //to calculate m  
  }
  maa+=ma;
}
m_cue=maa*invdenN; //value of m[mu] for each mu
}

///*-----------------------------------------------------------*/
///*computes fraction of units co-inactive, co-same-state-active, co-active*/
//void compute_C0_C1_C2(double a)
//{
//int i,j,k,i_c,i_cc;

//for(i_c=0;i_c<tot_num_trials;i_c++)
//{
  //for(i_cc=0; i_cc<tot_num_trials; i_cc++)
  //{
    //C0[i_c][i_cc]=0;
    //C1[i_c][i_cc]=0;
    //C2[i_c][i_cc]=0;
    
    //for(i=0;i<N;i++)
    //{
////       printf("i=%d\n", i);
      
      //C0[i_c][i_cc]+=(fin_state[i_c][i][S]==fin_state[i_cc][i][S])*(fin_state[i_cc][i][S]==1);
      //for(k=0;k<S;k++)
      //{
//// 	printf("k=%d\n", k);
	//C1[i_c][i_cc]+=(fin_state[i_c][i][k]==fin_state[i_cc][i][k])*(fin_state[i_cc][i][k]==1);
//// 	printf("%f %f %f \n",fin_state[i_c][i][k], fin_state[i_cc][i][k], C1[i_c][i_cc]);
	
	//for(j=0;j<S;j++)
	//{
	  //if (k != j)
	  //{
	    //C2[i_c][i_cc]+=(fin_state[i_c][i][k]==fin_state[i_cc][i][j])*(fin_state[i_cc][i][j]==1);
	  //}
	//}
      //}
    //}
    //C0[i_c][i_cc]=C0[i_c][i_cc]/double(N);//(N*(1.-a));
    //C1[i_c][i_cc]=C1[i_c][i_cc]/double(N);//(N*a);
    //C2[i_c][i_cc]=C2[i_c][i_cc]/double(N);//(N*a);
    
    //if (i_cc > i_c)
    //{
      //fprintf(C0_file, "%f	", C0[i_c][i_cc]);
      //fprintf(C1_file, "%f	", C1[i_c][i_cc]);
      //fprintf(C2_file, "%f	", C2[i_c][i_cc]);
      //fflush(C0_file);
      //fflush(C1_file);
      //fflush(C2_file);
    //}
  //}

//}
//fprintf(C0_file, "\n");
//fprintf(C1_file, "\n");
//fprintf(C2_file, "\n");
//}
/*-----------------------------------------------------------*/
/*computes fraction of units co-inactive, co-same-state-active, co-active*/
double compute_info()
{
int i,k,l;
double I;

//double I0,I1,I2,I3;
//double Cm0,Cr0,Cmr00; // fraction of units co-inactive
//Cm0=0.;
//Cr0=0.;
//Cmr00=0.;
//for(i=0;i<N;i++)
//{
  //Cm0+=(xi[mu][i]==S);
  //Cr0+=(s[i][S]==1);
  //Cmr00+=(xi[mu][i]==S)*(s[i][S]==1);
//}
//Cm0=Cm0/N;
//Cr0=Cr0/N;
//Cmr00=Cmr00/N;
////printf("%.3f \t %.3f \t %.3f\n\n", Cmr00,Cm0,Cr0);
//if(Cmr00 !=0) I0 = Cmr00*log(Cmr00/(Cm0*Cr0));

/////*-------------------------*/

//double Crk,Cmr0k; // fraction of units inactive in memory but activated in retrieved
//for(k=0;k<S;k++)
//{
	//Crk=0.;
	//Cmr0k=0.;
	//for(i=0;i<N;i++)
	//{
		//Crk+=(s[i][k]==1);
		//Cmr0k+=(xi[mu][i]==S)*(s[i][k]==1);
	//}
	//Crk=Crk/N;
	//Cmr0k=Cmr0k/N;
	////printf("%.3f \t %.3f \t %.3f\n", Cmr0k, Cm0, Crk);
	//if(Cmr0k !=0) I1+=Cmr0k*log(Cmr0k/(Cm0*Crk));
//}
//printf("\n");
/////*-------------------------*/

//double Cmk,Cmrk0; // fraction of units active in memory but inactive in retrieved
//for(k=0;k<S;k++)
//{
	//Cmk=0.;
	//Cmrk0=0.;
	//for(i=0;i<N;i++)
	//{
		//Cmk=(xi[mu][i]==k);
		//Cmrk0+=(xi[mu][i]==k)*(s[i][S]==1);
	//}
	//Cmk=Cmk/N;
	//Cmrk0=Cmrk0/N;
	////printf("%.3f \t %.3f \t %.3f\n",Cmrk0, Cmk, Cr0);
	//if(Cmrk0 != 0) I2+=Cmrk0*log(Cmrk0/(Cmk*Cr0));
//}
////printf("\n");
/////*-------------------------*/

//double Cmrkl, Crl;
//for(k=0;k<S;k++)
//{
	//for(l=0;l<S;l++)
	//{
		//Cmrkl=0.;
		//Cmk=0.;
		//Crl=0.;
		//for(i=0;i<N;i++)
		//{
			//Cmk+=(xi[mu][i]==k);
			//Crl+=(s[i][l]==1);
			//Cmrkl+=(xi[mu][i]==k)*(s[i][l]==1);
		//}
		//Cmk=Cmk/N;
		//Crl=Crl/N;
		//Cmrkl=Cmrkl/N;
		////printf("%.3f %.3f %.3f\n",Cmk, Crl, Cmrkl);

		//if(Cmrkl != 0) I3+=Cmrkl*log(Cmrkl/(Cmk*Crl));
	//}
//}

/*test*/
double C, Cr, Cm;
I=0.;
for(k=0;k<S+1;k++)
{
	for(l=0;l<S+1;l++)
	{
		C=0.;
		Cm=0.;
		Cr=0.;
		for(i=0;i<N;i++)
		{
			Cm+=(xi[cue][i]==k);
			Cr+=(s[i][l]==1);
			C+=(xi[cue][i]==k)*(s[i][l]==1);
		}
		C=C/N;
		Cm=Cm/N;
		Cr=Cr/N;
		if(C!=0) I+=C*log(C/(Cm*Cr));
	}
}
//printf("testI=%.3f\n", I);
return I/(log(2.));
}

/*-----------------------------------------------------------*/
/*constructs C and J*/
void construct_CJ()
{

int i, j, l, k, mu;
int i_c, x, new_one;
//for(i=0; i<N; i++)
//{
  //for(j=0; j<N; j++)
  //{
    //C[i][j] = N;
  //}
//}

/*Cij full*/
if (full_connectivity == 1)
{
  for (i=0; i<N; i++) number_connections[i]=Cm;

  printf("full_connectivity\n");

  for(i=0; i<N; i++)
  {
	//printf("i=%d %d\n", i, number_connections[i]);
	k=0;
    for(x=0; x<number_connections[i]; x++)
    {
      if (i != x) 
      {
		  C[i][k] = x;
		  //printf("%d  ", C[i][k]);
		  k++;
	  }
    }
  }
}

/*Cij randomly diluted*/
else if (random_dilution == 1)
{
  for (i=0; i<N; i++) number_connections[i]=Cm;
  
  printf("random_dilution\n");
  for(i=0; i<N; i++)
  {
    i_c = 0;  
    while(i_c<Cm)
    {
      new_one = 1;
      j = (int)((double)N*drand48());
      if(j==i) new_one = 0; 
      for(x=0; x<i_c; x++) 
      {
		if(C[i][x]==j) new_one = 0;
      }
      if(new_one)
      {
		C[i][i_c] = j;
		i_c++;
      }
    }
  }
}

/*Cij randomly diluted but variable*/
else if (variable_connectivity == 1)
{
  printf("variable_dilution\n");
	
  int i_c;
  double rand_num;
  
  for (i=0; i<N; i++) number_connections[i]=0;
  
  for(i=0; i<N; i++)
  {    
	i_c=0;  
    for(x=0;x<N;x++)
    {
      rand_num = drand48();
      if(rand_num <= (double)Cm/(double)N)
      {
        C[i][x] = (int)((double)N*drand48());
        //printf("%d\n", C[i][x]);
        i_c++;
      }
    }
    number_connections[i] = i_c;
    printf("num_conn=%d\n", number_connections[i]);
  }
}

///*Cij symmetrically randomly diluted*/
//else if (symmetric_random_dilution == 1)
//{
  //for (i=0; i<N; i++) number_connections[i]=Cm;

  //printf("symmetric_random_dilution\n");
  
  //igraph_integer_t nodes = N; /* number of nodes */
  //int edges = Cm; /* fixed degree*/

  //igraph_t graph;/* pointer to the main graph object */

  ///* initialization igraph vectors */
  //igraph_k_regular_game(&graph, nodes, edges, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

  ///*make adjacency matrix*/
  //igraph_matrix_t m;
  //igraph_matrix_init(&m, nodes, nodes);
  //igraph_matrix_null(&m);

  //igraph_get_adjacency(&graph, &m, IGRAPH_GET_ADJACENCY_BOTH, false);

  //for(i=0; i<nodes; i++)
  //{
    //for(j=0; j<nodes; j++)
    //{
      //C[i][j] = (int)igraph_matrix_e(&m,i,j);
    //}
  //}
  ///* destroy graph and free space */
  //igraph_destroy(&graph);
  
//}

/*-------------------------------Jixkl-----------------------------------*/
//printf("before Jixkl	\n");

for(i=0; i<N; i++)
{
  for(x=0; x<number_connections[i]; x++)
  {
    for(k=0; k<S; k++)
    {
      for(l=0; l<S; l++)
      {
	    J[i][x][k][l]=0;
        for(mu=0; mu<p; mu++)
        {
          /*standard covariance rule*/
	      J[i][x][k][l]+=((double)(xi[mu][i]==k)-a/(float)S)*((double)(xi[mu][C[i][x]]==l)-a/(float)S);
          /*popularity rule*/
          //J[i][x][k][l]+=((double)(xi[mu][i]==k)-a/(float)S)*((double)(xi[mu][C[i][x]]==l)-popularity[C[i][x]][l]);
          /*symmetric activity rule*/
          //J[i][x][k][l]+=((double)(xi[mu][i]==k)-sparsity_by_state[mu][k])*((double)(xi[mu][C[i][x]]==l)-sparsity_by_state[mu][l]);
          /*symmetric popularity rule*/
          //J[i][x][k][l]+=((double)(xi[mu][i]==k)-popularity[i][k])*((double)(xi[mu][C[i][x]]==l)-popularity[C[i][x]][l]);
          /*mixed rule*/
          //J[i][x][k][l]+=((double)(xi[mu][i]==k)-popularity[i][k])*((double)(xi[mu][C[i][x]]==l)-activity[mu][l]);
          /*mixed inverse rule*/
          //J[i][x][k][l]+=((double)(xi[mu][i]==k)-activity[mu][k])*((double)(xi[mu][C[i][x]]==l)-popularity[C[i][x]][l]);          
        } 
		//printf("%.1f\n",J[i][x][k][l]);	
		J[i][x][k][l]=J[i][x][k][l]/(a*(1.-a/S)*(double)Cm);
      }
    }
    //printf("x=%d\n", x);
  }
  //printf("i=%d", i);
}
printf("after Jixkl	\n");


}

/*-----------------------------------------------------------*/
/*initializes network*/
void initializing()
{
int i, j, k;
double rnumber, rrnumber;
/*set all states to zero, the exponential cue starts with first update*/
if (exp_cue == 1)
{
  for(i=0;i<N;i++)
  {
    units_receive_cue[i]=0.;
    if(drand48()<fraction_units_cued) {units_receive_cue[i]=1.;}
    for(j=0;j<S+1;j++)
    {
      s[i][j] = 0;
    }
  }
}
/*set all states of network into that of a given pattern and let network evolve*/
else if (full_cue == 1)
{
  //printf("this one\n");

  for(i=0;i<N;i++)
  {
	units_receive_cue[i]=1.;
    for(j=0;j<S+1;j++)
    {
      s[i][j] = 0;
    }
    s[i][xi[cue][i]] = 1.;
  }
  
}
/* to construct a cue from a given learned pattern such that SPARSITY IS CONSERVED and put net in this state PARTIAL CUE*/
else if (partial_cue == 1)
{

for(i=0;i<N;i++)
{
  units_receive_cue[i]=0.;
  if(drand48()<fraction_units_cued) {units_receive_cue[i]=1.;}
  
  /*first initialize all states to zero*/
  for(j=0;j<S+1;j++)
  {
    s[i][j] = 0;
  }
  /*if units_receive_cue = 1 put unit to state of pattern */
  if (units_receive_cue[i] == 1.)
  {
    s[i][xi[cue][i]] = 1.;
  }
  /*if units_receive_cue = 0 put unit with probability 1-a in inactive state and a/S in any active state*/
  else if (units_receive_cue[i] == 0)
  {
    rnumber = drand48();
    if (rnumber > a)
    {
      s[i][S] = 1;
    }
    else if (rnumber < a)
    {
      rrnumber = drand48();
      for (k=0; k<S; k++)
      {
	if (rrnumber >= k/(float)S and rrnumber < (k+1)/(float)S ) s[i][k] = 1;
      }
    } 
  }
}
}


//if (update_h == 1)
//{
  ///*compute overlap*/
  //compute_m(p);

  ///*initialize h*/
  ///*compute h*/
  //for(i=0; i<N; i++)
  //{
    ////printf("%d  ", i);
    //for(k=0;k<S;k++)
    //{
      //h[i][k]=0.;
      //for(x=0;x<number_connections[i];x++)
      //{
		//for(l=0;l<S;l++)
		//{
		  //h[i][k]+= J[i][x][k][l]*s[C[i][x]][l]; 
		//}
      //}
    //}
  //}
//}

//else if (update_r == 1)
//{
  ///*initialize s_0 and theta_0*/
  //for(i=0;i<N;i++)
  //{
    //for(k=0;k<S;k++)
    //{
    //s[i][k]=(-2*beta-2*exp(beta*U)-2*S+sqrt((2*beta+2*exp(beta*U)+2*S)*(2*beta+2*exp(beta*U)+2*S)+8*(-beta*beta-2*beta*S+2*beta*S*exp(beta*U))))/(2*(-beta*beta-2*beta*S+2*beta*S*exp(beta*U)));
    //}	
  //s[i][S]=1.-S*s[i][0];
  //theta[i][S]=1.-s[i][S];
  //}

  ///*initialize h and r*/
  ///*compute h*/
  //for(i=0; i<N; i++)
  //{
    //for(k=0;k<S;k++)
    //{
      //h[i][k]=0.;
      //for(x=0;x<number_connections[i];x++)
      //{
	//for(l=0;l<S;l++)
	//{
	  //h[i][k]+= J[i][x][k][l]*s[C[i][x]][l]; 
	//}
      //}
      //r[i][k]=h[i][k];
      //theta[i][k]=s[i][k];
    //}
  //}
//}


/*
H_interaction = 0.;
H_quad=0.;

for(v=0;v<N;v++)
{
  for(k=0;k<S;k++)
  {
    H_interaction+=-h[v][k]*s[v][k]/2. + U*s[v][k];
    for(x=0;x<Cm;x++)
    {
      for(l=0;l<S;l++)
      {
	H_quad+=-J[v][x][k][l]*s[C[v][x]][l]*s[v][k]/2.;
      }
    }
    H_quad+=U*s[v][k];
  }
}

fprintf(H_quad_file, "%d	%.4f\n", 0, H_quad);
fprintf(H_interaction_file, "%d	%.2f\n", 0, H_interaction);

fflush(H_quad_file);
fflush(H_interaction_file);
*/
}
/*------------------------------------------------------*/
/* update state */
void update_state(int trial, int i, int n)
{
int   k, x, l, kmax;
double Z, invZ, hmax, rmax, smax, rho, cue_field;	

if (update_h == 1)
  {

    /*compute h of neuron i*/
    hmax=0.;
    kmax=S+1;
    cue_field=g*exp(-(float)n/((float)tau));

    /*compute h of neuron i*/
    for(k=0;k<S;k++)
    {
      h[i][k]=0.;

      for(x=0;x<number_connections[i];x++)
      {
		for(l=0;l<S;l++)
		{
		  h[i][k]= h[i][k] + J[i][x][k][l]*s[C[i][x]][l];
		}
      }
    //   if(t >= 9.0){fprintf(field_file, "%f\n", h[i][k]);}
      //printf("%d\t", cue);
      h[i][k] = h[i][k] + (exp_cue==1)*cue_field*units_receive_cue[i]*(xi[cue][i]==k);
      /*find state receiving maximal field*/
      
      if(h[i][k]>hmax) 
      {
		  hmax=h[i][k]; 
		  kmax=k;
	  } 
    }
    /*------------------------ with discrete units -------------------------*/
    if (discrete == 1 and graded == 0)
    {
      /*check whether inactive gets more field with fixed input U*/
      if(U>hmax){kmax=S;} 
      else if(U==hmax) // if both threshold and largest field coincide, then pick one randomly
      {
		int random_number;
		random_number = drand48();
		if (random_number > 0.5){kmax=S;}
      }
      /*with a high temperature where units are essentially non graded*/
      for(k=0;k<S+1;k++)
      {
		s[i][k]=0;
      }
      s[i][kmax]=1;
    }
    /*------------------------ with graded units -------------------------*/
    else if (discrete == 0 and graded == 1)
    {
		/*compute normalization const*/
		Z=0.;
		for(k=0;k<S;k++)
		{ 
		  Z+=exp(beta*(h[i][k] - hmax)); 
		}
		Z+=exp(beta*(U-hmax)); 
		invZ=1./Z;

		/*compute each state of unit */
		for(k=0;k<S;k++)
		{
		  s[i][k]=invZ*exp(beta*(h[i][k]-hmax));
		}
		s[i][S]=invZ*exp(beta*(U-hmax));
    }
  }
  /*------------------------------------------------------------------------*/
  else if (update_r == 1)
  {  
    kmax=S+1;
    rmax=r[i][S];
    cue_field=g*exp(-(float)n/((float)tau));

    /*compute h of neuron i*/
    for(k=0;k<S;k++)
    {
      h[i][k]=0.;

      for(x=0;x<number_connections[i];x++)
      {
		for(l=0;l<S;l++)
		{
		  h[i][k] += J[i][x][k][l]*s[C[i][x]][l];
		}
      }
    //   if(t >= 9.0){fprintf(field_file, "%f\n", h[i][k]);}
      /*adds cue_field to h only if cue_i is in state k and global noise*/
      rho = 2*drand48()-1;
      h[i][k] += cue_field*units_receive_cue[i]*(xi[cue][i]==k) + nu*rho;
      /*update r and theta*/
      theta[i][k] += b2*(s[i][k]-theta[i][k]);
      r[i][k]+= b1*(h[i][k]-theta[i][k]-r[i][k]);
      /*find state receiving maximal field*/
      if(r[i][k]>rmax) {rmax=r[i][k]; kmax=k;} 
    }

    /*------------------------ with discrete units -------------------------*/
    if (discrete == 1 and graded == 0)
    {
      /*check whether inactive gets more field with fixed input U*/
      if(U>rmax){kmax=S;} 
      else if(U==rmax) // if both threshold and largest field coincide, then pick one randomly
      {
		int random_number;
		random_number = drand48();
		if (random_number > 0.5){kmax=S;}
      }
      /*with a high temperature where units are essentially non graded*/
      for(k=0;k<S+1;k++)
      {
		s[i][k]=0;
      }
      s[i][kmax]=1;
    }
    /*------------------------ with graded units -------------------------*/
    else if (discrete == 0 and graded == 1)
    {
    smax=0.0;
    kmax=0;
    for(k=0;k<S;k++)
    {
      if(s[i][k]>smax){kmax=k;}
    }

    theta[i][S] += b3*(1.-s[i][S]-r[i][S]);
    /*compute normalization const*/
    Z=0.;
    for(k=0;k<S;k++)
    { 
      Z+=exp(beta*(r[i][k]-rmax)); 
    }
    Z+=exp(beta*(theta[i][S]+U-rmax)); 
    invZ=1./Z;

    /*compute each state of unit */
    for(k=0;k<S;k++)
    {
      s[i][k]=invZ*exp(beta*(r[i][k]-rmax));
    }
    s[i][S]=invZ*exp(beta*(theta[i][S]+U-rmax));
    }
  }
/*------ compute H after every tempostampa number of neurons are updated---*/
/*
if((n%tempostampa)==0)
{

fprintf(state_file, "%.2f  ", t);
for(j=0;j<N;j++)
{
  count=0;
  for(k=0;k<S+1;k++)
  {
    if (s[j][k] > 0.95)
    {
      count++;
    }    
  } 
  if (count == 0)
  {
    printf("%d=	", j);
    fprintf(state_file, "%d  ", S+1);
  }
  else if (count != 0)
  {
    printf("%d=	", j);
    for(k=0;k<S+1;k++)
    {
      if (s[j][k] > 0.95)
      {
	fprintf(state_file, "%d  ", k);
      }    
    }
  }
}
fprintf(state_file, "\n");

for(k=0;k<S;k++)
{
  for(j=0;j<Cm;j++)
  {
    for(l=0;l<S;l++)
    {
      H_interaction=H_interaction-J[i][j][k][l]*(s[i][k]-sold[k])*s[C[i][j]][l];
    }
  }
  H_interaction=H_interaction+ U*(s[i][k]-sold[k]);
}


for(k=0;k<S;k++)
{
  for(l=0;l<S;l++)
  {
    for(j=0;j<Cm;j++)
    {
	    H_interaction=H_interaction-J[i][j][k][l]*(s[i][k]-spast[k])*s[C[i][j]][l];
    }
  }
  H_interaction=H_interaction+U*(s[i][k]-spast[k]);
}

fprintf(H_interaction_file, "%d	%d	%f	%.4f\n", Retr, trial, t, H_interaction);
fflush(H_interaction_file);

H_quad=0.;
for(v=0;v<N;v++)
{
  for(k=0;k<S;k++)
  {
    for(x=0;x<Cm;x++)
    {
      for(l=0;l<S;l++)
      {
	H_quad+=-J[v][x][k][l]*s[C[v][x]][l]*s[v][k]/2.;
      }
    }
    H_quad+=U*s[v][k];
  }
}
fprintf(H_quad_file, "%d	%d	%f	%.4f\n", Retr, trial, t, H_quad);
fflush(H_quad_file);

}
*/
}

/*--------------------------------------------------------*/
/* get memory */
void getmemory()
{
int i, j, x, z;

m = new double[p];
m_sorted = new double[p];
  
number_connections = new int[N];

sparsity_by_state = new double*[p];
for(i=0;i<p;i++)
  sparsity_by_state[i]=new double[S];
  
units_receive_cue = new double[N];

Permut= new int*[N];
for(i=0; i<N; i++)
  Permut[i]=new int[NumSet];

xi= new int*[p];
for(i=0; i<p; i++)
  xi[i]=new int[N];

s= new double*[N];
for(i=0; i<N; i++)
  s[i]=new double[S+1];

h= new double*[N];
for(i=0; i<N; i++)
  h[i]=new double[S];

r= new double*[N];
for(i=0; i<N; i++)
  r[i]=new double[S+1];

theta= new double*[N];
for(i=0; i<N; i++)
  theta[i]=new double[S+1];

J= new double***[N];
for(i=0; i<N; i++)
{
  J[i]=new double**[N];
  for(x=0; x<N; x++)
  {
    J[i][x]=new double*[S];
    for(z=0; z<S; z++)
    {
      J[i][x][z]=new double[S];
    }
  }
}

C= new int*[N];
for(i=0; i<N; i++)
  C[i]=new int[N];
  
C0= new double*[tot_num_trials];
for(i=0; i<tot_num_trials; i++)
  C0[i]=new double[tot_num_trials];

C1= new double*[tot_num_trials];
for(i=0; i<tot_num_trials; i++)
  C1[i]=new double[tot_num_trials];

C2= new double*[tot_num_trials];
for(i=0; i<tot_num_trials; i++)
  C2[i]=new double[tot_num_trials];

fin_state = new double **[tot_num_trials];

for(i=0;i<tot_num_trials;i++)
{
  fin_state[i] = new double *[N];
  for(j=0; j<N; j++)
  {
    fin_state[i][j] = new double [S+1];
  }
} 

}

//------------------------------------------------------------------------------//


void deletememory()
{
int i,x,z;

delete m;
delete m_sorted;

delete number_connections;

sparsity_by_state = new double*[p];
for(i=0;i<p;i++)
  sparsity_by_state[i] = new double[S];

delete units_receive_cue;

for(i=0; i<N; i++) delete (Permut[i]);
delete(Permut);

for(i=0; i<p; i++) delete (xi[i]);
delete(xi);

for(i=0; i<N; i++) delete (s[i]);
delete(s);

for(i=0; i<N; i++) delete (h[i]);
delete(h);

for(i=0; i<N; i++) delete (r[i]);
delete(r);

for(i=0; i<N; i++) delete (theta[i]);
delete(theta);

 for(i=0; i<N; i++)
 {
   for(x=0; x<N; x++)
   {
     for(z=0; z<S; z++)
     {
       delete J[i][x][z];
     }
   }
 }
 delete J;


for(i=0; i<N; i++) delete (C[i]);
delete(C);

for(i=0; i<tot_num_trials; i++) delete (C0[i]);
delete(C0);

for(i=0; i<tot_num_trials; i++) delete (C1[i]);
delete(C1);

for(i=0; i<tot_num_trials; i++) delete (C2[i]);
delete(C2);

for(i=0; i<tot_num_trials; i++) 
{
	for(x=0; x<N; x++)
	{
		delete fin_state[i][x];
	}
}
delete fin_state;

}

//------------------------------------------------------------------------------//


void SetUpTables() 
{
int item, jtem, info;
int fatto, kk;  

//srand48(time(0));
srand48(6937);
for(kk=0; kk<NumSet; kk++)
{
  item = 0; 
  while(item<N)
  {
    info = (int)((double)N*drand48());
    fatto=0;
    while(fatto==0)
    {
      fatto=1;
      for(jtem=0; jtem<item; jtem++)
      {
	if( Permut[jtem][kk] == info  )
	{
	  info=(info +1)-(int)((info+1)/N)*N;
	  jtem=item;
	  fatto=0;
	}
      }
    }
    Permut[item][kk]= info;
    item++;
  }
}
}
