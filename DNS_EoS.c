/* set EoS in terms of basic variable q */
/* Wolfgang Tichy, Tim Dietrich 5/2014 */

#include "sgrid.h"
#include "DNSdata.h"

/* functions in this file 
   DNS_polytrope_EoS_of_hm1
   DNS_polytrope_rho0_of_hm1
   DNS_polytrope_P_of_hm1
   DNS_polytrope_hm1_of_P
   DNS_pwp_init_...
 */

#define Power pow

/* global vars the describe polytrope pieces */
/* we switch from piece i to piece i+1 if rho0>=EoS_rho0[i+1] */
int     EoS_n_pieces = 0; /* number of pieces, i.e. how many different n */
double *EoS_rho0 = NULL;  /* EoS_rho0[0]=0 always,  EoS_rho0[i] is rho0 where we switch */
double *EoS_kappa = NULL; /* EoS_kappa[i] is kappa in piece i */
double *EoS_n = NULL;     /* EoS_n[i] is n in piece i */
double *EoS_k = NULL;     /* EoS_k[i] is k in piece i */
double *EoS_q = NULL;     /* EoS_q[i] = q(EoS_rho0[i]),  q=h-1 */
double *EoS_P = NULL;     /* EoS_P[i] = P(EoS_rho0[i]) */

/* allocate memory for global EoS variables. */
void alloc_EoS_globals(int pieces)
{
  int maxpieces = pieces + 3;
  EoS_rho0  = calloc(maxpieces, sizeof(*EoS_rho0));
  EoS_kappa = calloc(maxpieces, sizeof(*EoS_kappa));
  EoS_n     = calloc(maxpieces, sizeof(*EoS_n));
  EoS_k     = calloc(maxpieces, sizeof(*EoS_k));
  EoS_q     = calloc(maxpieces, sizeof(*EoS_q));
  EoS_P     = calloc(maxpieces, sizeof(*EoS_P));
}
/* free  memory for global EoS variables. */
void free_EoS_globals(void)
{
  free(EoS_rho0);
  free(EoS_kappa);
  free(EoS_n);
  free(EoS_k);
  free(EoS_q);
  free(EoS_P);
  EoS_rho0 = EoS_kappa = EoS_n = EoS_k = EoS_q = EoS_P = NULL;
}


/* for pwp */
/* find piece in which q=hm1 lies
   it is inside piece m, if EoS_q[m]<hm1<EoS_q[m+1] */
void DNS_select_polytrope_n_kappa_k_of_hm1(double hm1,
                                           double *n, double *kappa, double *k)
{
  int m;

  /* find m such that EoS_q[m]<hm1<EoS_q[m+1] */
  for(m=0; m<EoS_n_pieces-1; m++)
    if(hm1 < EoS_q[m+1]) break;
  /* Note, if EoS_n_pieces=1 we always get m=0 */

  *n     = EoS_n[m];
  *kappa = EoS_kappa[m];
  *k     = EoS_k[m];
}

/* for pwp select based on P */
void DNS_select_polytrope_n_kappa_k_of_P(double P,
                                         double *n, double *kappa, double *k)
{
  int m;

  /* find m such that EoS_P[m]<P<EoS_P[m+1] */
  for(m=0; m<EoS_n_pieces-1; m++)
    if(P < EoS_P[m+1]) break;
  /* Note, if EoS_n_pieces=1 we always get m=0 */

  *n     = EoS_n[m];
  *kappa = EoS_kappa[m];
  *k     = EoS_k[m];
}

/* EoS */
/* inputs: hm1 := h-1 , n=1/(Gamma-1) , kappa, constant k
   here: P = kappa rho0^Gamma 
         rho0 espilon = P/(Gamma-1) + k P^(1/Gamma)
   outputs: rest mass density rho0, pressure P, total energy density rhoE and polytropic index n*/
void DNS_polytrope_EoS_of_hm1(double hm1,
                              double *rho0, double *P, double *rhoE,
                              double *drho0dhm1)
{
  double hmk, n,kappa,k;

  DNS_select_polytrope_n_kappa_k_of_hm1(hm1, &n, &kappa, &k);
 	 
  hmk = hm1 + (1.0 - k);  /* hmk = hm1 + 1 - k */
  *rho0 = pow(hmk/(kappa*(n+1.0)), n);
  *P    = (*rho0) * hmk/(n+1.0);
//  *rhoE = (hm1+1.0) * (*rho0) - (*P);
  *rhoE = n*(*P) + k*(*rho0); 
  /* drho0/dhm1 = (n/(kappa*(n+1))) [hmk/(kappa*(n+1))]^{n-1} */
  if(*P>0.0)
    *drho0dhm1 = (*rho0)*(*rho0)*n/((*P)*(n+1.0));
  else /* limit for hmk->0 */
  {
    if(n<1.0) *drho0dhm1 = (n/(kappa*(n+1.0)))*pow(hmk/(kappa*(n+1.0)), n-1.0);
    else if(n==1.0) *drho0dhm1 = (n/(kappa*(n+1.0)));
    else            *drho0dhm1 = 0.0;
  }
}

/* get rho0 only from DNS_polytrope_EoS_of_hm1 */
double DNS_polytrope_rho0_of_hm1(double hm1)
{
  double rho0, dummy;
  DNS_polytrope_EoS_of_hm1(hm1, &rho0, &dummy, &dummy, &dummy);
  return rho0;
}
/* get P only from DNS_polytrope_EoS_of_hm1 */
double DNS_polytrope_P_of_hm1(double hm1)
{
  double rho0, P, dummy;
  DNS_polytrope_EoS_of_hm1(hm1, &rho0, &P, &dummy, &dummy);
  return P;
}

/* get hm1=h-1 from P */
/* P = kappa^{-n} [hmk/(n+1)]^{n+1}
   h-k = hmk = (n+1) kappa (P/kappa)^{1/(n+1)}  */
double DNS_polytrope_hm1_of_P(double P)
{
	
  double n, k, kappa, hmk; 

  DNS_select_polytrope_n_kappa_k_of_P(P, &n, &kappa, &k);

  hmk = (n+1.0) * kappa * pow(P/kappa, 1.0/(n+1.0));
  return hmk + (k-1.0);
}

/* get rho0, rhoE from P */
/* rho0 = pow(P/kappa, n/(n+1) 
   rhoE = n*P + k*rho0          */
void DNS_polytrope_rho0_rhoE_of_P(double P, double *rho0, double *rhoE)
{
  double n, k, kappa; 
  int m;

  DNS_select_polytrope_n_kappa_k_of_P(P, &n, &kappa, &k);

  *rho0 = pow(P/kappa, n/(n+1));
  *rhoE = n*P + k*(*rho0);
}

/* test if we get the correct n,kappa,k for different q and P */
void DNS_test_piecewise_n_kappa_k()
{
  int m;
  double hm1, P, n,kappa,k;

  printf("values of n,kappa,k for different test q and P:\n");
  for(m=-1; m<EoS_n_pieces; m++)
  {
    if(m<0) hm1 = P = -0.1;
    else
    {
      hm1 = 0.5*(EoS_q[m]+EoS_q[m+1]);
      P   = 0.5*(EoS_P[m]+EoS_P[m+1]);
    }
    DNS_select_polytrope_n_kappa_k_of_hm1(hm1, &n, &kappa, &k);
    printf("q=%e -> n=%e kappa=%e k=%e\n", hm1, n,kappa,k);
    DNS_select_polytrope_n_kappa_k_of_P(P, &n, &kappa, &k);
    printf("P=%e -> n=%e kappa=%e k=%e\n", P, n,kappa,k);
  }
}


/* initialize piecewise polytropes from files in Tim's format */
int DNS_pwp_init_file()
{
  int i;
  char *fname = Gets("DNSdata_EoS_file");
  char sline[1024];
  int counter = 0;

  printf(" reading PWP pars from file:\n    %s\n", fname); 

  // open file and read the data in
  FILE *ifp = fopen(fname, "r");  
  if (!ifp) errorexits(" could not open file %s", fname); 

  while (fgets(sline,256,ifp)) 
  {
    if (sline[0]=='#') {
      // skip comment
    }
    else if(counter==0){
      sscanf(sline,"%d", &EoS_n_pieces);
      free_EoS_globals();
      alloc_EoS_globals(EoS_n_pieces); /*get mem and set everything to zero */
      counter = counter + 1;
    }
    else if(!(counter==0)){
        sscanf(sline,"%le %le %le",&EoS_rho0[counter-1],&EoS_kappa[counter-1], &EoS_n[counter-1]); 
        counter = counter + 1;
    }
  }
  fclose (ifp);
 
  if(!((counter-1)==EoS_n_pieces)) errorexit("input EoS_file wrong");
  /* done with reading*/

  /* activate EoS_pwp*/
  EoS_pwp  = 1;
  EoS_k[0] = 1.;

  printf("DNS_pwp_init_file: using %d polytrope pieces \n", EoS_n_pieces);
  printf("rho0           q            P            kappa        n            k \n");

  /* compute the constant k, k = a+1 according to arXiv:0812.2163*/
  for(i=0;i<EoS_n_pieces;i++){
    if( i > 0) {    
      EoS_k[i] = ((EoS_k[i-1]*EoS_rho0[i]+EoS_n[i-1]
               * EoS_kappa[i-1]*pow(EoS_rho0[i],1.+1./EoS_n[i-1]))/EoS_rho0[i])
               - EoS_kappa[i]*EoS_n[i]*pow(EoS_rho0[i],1./EoS_n[i]);
    }
    EoS_q[i]  = EoS_k[i] + (EoS_n[i]+1)*EoS_kappa[i]*pow(EoS_rho0[i],1./EoS_n[i])-1.;
    EoS_P[i]  = EoS_kappa[i]*pow(EoS_rho0[i],1.+1./(EoS_n[i]));
    
    /* print rho, kappa, n, k in log-file*/ 
    printf("%e %e %e %e %e %e \n", 
    EoS_rho0[i], EoS_q[i], EoS_P[i], EoS_kappa[i], EoS_n[i], EoS_k[i]);
  }

  /* copy last column and set large end-value*/  
  EoS_rho0[EoS_n_pieces]  = EoS_rho0[(EoS_n_pieces-1)]  + 1.e10;
  EoS_q[EoS_n_pieces]     = EoS_q[(EoS_n_pieces-1)]     + 1.e10;
  EoS_P[EoS_n_pieces]     = EoS_P[(EoS_n_pieces-1)]     + 1.e10;
  EoS_kappa[EoS_n_pieces] = EoS_kappa[(EoS_n_pieces-1)];
  EoS_n[EoS_n_pieces]     = EoS_n[(EoS_n_pieces-1)];
  EoS_k[EoS_n_pieces]     = EoS_k[(EoS_n_pieces-1)];
  
  printf("%e %e %e %e %e %e \n", 
  EoS_rho0[EoS_n_pieces], EoS_q[EoS_n_pieces], EoS_P[EoS_n_pieces], 
  EoS_kappa[EoS_n_pieces], EoS_n[EoS_n_pieces], EoS_k[EoS_n_pieces]);
    

  /* Set parameters, 
     be careful that kappa is set correct, we do not recompute it, 
     when we use a file, but the DNSdata uses only kappa[0] 
     best: use notebook in /DNSdata/pwpfits */
  char par[1000];

  Setd("DNSdata_kappa",EoS_kappa[0]);

  sprintf(par,"%15.14f",EoS_n[0]);
  for(i=1;i<EoS_n_pieces;i++) sprintf(par,"%s %15.14f",par, EoS_n[i]);
  Sets("DNSdata_n",par);

  /* we do not use:  sprintf(par,"%15.14f", EoS_rho0[0]); 
     it's redundant. */
  sprintf(par,"%15.14f", EoS_rho0[1]);
  for(i=2;i<EoS_n_pieces;i++) sprintf(par,"%s %15.14f",par, EoS_rho0[i]);
  Sets("DNSdata_pwp_rho0",par);

  return 0;
}

/* initialize from parameters in sgrid parfile */
int DNS_pwp_init_parameter()
{
  int i;
  int rho0count, ncount;
  char *str, *s;

  /* count words in par DNSdata_n */
  str = strdup(Gets("DNSdata_n"));
  ncount = 0;
  for(s=strtok(str, " "); s!=NULL; s=strtok(NULL, " "))
    ncount++;
  free(str);

  /* initialize mem, and set everything to zero */
  free_EoS_globals();
  alloc_EoS_globals(ncount);

  /* activate EoS_pwp*/
  EoS_pwp  = 1;
  EoS_k[0] = 1.;

  /* loop over pars we read */
  for (i=1; i<3; i++)
  {
    int count, start;
    char str[200];
    char *par;
    double val;

    start = 0;
    count = 0;
    if(i==1) par = Gets("DNSdata_n");
    if(i==2) par = Gets("DNSdata_pwp_rho0");
    while(sscanf(par+start, "%s", str)==1)
    {
      start += strlen(str);
      if(par[start]==' ') start++;
      val = atof(str);
      if(i==1) EoS_n[count] = val;
      if(i==2)
      {
        /* if there is no leading 0 in DNSdata_pwp_rho0, skip to next EoS_rho0 */
        if(count==0 && val!=0.0) count++;
        EoS_rho0[count] = val;
      }
      count++;
    }
    if(i==1)  ncount    = count;
    if(i==2)  rho0count = count;
  }

  EoS_n_pieces = ncount;
  EoS_kappa[0] = Getd("DNSdata_kappa");

  printf("DNS_pwp_init_parameter: using %d polytrope pieces \n", EoS_n_pieces);
  printf("rho0           q            P            kappa        n            k \n");
 
  /* compute the constant k, k = a+1 according to arXiv:0812.2163*/
  for(i=0;i<EoS_n_pieces;i++)
  {
    if(i > 0)
    {  
      EoS_kappa[i] = EoS_kappa[i-1]*pow(EoS_rho0[i],(EoS_n[i]-EoS_n[i-1])/
                     (EoS_n[i-1]*EoS_n[i]));

      EoS_k[i] = ((EoS_k[i-1]*EoS_rho0[i]+EoS_n[i-1]
               * EoS_kappa[i-1]*pow(EoS_rho0[i],1.+1./EoS_n[i-1]))/EoS_rho0[i])
               - EoS_kappa[i]*EoS_n[i]*pow(EoS_rho0[i],1./EoS_n[i]);
    }
    EoS_q[i]  = EoS_k[i] + (EoS_n[i]+1)*EoS_kappa[i]*pow(EoS_rho0[i],1./EoS_n[i])-1.;
    EoS_P[i]  = EoS_kappa[i]*pow(EoS_rho0[i],1.+1./(EoS_n[i]));
    
    /* print rho, kappa, n, k in log-file*/ 
    printf("%e %e %e %e %e %e \n", 
    EoS_rho0[i], EoS_q[i], EoS_P[i], EoS_kappa[i], EoS_n[i], EoS_k[i]);
  }

  /* Copy last column and set large end-value.
     This is probably not needed at all. */  
  EoS_rho0[EoS_n_pieces]  = EoS_rho0[(EoS_n_pieces-1)]  + 1.e10;
  EoS_q[EoS_n_pieces]     = EoS_q[(EoS_n_pieces-1)]     + 1.e10;
  EoS_P[EoS_n_pieces]     = EoS_P[(EoS_n_pieces-1)]     + 1.e10;
  EoS_kappa[EoS_n_pieces] = EoS_kappa[(EoS_n_pieces-1)];
  EoS_n[EoS_n_pieces]     = EoS_n[(EoS_n_pieces-1)];
  EoS_k[EoS_n_pieces]     = EoS_k[(EoS_n_pieces-1)];
  /*
  printf("after last piece:\n");
  printf("%e %e %e %e %e %e\n",
  EoS_rho0[EoS_n_pieces], EoS_q[EoS_n_pieces], EoS_P[EoS_n_pieces], 
  EoS_kappa[EoS_n_pieces], EoS_n[EoS_n_pieces], EoS_k[EoS_n_pieces]);
  */
  DNS_test_piecewise_n_kappa_k();

  printf("\n");
  return 0;
}

/* intialize if we just use one regular polytrope */
int DNS_poly_init()
{
  EoS_pwp      = 0;
  EoS_n_pieces = 1;
  free_EoS_globals();
  alloc_EoS_globals(EoS_n_pieces);
  EoS_n[0]     = Getd("DNSdata_n");
  EoS_kappa[0] = Getd("DNSdata_kappa");
  EoS_k[0]     = 1;
  
  return 0;
}
