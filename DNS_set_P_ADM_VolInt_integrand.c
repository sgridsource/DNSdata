/* DNS_set_P_ADM_VolInt_integrand.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 28.2.2020 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "DNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))


extern tEoS EoS[1];


void DNS_set_P_ADM_VolInt_integrand_Om_xcm(tGrid *grid,                                                  int iIntegx, int iIntegy, int iIntegz,                                                  double Om, double xcm)
{
int VwApprox1 = Getv("DNSdata_rotationstate1","VwApproximation");
int VwApprox2 = Getv("DNSdata_rotationstate2","VwApproximation");
int corot1 = VwApprox1 || Getv("DNSdata_rotationstate1","corotation");
int corot2 = VwApprox2 || Getv("DNSdata_rotationstate2","corotation");
int VwApprox, corot;
double Omega = Getd("DNSdata_Omega");
double xCM = Getd("DNSdata_x_CM");
double ecc = Getd("DNSdata_ecc");
double rdot = Getd("DNSdata_rdot");
double xmax1 = Getd("DNSdata_actual_xmax1");
double xmax2 = Getd("DNSdata_actual_xmax2");
double omegax1 = Getd("DNSdata_omegax1");
double omegay1 = Getd("DNSdata_omegay1");
double omegaz1 = Getd("DNSdata_omegaz1");
double omegax2 = Getd("DNSdata_omegax2");
double omegay2 = Getd("DNSdata_omegay2");
double omegaz2 = Getd("DNSdata_omegaz2");

int bi;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;

int isSTAR1 = (box->SIDE == STAR1);



int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
int index_z = Ind("z");
double *z = box->v[index_z + 0];
int index_DNSdata_Psi = Ind("DNSdata_Psi");
double *Psi = box->v[index_DNSdata_Psi + 0];
int index_DNSdata_Bx = Ind("DNSdata_Bx");
double *B1 = box->v[index_DNSdata_Bx + 0];
double *B2 = box->v[index_DNSdata_Bx + 1];
double *B3 = box->v[index_DNSdata_Bx + 2];
int index_DNSdata_alphaP = Ind("DNSdata_alphaP");
double *alphaP = box->v[index_DNSdata_alphaP + 0];
int index_DNSdata_Sigma = Ind("DNSdata_Sigma");
double *Sigma = box->v[index_DNSdata_Sigma + 0];
int index_DNSdata_q = Ind("DNSdata_q");
double *q = box->v[index_DNSdata_q + 0];
int index_DNSdata_wBx = Ind("DNSdata_wBx");
double *wB1 = box->v[index_DNSdata_wBx + 0];
double *wB2 = box->v[index_DNSdata_wBx + 1];
double *wB3 = box->v[index_DNSdata_wBx + 2];
int index_DNSdata_Sigmax = Ind("DNSdata_Sigmax");
double *dSigma1 = box->v[index_DNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_DNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_DNSdata_Sigmax + 2];
double *Integx = box->v[iIntegx+0];
double *Integy = box->v[iIntegy+0];
double *Integz = box->v[iIntegz+0];


/* Jetzt geht's los! */
double alpha;
double alpha2;
double bet1;
double bet2;
double bet3;
double beta1;
double beta2;
double beta3;
double DSigmaUp1;
double DSigmaUp2;
double DSigmaUp3;
double h;
double h2;
double jup1;
double jup2;
double jup3;
double OmCrossR1;
double OmCrossR2;
double OmCrossR3;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double omegMOmeg1;
double omegMOmeg2;
double omegMOmeg3;
double oouzerosqr;
double Psi2;
double Psi4;
double Psim2;
double Psim4;
double Psim6;
double rdotor;
double uzero;
double uzerosqr;
double vR1;
double vR2;
double vR3;
double w1;
double w2;
double w3;
double wBDown1;
double wBDown2;
double wBDown3;
double wDown1;
double wDown2;
double wDown3;
double xc;
double xC;
double xmax;
double xMxmax1;
double xMxmax2;
double xMxmax3;
double xrd1;
double xrd2;
double xrd3;
double xrdotor1;
double xrdotor2;
double xrdotor3;




FirstDerivsOf_S(box, Ind("DNSdata_Sigma"),     Ind("DNSdata_Sigmax")); 


VwApprox = corot = 0; 



/* conditional */
if (isSTAR1) {


if(corot1)    corot = 1; 


if(VwApprox1) VwApprox = 1; 

xmax
=
xmax1
;

omegMOmeg1
=
omegax1
;

omegMOmeg2
=
omegay1
;

omegMOmeg3
=
-Omega + omegaz1
;


} else { /* if (!isSTAR1) */


if(corot2)    corot = 1; 


if(VwApprox2) VwApprox = 1; 

xmax
=
xmax2
;

omegMOmeg1
=
omegax2
;

omegMOmeg2
=
omegay2
;

omegMOmeg3
=
-Omega + omegaz2
;

}
/* if (isSTAR1) */



forallpoints(box, ijk) { 


double rho0, P, rhoE, drho0dhm1; /* declare them */ 

alpha
=
alphaP[ijk]/Psi[ijk]
;

alpha2
=
pow2(alpha)
;

Psi2
=
pow2(Psi[ijk])
;

Psi4
=
pow2(Psi2)
;

xC
=
xCM + ecc*(-xCM + xmax)
;

OmegaCrossR1
=
-(Omega*y[ijk])
;

OmegaCrossR2
=
Omega*(-xC + x[ijk])
;

OmegaCrossR3
=
0
;

xc
=
xcm + ecc*(-xcm + xmax)
;

OmCrossR1
=
-(Om*y[ijk])
;

OmCrossR2
=
Om*(-xc + x[ijk])
;

OmCrossR3
=
0
;

rdotor
=
rdot/(xmax1 - xmax2)
;

xrdotor1
=
rdotor*(-xCM + x[ijk])
;

xrdotor2
=
rdotor*y[ijk]
;

xrdotor3
=
rdotor*z[ijk]
;

rdotor
=
rdot/(xmax1 - xmax2)
;

xrd1
=
rdotor*(-xcm + x[ijk])
;

xrd2
=
rdotor*y[ijk]
;

xrd3
=
rdotor*z[ijk]
;

beta1
=
OmegaCrossR1 + xrdotor1 + B1[ijk]
;

beta2
=
OmegaCrossR2 + xrdotor2 + B2[ijk]
;

beta3
=
OmegaCrossR3 + xrdotor3 + B3[ijk]
;

bet1
=
OmCrossR1 + xrd1 + B1[ijk]
;

bet2
=
OmCrossR2 + xrd2 + B2[ijk]
;

bet3
=
OmCrossR3 + xrd3 + B3[ijk]
;



/* conditional */
if (corot) {



/* conditional */
if (VwApprox) {

xMxmax1
=
-xmax + x[ijk]
;

xMxmax2
=
y[ijk]
;

xMxmax3
=
z[ijk]
;

vR1
=
-(omegMOmeg3*xMxmax2) + omegMOmeg2*xMxmax3
;

vR2
=
omegMOmeg3*xMxmax1 - omegMOmeg1*xMxmax3
;

vR3
=
-(omegMOmeg2*xMxmax1) + omegMOmeg1*xMxmax2
;


} else { /* if (!VwApprox) */

vR1
=
0
;

vR2
=
0
;

vR3
=
0
;

}
/* if (VwApprox) */


oouzerosqr
=
alpha2 - Psi4*(2.*(beta1*vR1 + beta2*vR2 + beta3*vR3) + pow2(beta1) + 
     pow2(beta2) + pow2(beta3) + pow2(vR1) + pow2(vR2) + pow2(vR3))
;

uzerosqr
=
1./oouzerosqr
;


} else { /* if (!VwApprox) */

Psim2
=
1/Psi2
;

Psim4
=
pow2(Psim2)
;

Psim6
=
Psim2*Psim4
;

DSigmaUp1
=
Psim4*dSigma1[ijk]
;

DSigmaUp2
=
Psim4*dSigma2[ijk]
;

DSigmaUp3
=
Psim4*dSigma3[ijk]
;

w1
=
Psim6*wB1[ijk]
;

w2
=
Psim6*wB2[ijk]
;

w3
=
Psim6*wB3[ijk]
;

wBDown1
=
wB1[ijk]
;

wBDown2
=
wB2[ijk]
;

wBDown3
=
wB3[ijk]
;

wDown1
=
Psim2*wBDown1
;

wDown2
=
Psim2*wBDown2
;

wDown3
=
Psim2*wBDown3
;

h
=
1. + q[ijk]
;

h2
=
pow2(h)
;

uzerosqr
=
1/alpha2 + ((DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) + 
     (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) + 
     (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk]))/(alpha2*h2)
;

uzero
=
sqrt(uzerosqr)
;

vR1
=
-beta1 + (DSigmaUp1 + w1)/(h*uzero)
;

vR2
=
-beta2 + (DSigmaUp2 + w2)/(h*uzero)
;

vR3
=
-beta3 + (DSigmaUp3 + w3)/(h*uzero)
;

}
/* if (VwApprox) */



EoS->vars_from_hm1(q[ijk], &rho0, &P, &rhoE, &drho0dhm1); 



/* conditional */
if (q[ijk] == 0) {

uzerosqr
=
0
;

}
/* if (q[ijk] == 0) */


jup1
=
alpha*(P + rhoE)*uzerosqr*(bet1 + vR1)
;

jup2
=
alpha*(P + rhoE)*uzerosqr*(bet2 + vR2)
;

jup3
=
alpha*(P + rhoE)*uzerosqr*(bet3 + vR3)
;

Integx[ijk]
=
jup1*Psi2*pow2(Psi4)
;

Integy[ijk]
=
jup2*Psi2*pow2(Psi4)
;

Integz[ijk]
=
jup3*Psi2*pow2(Psi4)
;


} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

void DNS_set_P_ADM_VolInt_integrand(tGrid *grid, int iIntegx, int iIntegy, int iIntegz)
{
double Omega = Getd("DNSdata_Omega");
double xCM = Getd("DNSdata_x_CM");
DNS_set_P_ADM_VolInt_integrand_Om_xcm(grid, iIntegx, iIntegy, iIntegz ,Omega, xCM);
}  /* end of 2nd function */

/* DNS_set_P_ADM_VolInt_integrand.c */
/* nvars = 19, n* = 119,  n/ = 51,  n+ = 117, n = 287, O = 1 */
