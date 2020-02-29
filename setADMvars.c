/* setADMvars.c */
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


void setADMvars(tGrid *grid)
{
int VwApprox1 = Getv("DNSdata_rotationstate1","VwApproximation");
int VwApprox2 = Getv("DNSdata_rotationstate2","VwApproximation");
int corot1 = VwApprox1 || Getv("DNSdata_rotationstate1","corotation");
int corot2 = VwApprox2 || Getv("DNSdata_rotationstate2","corotation");
int VwApprox, corot;
int shiftver1 = Getv("DNSdata_ADMshift","B^i+phidotphi^i+rdotor0r^i");
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
int index_DNSdata_VRx = Ind("DNSdata_VRx");
double *VR1 = box->v[index_DNSdata_VRx + 0];
double *VR2 = box->v[index_DNSdata_VRx + 1];
double *VR3 = box->v[index_DNSdata_VRx + 2];
int index_DNSdata_rhobar = Ind("DNSdata_rhobar");
double *rhobar = box->v[index_DNSdata_rhobar + 0];
int index_DNSdata_Bxx = Ind("DNSdata_Bxx");
double *dB11 = box->v[index_DNSdata_Bxx + 0];
double *dB12 = box->v[index_DNSdata_Bxx + 1];
double *dB13 = box->v[index_DNSdata_Bxx + 2];
double *dB21 = box->v[index_DNSdata_Bxx + 3];
double *dB22 = box->v[index_DNSdata_Bxx + 4];
double *dB23 = box->v[index_DNSdata_Bxx + 5];
double *dB31 = box->v[index_DNSdata_Bxx + 6];
double *dB32 = box->v[index_DNSdata_Bxx + 7];
double *dB33 = box->v[index_DNSdata_Bxx + 8];
int index_DNSdata_Sigmax = Ind("DNSdata_Sigmax");
double *dSigma1 = box->v[index_DNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_DNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_DNSdata_Sigmax + 2];
int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
int index_z = Ind("z");
double *z = box->v[index_z + 0];
int index_gxx = Ind("gxx");
double *g11 = box->v[index_gxx + 0];
double *g12 = box->v[index_gxx + 1];
double *g13 = box->v[index_gxx + 2];
double *g22 = box->v[index_gxx + 3];
double *g23 = box->v[index_gxx + 4];
double *g33 = box->v[index_gxx + 5];
int index_psi = Ind("psi");
double *psi = box->v[index_psi + 0];
int index_Kxx = Ind("Kxx");
double *K11 = box->v[index_Kxx + 0];
double *K12 = box->v[index_Kxx + 1];
double *K13 = box->v[index_Kxx + 2];
double *K22 = box->v[index_Kxx + 3];
double *K23 = box->v[index_Kxx + 4];
double *K33 = box->v[index_Kxx + 5];
int index_alpha = Ind("alpha");
double *alpha = box->v[index_alpha + 0];
int index_betax = Ind("betax");
double *beta1 = box->v[index_betax + 0];
double *beta2 = box->v[index_betax + 1];
double *beta3 = box->v[index_betax + 2];
int index_rho = Ind("rho");
double *rho = box->v[index_rho + 0];
int index_jx = Ind("jx");
double *jdo1 = box->v[index_jx + 0];
double *jdo2 = box->v[index_jx + 1];
double *jdo3 = box->v[index_jx + 2];
int index_Sxx = Ind("Sxx");
double *Sdo11 = box->v[index_Sxx + 0];
double *Sdo12 = box->v[index_Sxx + 1];
double *Sdo13 = box->v[index_Sxx + 2];
double *Sdo22 = box->v[index_Sxx + 3];
double *Sdo23 = box->v[index_Sxx + 4];
double *Sdo33 = box->v[index_Sxx + 5];


/* Jetzt geht's los! */
double alpha2;
double DSigmaUp1;
double DSigmaUp2;
double DSigmaUp3;
double gdB;
double h;
double h2;
double jup1;
double jup2;
double jup3;
double LB11;
double LB12;
double LB13;
double LB22;
double LB23;
double LB33;
double LBdo11;
double LBdo12;
double LBdo13;
double LBdo21;
double LBdo22;
double LBdo23;
double LBdo31;
double LBdo32;
double LBdo33;
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
double vRplusbetado1;
double vRplusbetado2;
double vRplusbetado3;
double w1;
double w2;
double w3;
double wBDown1;
double wBDown2;
double wBDown3;
double wDown1;
double wDown2;
double wDown3;
double xC;
double xmax;
double xMxmax1;
double xMxmax2;
double xMxmax3;
double xrdotor1;
double xrdotor2;
double xrdotor3;




FirstDerivsOf_Sa(box, Ind("DNSdata_Bx"),       Ind("DNSdata_Bxx")); 


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

alpha[ijk]
=
alphaP[ijk]/Psi[ijk]
;

alpha2
=
pow2(alpha[ijk])
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

beta1[ijk]
=
OmegaCrossR1 + xrdotor1 + B1[ijk]
;

beta2[ijk]
=
OmegaCrossR2 + xrdotor2 + B2[ijk]
;

beta3[ijk]
=
OmegaCrossR3 + xrdotor3 + B3[ijk]
;

gdB
=
dB11[ijk] + dB22[ijk] + dB33[ijk]
;

LB11
=
-0.66666666666666666667*gdB + 2.*dB11[ijk]
;

LB12
=
dB12[ijk] + dB21[ijk]
;

LB13
=
dB13[ijk] + dB31[ijk]
;

LB22
=
-0.66666666666666666667*gdB + 2.*dB22[ijk]
;

LB23
=
dB23[ijk] + dB32[ijk]
;

LB33
=
-0.66666666666666666667*gdB + 2.*dB33[ijk]
;

LBdo11
=
LB11
;

LBdo12
=
LB12
;

LBdo13
=
LB13
;

LBdo21
=
LB12
;

LBdo22
=
LB22
;

LBdo23
=
LB23
;

LBdo31
=
LB13
;

LBdo32
=
LB23
;

LBdo33
=
LB33
;

psi[ijk]
=
1.
;

g11[ijk]
=
Psi4
;

g12[ijk]
=
0
;

g13[ijk]
=
0
;

g22[ijk]
=
Psi4
;

g23[ijk]
=
0
;

g33[ijk]
=
Psi4
;

K11[ijk]
=
(0.5*LBdo11*Psi4)/alpha[ijk]
;

K12[ijk]
=
(0.5*LBdo12*Psi4)/alpha[ijk]
;

K12[ijk]
=
(0.5*LBdo21*Psi4)/alpha[ijk]
;

K13[ijk]
=
(0.5*LBdo13*Psi4)/alpha[ijk]
;

K13[ijk]
=
(0.5*LBdo31*Psi4)/alpha[ijk]
;

K22[ijk]
=
(0.5*LBdo22*Psi4)/alpha[ijk]
;

K23[ijk]
=
(0.5*LBdo23*Psi4)/alpha[ijk]
;

K23[ijk]
=
(0.5*LBdo32*Psi4)/alpha[ijk]
;

K33[ijk]
=
(0.5*LBdo33*Psi4)/alpha[ijk]
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
alpha2 - Psi4*(2.*(vR1*beta1[ijk] + vR2*beta2[ijk] + vR3*beta3[ijk]) + 
     pow2(vR1) + pow2(vR2) + pow2(vR3) + pow2(beta1[ijk]) + 
     pow2(beta2[ijk]) + pow2(beta3[ijk]))
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
(DSigmaUp1 + w1)/(h*uzero) - beta1[ijk]
;

vR2
=
(DSigmaUp2 + w2)/(h*uzero) - beta2[ijk]
;

vR3
=
(DSigmaUp3 + w3)/(h*uzero) - beta3[ijk]
;

}
/* if (VwApprox) */


VR1[ijk]
=
vR1
;

VR2[ijk]
=
vR2
;

VR3[ijk]
=
vR3
;


EoS->vars_from_hm1(q[ijk], &rho0, &P, &rhoE, &drho0dhm1); 



/* conditional */
if (q[ijk] == 0) {

uzerosqr
=
0
;

}
/* if (q[ijk] == 0) */


rho[ijk]
=
alpha2*rhoE*uzerosqr + P*(-1. + alpha2*uzerosqr)
;

jup1
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR1 + beta1[ijk])
;

jup2
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR2 + beta2[ijk])
;

jup3
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR3 + beta3[ijk])
;

jdo1[ijk]
=
jup1*g11[ijk] + jup2*g12[ijk] + jup3*g13[ijk]
;

jdo2[ijk]
=
jup1*g12[ijk] + jup2*g22[ijk] + jup3*g23[ijk]
;

jdo3[ijk]
=
jup1*g13[ijk] + jup2*g23[ijk] + jup3*g33[ijk]
;

vRplusbetado1
=
(vR1 + beta1[ijk])*g11[ijk] + (vR2 + beta2[ijk])*g12[ijk] + 
  (vR3 + beta3[ijk])*g13[ijk]
;

vRplusbetado2
=
(vR1 + beta1[ijk])*g12[ijk] + (vR2 + beta2[ijk])*g22[ijk] + 
  (vR3 + beta3[ijk])*g23[ijk]
;

vRplusbetado3
=
(vR1 + beta1[ijk])*g13[ijk] + (vR2 + beta2[ijk])*g23[ijk] + 
  (vR3 + beta3[ijk])*g33[ijk]
;

Sdo11[ijk]
=
P*g11[ijk] + (P + rhoE)*uzerosqr*pow2(vRplusbetado1)
;

Sdo12[ijk]
=
(P + rhoE)*uzerosqr*vRplusbetado1*vRplusbetado2 + P*g12[ijk]
;

Sdo13[ijk]
=
(P + rhoE)*uzerosqr*vRplusbetado1*vRplusbetado3 + P*g13[ijk]
;

Sdo22[ijk]
=
P*g22[ijk] + (P + rhoE)*uzerosqr*pow2(vRplusbetado2)
;

Sdo23[ijk]
=
(P + rhoE)*uzerosqr*vRplusbetado2*vRplusbetado3 + P*g23[ijk]
;

Sdo33[ijk]
=
P*g33[ijk] + (P + rhoE)*uzerosqr*pow2(vRplusbetado3)
;



/* conditional */
if (shiftver1) {

OmegaCrossR1
=
(-1. + ecc)*Omega*y[ijk]
;

OmegaCrossR2
=
Omega*((-1. + ecc)*xCM + x[ijk] - ecc*x[ijk])
;

OmegaCrossR3
=
0
;

beta1[ijk]
=
OmegaCrossR1 + xrdotor1 + B1[ijk]
;

beta2[ijk]
=
OmegaCrossR2 + xrdotor2 + B2[ijk]
;

beta3[ijk]
=
OmegaCrossR3 + xrdotor3 + B3[ijk]
;

}
/* if (shiftver1) */


rhobar[ijk]
=
rho[ijk]*pow2(Psi4)
;


} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* setADMvars.c */
/* nvars = 56, n* = 218,  n/ = 61,  n+ = 234, n = 513, O = 1 */
