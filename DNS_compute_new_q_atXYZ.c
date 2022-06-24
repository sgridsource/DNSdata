/* DNS_compute_new_q_atXYZ.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 23.6.2022 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "DNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    (sqrt((double) (x)))
#define Abs(x)     (fabs((double) (x)))
#define Log(x)     (log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




double DNS_compute_new_q_atXYZ(tGrid *grid, int bi, double X, double Y, double Z)
{
int VwApprox1 = Getv("DNSdata_rotationstate1","VwApproximation");
int VwApprox2 = Getv("DNSdata_rotationstate2","VwApproximation");
int corot1 = VwApprox1 || Getv("DNSdata_rotationstate1","corotation");
int corot2 = VwApprox2 || Getv("DNSdata_rotationstate2","corotation");
int VwApprox, corot;
int qFromFields = Getv("DNSdata_new_q","FromFields");
int FlipSgnL2   = Getv("DNSdata_new_q","FlipSignUnderRootOfL2Eqn");
double SgnL2    = -1 + 2*FlipSgnL2; //-1 or +1
int wB0outside  = Getv("DNSdata_wB_outside","0");
double C1 = Getd("DNSdata_C1");
double C2 = Getd("DNSdata_C2");
double Omega = Getd("DNSdata_Omega");
double xCM = Getd("DNSdata_x_CM");
double ecc = Getd("DNSdata_ecc");
double rdot = Getd("DNSdata_rdot");
double xmax1 = Getd("DNSdata_actual_xmax1");
double xmax2 = Getd("DNSdata_actual_xmax2");
double omegax1,omegay1,omegaz1, omegax2,omegay2,omegaz2;

tBox *box = grid->box[bi];
int ijk;

int MATTRtouch  = (box->MATTR== TOUCH);



int index_DNSdata_Psi = Ind("DNSdata_Psi");
double *DNSPsi = box->v[index_DNSdata_Psi + 0];
int index_DNSdata_alphaP = Ind("DNSdata_alphaP");
double *DNSalphaP = box->v[index_DNSdata_alphaP + 0];
int index_DNSdata_Bx = Ind("DNSdata_Bx");
double *DNSB1 = box->v[index_DNSdata_Bx + 0];
double *DNSB2 = box->v[index_DNSdata_Bx + 1];
double *DNSB3 = box->v[index_DNSdata_Bx + 2];
int index_DNSdata_q = Ind("DNSdata_q");
double *DNSq = box->v[index_DNSdata_q + 0];
int index_DNSdata_wBx = Ind("DNSdata_wBx");
double *DNSwB1 = box->v[index_DNSdata_wBx + 0];
double *DNSwB2 = box->v[index_DNSdata_wBx + 1];
double *DNSwB3 = box->v[index_DNSdata_wBx + 2];
int index_DNSdata_Sigma = Ind("DNSdata_Sigma");
double *DNSSigma = box->v[index_DNSdata_Sigma + 0];
int index_DNSdata_Sigmax = Ind("DNSdata_Sigmax");
double *DNSdSigma1 = box->v[index_DNSdata_Sigmax + 0];
double *DNSdSigma2 = box->v[index_DNSdata_Sigmax + 1];
double *DNSdSigma3 = box->v[index_DNSdata_Sigmax + 2];
int index_DNSdata_qgold = Ind("DNSdata_qgold");
double *DNSqgold = box->v[index_DNSdata_qgold + 0];
int index_DNSdata_temp4 = Ind("DNSdata_temp4");
double *temp4 = box->v[index_DNSdata_temp4 + 0];
double Psi, B1,B2,B3, alphaP;
double Sigma, dSigma1,dSigma2,dSigma3, wB1,wB2,wB3, x,y,z, qgold;


/* Jetzt geht's los! */
double alpha;
double alpha2;
double bb;
double beta1;
double beta2;
double beta3;
double betadSigmaMinusCC;
double betaVRDown1;
double betaVRDown2;
double betaVRDown3;
double CC;
double DSigmaUp1;
double DSigmaUp2;
double DSigmaUp3;
double h;
double hOLD;
double hOLD2;
double L2;
double LOLD2;
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
double q;
double rdotor;
double twoalpha2wdSigmapw;
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
double xC;
double xmax;
double xMxmax1;
double xMxmax2;
double xMxmax3;
double xrdotor1;
double xrdotor2;
double xrdotor3;




FirstDerivsOf_S(box,index_DNSdata_Sigma,                                    Ind("DNSdata_Sigmax")); 



/* conditional */
if (MATTRtouch && !wB0outside) {


int biin = bi-6; /* works only for my CubSph setup */                         
     copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_Sigmax"), biin,bi);
     copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_Sigmay"), biin,bi);
     copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_Sigmaz"), biin,bi);
     copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_wBx"), biin,bi);
     copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_wBy"), biin,bi);
     copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_wBz"), biin,bi);}
/* if (MATTRtouch && !wB0outside) */



if(box->x_of_X[1] != NULL) { 


x = box->x_of_X[1]((void *) box, -1, X,Y,Z); 


y = box->x_of_X[2]((void *) box, -1, X,Y,Z); 


z = box->x_of_X[3]((void *) box, -1, X,Y,Z); 


} else { 


x = X; 


y = Y; 


z = Z; 


} 


spec_Coeffs(box, DNSPsi, temp4); 


Psi = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSB1, temp4); 


B1 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSB2, temp4); 


B2 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSB3, temp4); 


B3 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSalphaP, temp4); 


alphaP = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSSigma, temp4); 


Sigma = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSdSigma1, temp4); 


dSigma1 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSdSigma2, temp4); 


dSigma2 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSdSigma3, temp4); 


dSigma3 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSwB1, temp4); 


wB1 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSwB2, temp4); 


wB2 = spec_interpolate(box, temp4, X,Y,Z); 


spec_Coeffs(box, DNSwB3, temp4); 


wB3 = spec_interpolate(box, temp4, X,Y,Z); 


VwApprox = corot = 0; 


if(box->SIDE==STAR1) { 

CC
=
C1
;

xmax
=
xmax1
;


if(corot1)    corot = 1; 


if(VwApprox1) VwApprox = 1; 


if(VwApprox) {                        omegax1 = Getd("DNSdata_omegax1");                        omegay1 = Getd("DNSdata_omegay1");                        omegaz1 = Getd("DNSdata_omegaz1"); } 

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


} else { 

CC
=
C2
;

xmax
=
xmax2
;


if(corot2)    corot = 1; 


if(VwApprox2) VwApprox = 1; 


if(VwApprox) {                        omegax2 = Getd("DNSdata_omegax2");                        omegay2 = Getd("DNSdata_omegay2");                        omegaz2 = Getd("DNSdata_omegaz2"); } 

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

xC
=
xCM + ecc*(-xCM + xmax)
;

OmegaCrossR1
=
-(Omega*y)
;

OmegaCrossR2
=
Omega*(x - xC)
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
rdotor*(x - xCM)
;

xrdotor2
=
rdotor*y
;

xrdotor3
=
rdotor*z
;

beta1
=
B1 + OmegaCrossR1 + xrdotor1
;

beta2
=
B2 + OmegaCrossR2 + xrdotor2
;

beta3
=
B3 + OmegaCrossR3 + xrdotor3
;

alpha
=
alphaP/Psi
;

alpha2
=
pow2(alpha)
;

Psi2
=
pow2(Psi)
;

Psi4
=
pow2(Psi2)
;

Psim4
=
1/Psi4
;

Psim2
=
Psi2*Psim4
;



/* conditional */
if (corot) {



/* conditional */
if (VwApprox) {

xMxmax1
=
x - xmax
;

xMxmax2
=
y
;

xMxmax3
=
z
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

betaVRDown1
=
Psi4*(beta1 + vR1)
;

betaVRDown2
=
Psi4*(beta2 + vR2)
;

betaVRDown3
=
Psi4*(beta3 + vR3)
;

oouzerosqr
=
alpha2 - betaVRDown1*(beta1 + vR1) - betaVRDown2*(beta2 + vR2) - 
  betaVRDown3*(beta3 + vR3)
;



/* conditional */
if (oouzerosqr == 0) {

oouzerosqr
=
-1.
;

}
/* if (oouzerosqr == 0) */




/* conditional */
if (oouzerosqr < 0) {

uzero
=
-1.
;


} else { /* if (!oouzerosqr < 0) */

uzero
=
Sqrt(1/oouzerosqr)
;

}
/* if (oouzerosqr < 0) */


wBDown1
=
wB1
;

wBDown2
=
wB2
;

wBDown3
=
wB3
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
(uzero*(-CC + vR1*wDown1 + vR2*wDown2 + vR3*wDown3))/
  (1. + (betaVRDown1*vR1 + betaVRDown2*vR2 + betaVRDown3*vR3)/oouzerosqr)
;

q
=
-1. + h
;


} else { /* if (!oouzerosqr < 0) */

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

oouzerosqr
=
alpha2 - Psi4*(2.*(beta1*vR1 + beta2*vR2 + beta3*vR3) + pow2(beta1) + 
     pow2(beta2) + pow2(beta3) + pow2(vR1) + pow2(vR2) + pow2(vR3))
;



/* conditional */
if (oouzerosqr == 0) {

oouzerosqr
=
-1.
;

}
/* if (oouzerosqr == 0) */




/* conditional */
if (oouzerosqr < 0) {

uzero
=
-1.
;


} else { /* if (!oouzerosqr < 0) */

uzero
=
Sqrt(1/oouzerosqr)
;

}
/* if (oouzerosqr < 0) */


q
=
-1. - CC*uzero
;

}
/* if (oouzerosqr < 0) */



} else { /* if (!oouzerosqr < 0) */

Psim6
=
Psim2*Psim4
;

DSigmaUp1
=
dSigma1*Psim4
;

DSigmaUp2
=
dSigma2*Psim4
;

DSigmaUp3
=
dSigma3*Psim4
;

w1
=
Psim6*wB1
;

w2
=
Psim6*wB2
;

w3
=
Psim6*wB3
;

wBDown1
=
wB1
;

wBDown2
=
wB2
;

wBDown3
=
wB3
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



/* conditional */
if (qFromFields) {

twoalpha2wdSigmapw
=
2.*alpha2*(w1*(dSigma1 + wDown1) + w2*(dSigma2 + wDown2) + 
    w3*(dSigma3 + wDown3))
;

betadSigmaMinusCC
=
-CC + beta1*dSigma1 + beta2*dSigma2 + beta3*dSigma3
;

bb
=
twoalpha2wdSigmapw + pow2(betadSigmaMinusCC)
;

L2
=
(0.5*(bb + Sqrt(Abs(pow2(bb) + SgnL2*pow2(twoalpha2wdSigmapw)))))/alpha2
;

h
=
Sqrt(Abs(L2 - (DSigmaUp1 + w1)*(dSigma1 + wDown1) - 
    (DSigmaUp2 + w2)*(dSigma2 + wDown2) - (DSigmaUp3 + w3)*(dSigma3 + wDown3)\
))
;


} else { /* if (!qFromFields) */


spec_Coeffs(box, DNSqgold, temp4); 


qgold = spec_interpolate(box, temp4, X,Y,Z); 

hOLD
=
1. + qgold
;

hOLD2
=
pow2(hOLD)
;

LOLD2
=
hOLD2 + (DSigmaUp1 + w1)*(dSigma1 + wDown1) + 
  (DSigmaUp2 + w2)*(dSigma2 + wDown2) + (DSigmaUp3 + w3)*(dSigma3 + wDown3)
;

uzerosqr
=
LOLD2/(alpha2*hOLD2)
;

uzero
=
sqrt(uzerosqr)
;

vR1
=
-beta1 + (DSigmaUp1 + w1)/(hOLD*uzero)
;

vR2
=
-beta2 + (DSigmaUp2 + w2)/(hOLD*uzero)
;

vR3
=
-beta3 + (DSigmaUp3 + w3)/(hOLD*uzero)
;

h
=
-(uzero*(CC + dSigma1*vR1 + dSigma2*vR2 + dSigma3*vR3))
;

}
/* if (qFromFields) */


q
=
-1. + h
;

}
/* if (qFromFields) */



/* end of computation */ 

return q;
}  /* end of function */

/* DNS_compute_new_q_atXYZ.c */
/* nvars = 14, n* = 148,  n/ = 71,  n+ = 148, n = 367, O = 1 */
