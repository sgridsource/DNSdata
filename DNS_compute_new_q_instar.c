/* DNS_compute_new_q_instar.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 21.2.2018 */
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




void DNS_compute_new_q_instar(tGrid *grid, int star, int iq)
{
int VwApprox1 = Getv("DNSdata_rotationstate1","VwApproximation");
int VwApprox2 = Getv("DNSdata_rotationstate2","VwApproximation");
int corot1 = VwApprox1 || Getv("DNSdata_rotationstate1","corotation");
int corot2 = VwApprox2 || Getv("DNSdata_rotationstate2","corotation");
int VwApprox, corot;
int qFromFields = Getv("DNSdata_new_q","FromFields");
double C1 = Getd("DNSdata_C1");
double C2 = Getd("DNSdata_C2");
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

int MATTRtouch  = (box->MATTR== TOUCH);



double *q = box->v[iq+0];
int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
int index_z = Ind("z");
double *z = box->v[index_z + 0];
int index_DNSdata_Psi = Ind("DNSdata_Psi");
double *Psi = box->v[index_DNSdata_Psi + 0];
int index_DNSdata_alphaP = Ind("DNSdata_alphaP");
double *alphaP = box->v[index_DNSdata_alphaP + 0];
int index_DNSdata_Bx = Ind("DNSdata_Bx");
double *B1 = box->v[index_DNSdata_Bx + 0];
double *B2 = box->v[index_DNSdata_Bx + 1];
double *B3 = box->v[index_DNSdata_Bx + 2];
int index_DNSdata_wBx = Ind("DNSdata_wBx");
double *wB1 = box->v[index_DNSdata_wBx + 0];
double *wB2 = box->v[index_DNSdata_wBx + 1];
double *wB3 = box->v[index_DNSdata_wBx + 2];
int index_DNSdata_Sigma = Ind("DNSdata_Sigma");
double *Sigma = box->v[index_DNSdata_Sigma + 0];
int index_DNSdata_Sigmax = Ind("DNSdata_Sigmax");
double *dSigma1 = box->v[index_DNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_DNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_DNSdata_Sigmax + 2];
int index_DNSdata_qgold = Ind("DNSdata_qgold");
double *qgold = box->v[index_DNSdata_qgold + 0];


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




if(box->SIDE != star) continue; 


FirstDerivsOf_S(box,index_DNSdata_Sigma,                   
                                   Ind("DNSdata_Sigmax"));


/* conditional */
if (MATTRtouch) {


int biin = bi-6; /* works only for my CubSph setup */                        

    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_Sigmax"), biin,bi);

    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_Sigmay"), biin,bi);

    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_Sigmaz"), biin,bi);

    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_wBx"), biin,bi);

    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_wBy"), biin,bi);

    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind("DNSdata_wBz"), biin,bi);}
/* if (MATTRtouch) */



VwApprox = corot = 0; 


if(box->SIDE==STAR1) { 


if(corot1)    corot = 1; 


if(VwApprox1) VwApprox = 1; 

CC
=
C1
;

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


} else { 


if(corot2)    corot = 1; 


if(VwApprox2) VwApprox = 1; 

CC
=
C2
;

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


forallpoints(box, ijk) { 

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
(uzero*(-CC + vR1*wDown1 + vR2*wDown2 + vR3*wDown3))/
  (1. + (betaVRDown1*vR1 + betaVRDown2*vR2 + betaVRDown3*vR3)/oouzerosqr)
;

q[ijk]
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


q[ijk]
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



/* conditional */
if (qFromFields) {

twoalpha2wdSigmapw
=
2.*alpha2*(w1*(wDown1 + dSigma1[ijk]) + w2*(wDown2 + dSigma2[ijk]) + 
    w3*(wDown3 + dSigma3[ijk]))
;

betadSigmaMinusCC
=
-CC + beta1*dSigma1[ijk] + beta2*dSigma2[ijk] + beta3*dSigma3[ijk]
;

bb
=
twoalpha2wdSigmapw + pow2(betadSigmaMinusCC)
;

L2
=
(0.5*(bb + Sqrt(Abs(pow2(bb) + pow2(twoalpha2wdSigmapw)))))/alpha2
;

h
=
Sqrt(Abs(L2 - (DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) - 
    (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) - 
    (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk])))
;


} else { /* if (!qFromFields) */

hOLD
=
1. + qgold[ijk]
;

hOLD2
=
pow2(hOLD)
;

LOLD2
=
hOLD2 + (DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) + 
  (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) + 
  (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk])
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
-(uzero*(CC + vR1*dSigma1[ijk] + vR2*dSigma2[ijk] + vR3*dSigma3[ijk]))
;

}
/* if (qFromFields) */


q[ijk]
=
-1. + h
;

}
/* if (qFromFields) */



} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* DNS_compute_new_q_instar.c */
/* nvars = 17, n* = 147,  n/ = 71,  n+ = 145, n = 363, O = 1 */
