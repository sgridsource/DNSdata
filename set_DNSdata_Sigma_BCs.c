/* set_DNSdata_Sigma_BCs.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 22.2.2018 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "DNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     (log((double) (x)))
#define Sqrt(x)    (sqrt((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void set_DNSdata_Sigma_BC(tVarList *vlFu, tVarList *vlu,  
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, 		   int nonlin)
{
int VwApprox1 = Getv("DNSdata_rotationstate1","VwApproximation");
int VwApprox2 = Getv("DNSdata_rotationstate2","VwApproximation");
int corot1 = VwApprox1 || Getv("DNSdata_rotationstate1","corotation");
int corot2 = VwApprox2 || Getv("DNSdata_rotationstate2","corotation");
int dqFromqg = Getv("DNSdata_q_derivs","dqg");
int dQFromdlam = Getv("DNSdata_drho0_inBC","dlam");
int SigmaZeroAtPoint = Getv("DNSdata_Sigma_surface_BCs","ZeroAtPoint");
int AddInnerVolIntToBC = Getv("DNSdata_Sigma_surface_BCs","AddInnerVolIntToBC");
int InnerVolIntZero = Getv("DNSdata_Sigma_surface_BCs","InnerVolIntZero");
int AddInnerSumToBC = Getv("DNSdata_Sigma_surface_BCs","AddInnerSumToBC");
int InnerSumZero = Getv("DNSdata_Sigma_surface_BCs","InnerSumZero");
int SigmaZeroInOuterBoxes = Getv("DNSdata_Sigma_surface_BCs","ZeroInOuterBoxes");
int noBCs = Getv("DNSdata_Sigma_surface_BCs","none");
int KeepInnerSigma = Getv("DNSdata_KeepInnerSigma","yes");
int ImposeActualBC = !Getv("DNSdata_Sigma_surface_BCs","EllEqn");
double Omega = Getd("DNSdata_Omega");
double xCM = Getd("DNSdata_x_CM");
double ecc = Getd("DNSdata_ecc");
double rdot = Getd("DNSdata_rdot");
double xmax1 = Getd("DNSdata_actual_xmax1");
double xmax2 = Getd("DNSdata_actual_xmax2");
double cxmax1 = Getd("DNSdata_xmax1");
double cxmax2 = Getd("DNSdata_xmax2");
double VolAvSigma1 = Getd("DNSdata_desired_VolAvSigma1");
double VolAvSigma2 = Getd("DNSdata_desired_VolAvSigma2");
double VolAvSigma, VolAvlSigma;
double OuterSigmaTransitionD1 = 1.0;
double OuterSigmaTransitionD2 = 1.0;

tGrid *grid = vlu->grid;
int bi;
tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vldu->vlPars;


/* do nothing if noBCs, i.e. DNSdata_Sigma_surface_BCs = none */
if(noBCs) return;


/* parse some pars: */
/* check if DNSdata_InnerToOuterSigmaTransition is only C1 or C0 */
if(Getv("DNSdata_InnerToOuterSigmaTransition","C1"))
  OuterSigmaTransitionD2 = 0.0;
if(Getv("DNSdata_InnerToOuterSigmaTransition","C0"))
  OuterSigmaTransitionD2 = OuterSigmaTransitionD1 = 0.0;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;

int n1 = box->n1;
int n2 = box->n2;
int n3 = box->n3;
int i,j,k, pln;

int isSTAR1     = (box->SIDE == STAR1);

int MATTRinside = (box->MATTR== INSIDE);

int MATTRtouch  = (box->MATTR== TOUCH);

int MATTRaway   = (box->MATTR== AWAY);

int hasSSURF    = (box->BOUND== SSURF);

int isXinDom    = (box->CI->dom == box->SIDE - STAR1);

int isVolAvBox  = (MATTRinside && hasSSURF && isXinDom);



double *FPsi = vlldataptr(vlFu, box, 0);
double *FB1 = vlldataptr(vlFu, box, 1);
double *FB2 = vlldataptr(vlFu, box, 2);
double *FB3 = vlldataptr(vlFu, box, 3);
double *FalphaP = vlldataptr(vlFu, box, 4);
double *FSigma = vlldataptr(vlFu, box, 5);
double *Psi = vlldataptr(vlu, box, 0);
double *B1 = vlldataptr(vlu, box, 1);
double *B2 = vlldataptr(vlu, box, 2);
double *B3 = vlldataptr(vlu, box, 3);
double *alphaP = vlldataptr(vlu, box, 4);
double *Sigma = vlldataptr(vlu, box, 5);
double *FlPsi = vlldataptr(vlJdu, box, 0);
double *FlB1 = vlldataptr(vlJdu, box, 1);
double *FlB2 = vlldataptr(vlJdu, box, 2);
double *FlB3 = vlldataptr(vlJdu, box, 3);
double *FlalphaP = vlldataptr(vlJdu, box, 4);
double *FlSigma = vlldataptr(vlJdu, box, 5);
double *lPsi = vlldataptr(vldu, box, 0);
double *lB1 = vlldataptr(vldu, box, 1);
double *lB2 = vlldataptr(vldu, box, 2);
double *lB3 = vlldataptr(vldu, box, 3);
double *lalphaP = vlldataptr(vldu, box, 4);
double *lSigma = vlldataptr(vldu, box, 5);
double *dlPsi1 = vlldataptr(vlduDerivs, box, 0);
double *dlPsi2 = vlldataptr(vlduDerivs, box, 1);
double *dlPsi3 = vlldataptr(vlduDerivs, box, 2);
double *ddlPsi11 = vlldataptr(vlduDerivs, box, 3);
double *ddlPsi12 = vlldataptr(vlduDerivs, box, 4);
double *ddlPsi13 = vlldataptr(vlduDerivs, box, 5);
double *ddlPsi22 = vlldataptr(vlduDerivs, box, 6);
double *ddlPsi23 = vlldataptr(vlduDerivs, box, 7);
double *ddlPsi33 = vlldataptr(vlduDerivs, box, 8);
double *dlB11 = vlldataptr(vlduDerivs, box, 9);
double *dlB12 = vlldataptr(vlduDerivs, box, 10);
double *dlB13 = vlldataptr(vlduDerivs, box, 11);
double *dlB21 = vlldataptr(vlduDerivs, box, 12);
double *dlB22 = vlldataptr(vlduDerivs, box, 13);
double *dlB23 = vlldataptr(vlduDerivs, box, 14);
double *dlB31 = vlldataptr(vlduDerivs, box, 15);
double *dlB32 = vlldataptr(vlduDerivs, box, 16);
double *dlB33 = vlldataptr(vlduDerivs, box, 17);
double *ddlB111 = vlldataptr(vlduDerivs, box, 18);
double *ddlB112 = vlldataptr(vlduDerivs, box, 19);
double *ddlB113 = vlldataptr(vlduDerivs, box, 20);
double *ddlB122 = vlldataptr(vlduDerivs, box, 21);
double *ddlB123 = vlldataptr(vlduDerivs, box, 22);
double *ddlB133 = vlldataptr(vlduDerivs, box, 23);
double *ddlB211 = vlldataptr(vlduDerivs, box, 24);
double *ddlB212 = vlldataptr(vlduDerivs, box, 25);
double *ddlB213 = vlldataptr(vlduDerivs, box, 26);
double *ddlB222 = vlldataptr(vlduDerivs, box, 27);
double *ddlB223 = vlldataptr(vlduDerivs, box, 28);
double *ddlB233 = vlldataptr(vlduDerivs, box, 29);
double *ddlB311 = vlldataptr(vlduDerivs, box, 30);
double *ddlB312 = vlldataptr(vlduDerivs, box, 31);
double *ddlB313 = vlldataptr(vlduDerivs, box, 32);
double *ddlB322 = vlldataptr(vlduDerivs, box, 33);
double *ddlB323 = vlldataptr(vlduDerivs, box, 34);
double *ddlB333 = vlldataptr(vlduDerivs, box, 35);
double *dlalphaP1 = vlldataptr(vlduDerivs, box, 36);
double *dlalphaP2 = vlldataptr(vlduDerivs, box, 37);
double *dlalphaP3 = vlldataptr(vlduDerivs, box, 38);
double *ddlalphaP11 = vlldataptr(vlduDerivs, box, 39);
double *ddlalphaP12 = vlldataptr(vlduDerivs, box, 40);
double *ddlalphaP13 = vlldataptr(vlduDerivs, box, 41);
double *ddlalphaP22 = vlldataptr(vlduDerivs, box, 42);
double *ddlalphaP23 = vlldataptr(vlduDerivs, box, 43);
double *ddlalphaP33 = vlldataptr(vlduDerivs, box, 44);
double *dlSigma1 = vlldataptr(vlduDerivs, box, 45);
double *dlSigma2 = vlldataptr(vlduDerivs, box, 46);
double *dlSigma3 = vlldataptr(vlduDerivs, box, 47);
double *ddlSigma11 = vlldataptr(vlduDerivs, box, 48);
double *ddlSigma12 = vlldataptr(vlduDerivs, box, 49);
double *ddlSigma13 = vlldataptr(vlduDerivs, box, 50);
double *ddlSigma22 = vlldataptr(vlduDerivs, box, 51);
double *ddlSigma23 = vlldataptr(vlduDerivs, box, 52);
double *ddlSigma33 = vlldataptr(vlduDerivs, box, 53);
int index_Psi = (vlu)->index[0];
int index_B1 = (vlu)->index[1];
int index_B2 = (vlu)->index[2];
int index_B3 = (vlu)->index[3];
int index_alphaP = (vlu)->index[4];
int index_Sigma = (vlu)->index[5];
int index_lPsi = (vldu)->index[0];
int index_lB1 = (vldu)->index[1];
int index_lB2 = (vldu)->index[2];
int index_lB3 = (vldu)->index[3];
int index_lalphaP = (vldu)->index[4];
int index_lSigma = (vldu)->index[5];
int index_dlPsi1 = (vlduDerivs)->index[0];
int index_dlPsi2 = (vlduDerivs)->index[1];
int index_dlPsi3 = (vlduDerivs)->index[2];
int index_ddlPsi11 = (vlduDerivs)->index[3];
int index_ddlPsi12 = (vlduDerivs)->index[4];
int index_ddlPsi13 = (vlduDerivs)->index[5];
int index_ddlPsi22 = (vlduDerivs)->index[6];
int index_ddlPsi23 = (vlduDerivs)->index[7];
int index_ddlPsi33 = (vlduDerivs)->index[8];
int index_dlB11 = (vlduDerivs)->index[9];
int index_dlB12 = (vlduDerivs)->index[10];
int index_dlB13 = (vlduDerivs)->index[11];
int index_dlB21 = (vlduDerivs)->index[12];
int index_dlB22 = (vlduDerivs)->index[13];
int index_dlB23 = (vlduDerivs)->index[14];
int index_dlB31 = (vlduDerivs)->index[15];
int index_dlB32 = (vlduDerivs)->index[16];
int index_dlB33 = (vlduDerivs)->index[17];
int index_ddlB111 = (vlduDerivs)->index[18];
int index_ddlB112 = (vlduDerivs)->index[19];
int index_ddlB113 = (vlduDerivs)->index[20];
int index_ddlB122 = (vlduDerivs)->index[21];
int index_ddlB123 = (vlduDerivs)->index[22];
int index_ddlB133 = (vlduDerivs)->index[23];
int index_ddlB211 = (vlduDerivs)->index[24];
int index_ddlB212 = (vlduDerivs)->index[25];
int index_ddlB213 = (vlduDerivs)->index[26];
int index_ddlB222 = (vlduDerivs)->index[27];
int index_ddlB223 = (vlduDerivs)->index[28];
int index_ddlB233 = (vlduDerivs)->index[29];
int index_ddlB311 = (vlduDerivs)->index[30];
int index_ddlB312 = (vlduDerivs)->index[31];
int index_ddlB313 = (vlduDerivs)->index[32];
int index_ddlB322 = (vlduDerivs)->index[33];
int index_ddlB323 = (vlduDerivs)->index[34];
int index_ddlB333 = (vlduDerivs)->index[35];
int index_dlalphaP1 = (vlduDerivs)->index[36];
int index_dlalphaP2 = (vlduDerivs)->index[37];
int index_dlalphaP3 = (vlduDerivs)->index[38];
int index_ddlalphaP11 = (vlduDerivs)->index[39];
int index_ddlalphaP12 = (vlduDerivs)->index[40];
int index_ddlalphaP13 = (vlduDerivs)->index[41];
int index_ddlalphaP22 = (vlduDerivs)->index[42];
int index_ddlalphaP23 = (vlduDerivs)->index[43];
int index_ddlalphaP33 = (vlduDerivs)->index[44];
int index_dlSigma1 = (vlduDerivs)->index[45];
int index_dlSigma2 = (vlduDerivs)->index[46];
int index_dlSigma3 = (vlduDerivs)->index[47];
int index_ddlSigma11 = (vlduDerivs)->index[48];
int index_ddlSigma12 = (vlduDerivs)->index[49];
int index_ddlSigma13 = (vlduDerivs)->index[50];
int index_ddlSigma22 = (vlduDerivs)->index[51];
int index_ddlSigma23 = (vlduDerivs)->index[52];
int index_ddlSigma33 = (vlduDerivs)->index[53];
int index_DNSdata_Sigmax = Ind("DNSdata_Sigmax");
double *dSigma1 = box->v[index_DNSdata_Sigmax + 0];
double *dSigma2 = box->v[index_DNSdata_Sigmax + 1];
double *dSigma3 = box->v[index_DNSdata_Sigmax + 2];
int index_DNSdata_Sigmaxx = Ind("DNSdata_Sigmaxx");
double *ddSigma11 = box->v[index_DNSdata_Sigmaxx + 0];
double *ddSigma12 = box->v[index_DNSdata_Sigmaxx + 1];
double *ddSigma13 = box->v[index_DNSdata_Sigmaxx + 2];
double *ddSigma22 = box->v[index_DNSdata_Sigmaxx + 3];
double *ddSigma23 = box->v[index_DNSdata_Sigmaxx + 4];
double *ddSigma33 = box->v[index_DNSdata_Sigmaxx + 5];
int index_x = Ind("x");
double *x = box->v[index_x + 0];
int index_y = Ind("y");
double *y = box->v[index_y + 0];
int index_z = Ind("z");
double *z = box->v[index_z + 0];
int index_DNSdata_q = Ind("DNSdata_q");
double *q = box->v[index_DNSdata_q + 0];
int index_DNSdata_wBx = Ind("DNSdata_wBx");
double *wB1 = box->v[index_DNSdata_wBx + 0];
double *wB2 = box->v[index_DNSdata_wBx + 1];
double *wB3 = box->v[index_DNSdata_wBx + 2];
int index_DNSdata_qx = Ind("DNSdata_qx");
double *dq1 = box->v[index_DNSdata_qx + 0];
double *dq2 = box->v[index_DNSdata_qx + 1];
double *dq3 = box->v[index_DNSdata_qx + 2];
int index_dXdx = Ind("dXdx");
double *dlam1 = box->v[index_dXdx + 0];
double *dlam2 = box->v[index_dXdx + 1];
double *dlam3 = box->v[index_dXdx + 2];
int index_DNSdata_SigmaX = Ind("DNSdata_SigmaX");
double *dSigmadlam = box->v[index_DNSdata_SigmaX + 0];
int index_DNSdata_SigmaXX = Ind("DNSdata_SigmaXX");
double *ddSigmadlam2 = box->v[index_DNSdata_SigmaXX + 0];
int index_DNSdata_SigmaXXX = Ind("DNSdata_SigmaXXX");
double *dddSigmadlam3 = box->v[index_DNSdata_SigmaXXX + 0];
int index_DNSdata_lSigmaX = Ind("DNSdata_lSigmaX");
double *dlSigmadlam = box->v[index_DNSdata_lSigmaX + 0];
int index_DNSdata_lSigmaXX = Ind("DNSdata_lSigmaXX");
double *ddlSigmadlam2 = box->v[index_DNSdata_lSigmaXX + 0];
int index_DNSdata_lSigmaXXX = Ind("DNSdata_lSigmaXXX");
double *dddlSigmadlam3 = box->v[index_DNSdata_lSigmaXXX + 0];


/* Jetzt geht's los! */
double alpha;
double alpha2;
double beta1;
double beta2;
double beta3;
double ddlSig11;
double ddlSig12;
double ddlSig13;
double ddlSig22;
double ddlSig23;
double ddlSig33;
double ddlSigin11;
double ddlSigin12;
double ddlSigin13;
double ddlSigin22;
double ddlSigin23;
double ddlSigin33;
double ddSig11;
double ddSig12;
double ddSig13;
double ddSig22;
double ddSig23;
double ddSig33;
double ddSigin11;
double ddSigin12;
double ddSigin13;
double ddSigin22;
double ddSigin23;
double ddSigin33;
double dlQ1;
double dlQ2;
double dlQ3;
double dlSig1;
double dlSig2;
double dlSig3;
double dlSigin1;
double dlSigin2;
double dlSigin3;
double dlSigmaUp1;
double dlSigmaUp2;
double dlSigmaUp3;
double dlwB11;
double dlwB12;
double dlwB13;
double dlwB21;
double dlwB22;
double dlwB23;
double dlwB31;
double dlwB32;
double dlwB33;
double dQ1;
double dQ2;
double dQ3;
double dSig1;
double dSig2;
double dSig3;
double dSigin1;
double dSigin2;
double dSigin3;
double dSigmaUp1;
double DSigmaUp1;
double dSigmaUp2;
double DSigmaUp2;
double dSigmaUp3;
double DSigmaUp3;
double h;
double h2;
double L2;
double lalpha;
double lh;
double lhuzeroPsi4beta1;
double lhuzeroPsi4beta2;
double lhuzeroPsi4beta3;
double lL2;
double lLnh;
double lq;
double luzero;
double luzerosqr;
double lwB1;
double lwB2;
double lwB3;
double nv1;
double nv2;
double nv3;
double nvm;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double Psi2;
double Psi3;
double Psi4;
double Psim2;
double Psim3;
double Psim4;
double Psim5;
double Psim6;
double Psim7;
double Psim8;
double Psim9;
double rdotor;
double uzero;
double uzerosqr;
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
double xrdotor1;
double xrdotor2;
double xrdotor3;





/* conditional */
if (!hasSSURF) {


continue; 

}
/* if (!hasSSURF) */



if(blkinfo != NULL) {                                               
                     if(blkinfo->bi   != bi) continue;
                     if(blkinfo->vari != index_lSigma) continue;  }

/* conditional */
if ((isSTAR1 && corot1) || (!isSTAR1 && corot2)) {



/* conditional */
if (nonlin) {


forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k); 

FSigma[ijk]
=
Sigma[ijk]
;


} /* end forplane1 */ 


forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k); 

FSigma[ijk]
=
Sigma[ijk]
;


} /* end forplane1 */ 


} else { /* if (!nonlin) */


forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k); 

FlSigma[ijk]
=
lSigma[ijk]
;


} /* end forplane1 */ 


forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k); 

FlSigma[ijk]
=
lSigma[ijk]
;


} /* end forplane1 */ 

}
/* if (nonlin) */



continue; /* for corot we are done with this box */ 

}
/* if (nonlin) */




/* conditional */
if (MATTRtouch) {


int    ind0, indin, biin, n1in,n2in,n3in; 


double Sig, Sigin, lSig, lSigin; 


double *Sigmain; 


double *lSigmain; 


double *dSigmain1, *dSigmain2, *dSigmain3; 


double *dlSigmain1, *dlSigmain2, *dlSigmain3; 


double *ddSigmain11, *ddSigmain12, *ddSigmain13,                      
                            *ddSigmain22, *ddSigmain23, *ddSigmain33;

double *ddlSigmain11, *ddlSigmain12, *ddlSigmain13,                      
                            *ddlSigmain22, *ddlSigmain23, *ddlSigmain33;

biin = bi - 6; /* only if CubSp are arranged my way */ 


n1in = grid->box[biin]->n1;                                              
                     n2in = grid->box[biin]->n2;
                     n3in = grid->box[biin]->n3;

               Sigmain    = grid->box[biin]->v[index_Sigma];
               lSigmain   = grid->box[biin]->v[index_lSigma];
               dSigmain1  = grid->box[biin]->v[index_DNSdata_Sigmax+0];
               dSigmain2  = grid->box[biin]->v[index_DNSdata_Sigmax+1];
               dSigmain3  = grid->box[biin]->v[index_DNSdata_Sigmax+2];
              dlSigmain1  = grid->box[biin]->v[index_dlSigma1];
              dlSigmain2  = grid->box[biin]->v[index_dlSigma2];
              dlSigmain3  = grid->box[biin]->v[index_dlSigma3];
              ddSigmain11 = grid->box[biin]->v[index_DNSdata_Sigmaxx+0];
              ddSigmain12 = grid->box[biin]->v[index_DNSdata_Sigmaxx+1];
              ddSigmain13 = grid->box[biin]->v[index_DNSdata_Sigmaxx+2];
              ddSigmain22 = grid->box[biin]->v[index_DNSdata_Sigmaxx+3];
              ddSigmain23 = grid->box[biin]->v[index_DNSdata_Sigmaxx+4];
              ddSigmain33 = grid->box[biin]->v[index_DNSdata_Sigmaxx+5];
             ddlSigmain11 = grid->box[biin]->v[index_ddlSigma11];
             ddlSigmain12 = grid->box[biin]->v[index_ddlSigma12];
             ddlSigmain13 = grid->box[biin]->v[index_ddlSigma13];
             ddlSigmain22 = grid->box[biin]->v[index_ddlSigma22];
             ddlSigmain23 = grid->box[biin]->v[index_ddlSigma23];
             ddlSigmain33 = grid->box[biin]->v[index_ddlSigma33];
if(n2in!=n2 || n3in!=n3) errorexit("we need n2in=n2 and n3in=n3"); 



/* conditional */
if (nonlin) {


spec_Deriv2(box, 1, Sigma, ddSigmadlam2); 


spec_Deriv1(box, 1, ddSigmadlam2, dddSigmadlam3); 


spec_Deriv1(box, 1, Sigma, dSigmadlam); 


forplane1(i,j,k, n1,n2,n3, 1){ ijk=Index(i,j,k); 


ind0   = Ind_n1n2(0,j,k,n1,n2);                                 
                       indin  = Ind_n1n2(n1in-1,j,k,n1in,n2in);
                       Sig    = Sigma[ind0];
                       Sigin  = Sigmain[indin];FSigma[ijk]
=
Sig - Sigin
;


} /* endfor */ 



/* conditional */
if (isSTAR1) {

xc
=
cxmax1
;


} else { /* if (!isSTAR1) */

xc
=
cxmax2
;

}
/* if (isSTAR1) */



forplane1(i,j,k, n1,n2,n3, 0){ 


ijk=Index(0,j,k); /* set index to i=0 */ 

nv1
=
-xc + x[ijk]
;

nv2
=
y[ijk]
;

nv3
=
z[ijk]
;

nvm
=
Sqrt(pow2(nv1) + pow2(nv2) + pow2(nv3))
;

nv1
=
nv1/nvm
;

nv2
=
nv2/nvm
;

nv3
=
nv3/nvm
;

dSig1
=
dSigma1[ijk]
;

dSig2
=
dSigma2[ijk]
;

dSig3
=
dSigma3[ijk]
;


ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */ 

dSigin1
=
dSigmain1[ijk]
;

dSigin2
=
dSigmain2[ijk]
;

dSigin3
=
dSigmain3[ijk]
;


ijk=Index(0,j,k); /* set index to i=0 */ 

FSigma[ijk]
=
nv1*(dSig1 - dSigin1*OuterSigmaTransitionD1) + 
  nv2*(dSig2 - dSigin2*OuterSigmaTransitionD1) + 
  nv3*(dSig3 - dSigin3*OuterSigmaTransitionD1)
;


} /* endfor */ 



/* conditional */
if (isSTAR1) {

xc
=
cxmax1
;


} else { /* if (!isSTAR1) */

xc
=
cxmax2
;

}
/* if (isSTAR1) */



forplane1(i,j,k, n1,n2,n3, 2){ 


ijk=Index(0,j,k); /* set index to i=0 */ 

nv1
=
-xc + x[ijk]
;

nv2
=
y[ijk]
;

nv3
=
z[ijk]
;

nvm
=
Sqrt(pow2(nv1) + pow2(nv2) + pow2(nv3))
;

nv1
=
nv1/nvm
;

nv2
=
nv2/nvm
;

nv3
=
nv3/nvm
;

ddSig11
=
ddSigma11[ijk]
;

ddSig12
=
ddSigma12[ijk]
;

ddSig13
=
ddSigma13[ijk]
;

ddSig22
=
ddSigma22[ijk]
;

ddSig23
=
ddSigma23[ijk]
;

ddSig33
=
ddSigma33[ijk]
;


ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */ 

ddSigin11
=
ddSigmain11[ijk]
;

ddSigin12
=
ddSigmain12[ijk]
;

ddSigin13
=
ddSigmain13[ijk]
;

ddSigin22
=
ddSigmain22[ijk]
;

ddSigin23
=
ddSigmain23[ijk]
;

ddSigin33
=
ddSigmain33[ijk]
;


ijk=Index(2,j,k); /* set index to i=2 */ 

FSigma[ijk]
=
2.*(ddSig23*nv2*nv3 + nv1*(ddSig12*nv2 + ddSig13*nv3) - 
     (ddSigin23*nv2*nv3 + nv1*(ddSigin12*nv2 + ddSigin13*nv3))*
      OuterSigmaTransitionD2) + 
  (ddSig11 - ddSigin11*OuterSigmaTransitionD2)*pow2(nv1) + 
  (ddSig22 - ddSigin22*OuterSigmaTransitionD2)*pow2(nv2) + 
  (ddSig33 - ddSigin33*OuterSigmaTransitionD2)*pow2(nv3)
;


} /* endfor */ 


} else { /* if (!isSTAR1) */


spec_Deriv2(box, 1, lSigma, ddlSigmadlam2); 


spec_Deriv1(box, 1, ddlSigmadlam2, dddlSigmadlam3); 


spec_Deriv1(box, 1, lSigma, dlSigmadlam); 


forplane1(i,j,k, n1,n2,n3, 1){ ijk=Index(i,j,k); 


ind0   = Ind_n1n2(0,j,k,n1,n2);                                 
                       indin  = Ind_n1n2(n1in-1,j,k,n1in,n2in);
                       lSig   = lSigma[ind0];
                       lSigin = lSigmain[indin];FlSigma[ijk]
=
lSig - lSigin
;


} /* endfor */ 



/* conditional */
if (isSTAR1) {

xc
=
cxmax1
;


} else { /* if (!isSTAR1) */

xc
=
cxmax2
;

}
/* if (isSTAR1) */



forplane1(i,j,k, n1,n2,n3, 0){ 


ijk=Index(0,j,k); /* set index to i=0 */ 

nv1
=
-xc + x[ijk]
;

nv2
=
y[ijk]
;

nv3
=
z[ijk]
;

nvm
=
Sqrt(pow2(nv1) + pow2(nv2) + pow2(nv3))
;

nv1
=
nv1/nvm
;

nv2
=
nv2/nvm
;

nv3
=
nv3/nvm
;

dlSig1
=
dlSigma1[ijk]
;

dlSig2
=
dlSigma2[ijk]
;

dlSig3
=
dlSigma3[ijk]
;


ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */ 

dlSigin1
=
dlSigmain1[ijk]
;

dlSigin2
=
dlSigmain2[ijk]
;

dlSigin3
=
dlSigmain3[ijk]
;


ijk=Index(0,j,k); /* set index to i=0 */ 

FlSigma[ijk]
=
nv1*(dlSig1 - dlSigin1*OuterSigmaTransitionD1) + 
  nv2*(dlSig2 - dlSigin2*OuterSigmaTransitionD1) + 
  nv3*(dlSig3 - dlSigin3*OuterSigmaTransitionD1)
;


} /* endfor */ 



/* conditional */
if (isSTAR1) {

xc
=
cxmax1
;


} else { /* if (!isSTAR1) */

xc
=
cxmax2
;

}
/* if (isSTAR1) */



forplane1(i,j,k, n1,n2,n3, 2){ 


ijk=Index(0,j,k); /* set index to i=0 */ 

nv1
=
-xc + x[ijk]
;

nv2
=
y[ijk]
;

nv3
=
z[ijk]
;

nvm
=
Sqrt(pow2(nv1) + pow2(nv2) + pow2(nv3))
;

nv1
=
nv1/nvm
;

nv2
=
nv2/nvm
;

nv3
=
nv3/nvm
;

ddlSig11
=
ddlSigma11[ijk]
;

ddlSig12
=
ddlSigma12[ijk]
;

ddlSig13
=
ddlSigma13[ijk]
;

ddlSig22
=
ddlSigma22[ijk]
;

ddlSig23
=
ddlSigma23[ijk]
;

ddlSig33
=
ddlSigma33[ijk]
;


ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */ 

ddlSigin11
=
ddlSigmain11[ijk]
;

ddlSigin12
=
ddlSigmain12[ijk]
;

ddlSigin13
=
ddlSigmain13[ijk]
;

ddlSigin22
=
ddlSigmain22[ijk]
;

ddlSigin23
=
ddlSigmain23[ijk]
;

ddlSigin33
=
ddlSigmain33[ijk]
;


ijk=Index(2,j,k); /* set index to i=2 */ 

FlSigma[ijk]
=
2.*(ddlSig23*nv2*nv3 + nv1*(ddlSig12*nv2 + ddlSig13*nv3) - 
     (ddlSigin23*nv2*nv3 + nv1*(ddlSigin12*nv2 + ddlSigin13*nv3))*
      OuterSigmaTransitionD2) + 
  (ddlSig11 - ddlSigin11*OuterSigmaTransitionD2)*pow2(nv1) + 
  (ddlSig22 - ddlSigin22*OuterSigmaTransitionD2)*pow2(nv2) + 
  (ddlSig33 - ddlSigin33*OuterSigmaTransitionD2)*pow2(nv3)
;


} /* endfor */ 

}
/* if (isSTAR1) */


}
/* if (isSTAR1) */




/* conditional */
if (MATTRinside) {



/* conditional */
if (dqFromqg) {


FirstDerivsOf_S(box, Ind("DNSdata_qg"), 					Ind("DNSdata_qx")); 


} else { /* if (!dqFromqg) */


FirstDerivsOf_S(box, Ind("DNSdata_q"), 					Ind("DNSdata_qx")); 

}
/* if (dqFromqg) */




/* conditional */
if (nonlin) {


FirstDerivsOf_S(box, index_Sigma,                                        Ind("DNSdata_Sigmax")); 



/* conditional */
if (AddInnerVolIntToBC || InnerVolIntZero) {


VolAvSigma = 0.0; 


if(isVolAvBox) VolAvSigma =                                             
                         BoxVolumeIntegral(grid->box[bi], index_Sigma);
}
/* if (AddInnerVolIntToBC || InnerVolIntZero) */




/* conditional */
if (AddInnerSumToBC || InnerSumZero) {


VolAvSigma = 0.0; 


if(isVolAvBox) forallpoints(box, ijk) { 


VolAvSigma += Sigma[ijk]; 


} /* endfor */ 

}
/* if (AddInnerSumToBC || InnerSumZero) */




/* conditional */
if (isVolAvBox) {


VolAvSigma = VolAvSigma - VolAvSigma1; 


} else { /* if (!isVolAvBox) */


VolAvSigma = VolAvSigma - VolAvSigma2; 

}
/* if (isVolAvBox) */



//printf("VolAvSigma=%g\n",VolAvSigma); 



/* conditional */
if (isSTAR1) {

xmax
=
xmax1
;


} else { /* if (!isSTAR1) */

xmax
=
xmax2
;

}
/* if (isSTAR1) */



forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k); 

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

Psi2
=
pow2(Psi[ijk])
;

Psi4
=
pow2(Psi2)
;

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

alpha
=
alphaP[ijk]/Psi[ijk]
;

alpha2
=
pow2(alpha)
;

h
=
1. + q[ijk]
;

h2
=
pow2(h)
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

dSigmaUp1
=
dSigma1[ijk]
;

dSigmaUp2
=
dSigma2[ijk]
;

dSigmaUp3
=
dSigma3[ijk]
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

L2
=
h2 + (DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) + 
  (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) + 
  (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk])
;

uzerosqr
=
L2/(alpha2*h2)
;

uzero
=
sqrt(uzerosqr)
;



/* conditional */
if (dQFromdlam) {

dQ1
=
dlam1[ijk]
;

dQ2
=
dlam2[ijk]
;

dQ3
=
dlam3[ijk]
;


} else { /* if (!dQFromdlam) */

dQ1
=
dq1[ijk]
;

dQ2
=
dq2[ijk]
;

dQ3
=
dq3[ijk]
;

}
/* if (dQFromdlam) */




/* conditional */
if (ImposeActualBC) {

FSigma[ijk]
=
dQ1*(dSigmaUp1 - beta1*h*Psi4*uzero) + 
  dQ2*(dSigmaUp2 - beta2*h*Psi4*uzero) + dQ3*(dSigmaUp3 - beta3*h*Psi4*uzero)
;

FSigma[ijk]
=
FSigma[ijk] + Psim2*(dQ1*wB1[ijk] + dQ2*wB2[ijk] + dQ3*wB3[ijk])
;

}
/* if (ImposeActualBC) */




/* conditional */
if (AddInnerVolIntToBC || AddInnerSumToBC) {

FSigma[ijk]
=
VolAvSigma + FSigma[ijk]
;

}
/* if (AddInnerVolIntToBC || AddInnerSumToBC) */



} /* end forplane1 */ 


if(isVolAvBox) { 



/* conditional */
if (SigmaZeroAtPoint) {


ijk=Index(n1-1,0,0); 

FSigma[ijk]
=
Sigma[ijk]
;

}
/* if (SigmaZeroAtPoint) */




/* conditional */
if (InnerVolIntZero || InnerSumZero) {


ijk=Index(n1-1,0,0); 

FSigma[ijk]
=
VolAvSigma
;

}
/* if (InnerVolIntZero || InnerSumZero) */



} /* end if(isVolAvBox) */ 


} else { /* if (!InnerVolIntZero || InnerSumZero) */


FirstDerivsOf_S(box, index_lSigma, index_dlSigma1); 



/* conditional */
if (AddInnerVolIntToBC || InnerVolIntZero) {


VolAvlSigma = 0.0; 


if(isVolAvBox) VolAvlSigma =                              
          BoxVolumeIntegral(grid->box[bi], index_lSigma);
}
/* if (AddInnerVolIntToBC || InnerVolIntZero) */




/* conditional */
if (AddInnerSumToBC || InnerSumZero) {


VolAvlSigma = 0.0; 


if(isVolAvBox) forallpoints(box, ijk) { 


VolAvlSigma += lSigma[ijk]; 


} /* endfor */ 

}
/* if (AddInnerSumToBC || InnerSumZero) */



//if(VolAvlSigma!=0.0) printf("box->b=%d VolAvlSigma=%g\n",box->b,VolAvlSigma); 



/* conditional */
if (isSTAR1) {

xmax
=
xmax1
;


} else { /* if (!isSTAR1) */

xmax
=
xmax2
;

}
/* if (isSTAR1) */



forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k); 

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

Psi2
=
pow2(Psi[ijk])
;

Psi3
=
Psi2*Psi[ijk]
;

Psi4
=
pow2(Psi2)
;

Psim2
=
1/Psi2
;

Psim3
=
1/Psi3
;

Psim4
=
pow2(Psim2)
;

Psim6
=
Psim2*Psim4
;

Psim8
=
pow2(Psim4)
;

Psim5
=
Psim6*Psi[ijk]
;

Psim7
=
Psim2*Psim5
;

Psim9
=
Psim2*Psim7
;

alpha
=
alphaP[ijk]/Psi[ijk]
;

alpha2
=
pow2(alpha)
;

h
=
1. + q[ijk]
;

h2
=
pow2(h)
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

dSigmaUp1
=
dSigma1[ijk]
;

dSigmaUp2
=
dSigma2[ijk]
;

dSigmaUp3
=
dSigma3[ijk]
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

L2
=
h2 + (DSigmaUp1 + w1)*(wDown1 + dSigma1[ijk]) + 
  (DSigmaUp2 + w2)*(wDown2 + dSigma2[ijk]) + 
  (DSigmaUp3 + w3)*(wDown3 + dSigma3[ijk])
;

uzerosqr
=
L2/(alpha2*h2)
;

uzero
=
sqrt(uzerosqr)
;

lq
=
0
;

dlQ1
=
0
;

dlQ2
=
0
;

dlQ3
=
0
;

lh
=
0
;

lLnh
=
0
;

lwB1
=
0
;

lwB2
=
0
;

lwB3
=
0
;

dlwB11
=
0
;

dlwB12
=
0
;

dlwB13
=
0
;

dlwB21
=
0
;

dlwB22
=
0
;

dlwB23
=
0
;

dlwB31
=
0
;

dlwB32
=
0
;

dlwB33
=
0
;

lalpha
=
-((alphaP[ijk]*lPsi[ijk])/Psi2) + lalphaP[ijk]/Psi[ijk]
;

dlSigmaUp1
=
dlSigma1[ijk]
;

dlSigmaUp2
=
dlSigma2[ijk]
;

dlSigmaUp3
=
dlSigma3[ijk]
;

lL2
=
4.*h2*lLnh - lPsi[ijk]*(8.*Psim5*
      (dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
        dSigmaUp3*dSigma3[ijk]) + 
     (16.*Psim9*wBDown1 + 24.*Psim7*dSigma1[ijk])*wB1[ijk] + 
     (16.*Psim9*wBDown2 + 24.*Psim7*dSigma2[ijk])*wB2[ijk] + 
     (16.*Psim9*wBDown3 + 24.*Psim7*dSigma3[ijk])*wB3[ijk]) + 
  2.*(Psim8*(lwB1*wBDown1 + lwB2*wBDown2 + lwB3*wBDown3) + 
     Psim4*(dSigmaUp1*dlSigma1[ijk] + dSigmaUp2*dlSigma2[ijk] + 
        dSigmaUp3*dlSigma3[ijk]) + 
     Psim6*(lwB1*dSigma1[ijk] + lwB2*dSigma2[ijk] + lwB3*dSigma3[ijk] + 
        dlSigma1[ijk]*wB1[ijk] + dlSigma2[ijk]*wB2[ijk] + 
        dlSigma3[ijk]*wB3[ijk]))
;

luzerosqr
=
(lL2 - 2.*L2*(lalpha/alpha + lLnh))/(alpha2*h2)
;

luzero
=
(0.5*luzerosqr)/uzero
;

lhuzeroPsi4beta1
=
h*Psi4*uzero*lB1[ijk] + beta1*(Psi4*(h*luzero + lh*uzero) + 
     4.*h*Psi3*uzero*lPsi[ijk])
;

lhuzeroPsi4beta2
=
h*Psi4*uzero*lB2[ijk] + beta2*(Psi4*(h*luzero + lh*uzero) + 
     4.*h*Psi3*uzero*lPsi[ijk])
;

lhuzeroPsi4beta3
=
h*Psi4*uzero*lB3[ijk] + beta3*(Psi4*(h*luzero + lh*uzero) + 
     4.*h*Psi3*uzero*lPsi[ijk])
;



/* conditional */
if (dQFromdlam) {

dQ1
=
dlam1[ijk]
;

dQ2
=
dlam2[ijk]
;

dQ3
=
dlam3[ijk]
;


} else { /* if (!dQFromdlam) */

dQ1
=
dq1[ijk]
;

dQ2
=
dq2[ijk]
;

dQ3
=
dq3[ijk]
;

}
/* if (dQFromdlam) */




/* conditional */
if (ImposeActualBC) {

FlSigma[ijk]
=
dQ1*(dlSigmaUp1 - lhuzeroPsi4beta1) + dQ2*(dlSigmaUp2 - lhuzeroPsi4beta2) + 
  dQ3*(dlSigmaUp3 - lhuzeroPsi4beta3) + 
  dlQ1*(dSigmaUp1 - beta1*h*Psi4*uzero) + 
  dlQ2*(dSigmaUp2 - beta2*h*Psi4*uzero) + 
  dlQ3*(dSigmaUp3 - beta3*h*Psi4*uzero)
;

FlSigma[ijk]
=
FlSigma[ijk] + Psim2*(dQ1*lwB1 + dQ2*lwB2 + dQ3*lwB3 + dlQ1*wB1[ijk] + 
     dlQ2*wB2[ijk] + dlQ3*wB3[ijk]) - 
  2.*Psim3*lPsi[ijk]*(dQ1*wB1[ijk] + dQ2*wB2[ijk] + dQ3*wB3[ijk])
;

}
/* if (ImposeActualBC) */




/* conditional */
if (AddInnerVolIntToBC || AddInnerSumToBC) {

FlSigma[ijk]
=
VolAvlSigma + FlSigma[ijk]
;

}
/* if (AddInnerVolIntToBC || AddInnerSumToBC) */



} /* end forplane1 */ 


if(isVolAvBox) { 



/* conditional */
if (SigmaZeroAtPoint) {


ijk=Index(n1-1,0,0); 

FlSigma[ijk]
=
lSigma[ijk]
;

}
/* if (SigmaZeroAtPoint) */




/* conditional */
if (InnerVolIntZero || InnerSumZero) {


ijk=Index(n1-1,0,0); 

FlSigma[ijk]
=
VolAvlSigma
;

}
/* if (InnerVolIntZero || InnerSumZero) */



} /* end if(isVolAvBox) */ 

}
/* if (InnerVolIntZero || InnerSumZero) */


}
/* if (InnerVolIntZero || InnerSumZero) */




/* conditional */
if (MATTRinside && KeepInnerSigma) {


int bb;   int star=grid->box[bi]->SIDE; 


forallboxes(grid,bb) { 


tBox *bo = grid->box[bb]; 


if(bo->SIDE != star) continue; 


if(bo->MATTR != INSIDE) continue; 



/* conditional */
if (nonlin) {


double *FSigma_bb = vlldataptr(vlFu, bo, 5); 


forallpoints(bo, ijk)                            
                           FSigma_bb[ijk] = 0.0;

} else { /* if (!nonlin) */


double *FlSigma_bb = vlldataptr(vlJdu, bo, 5); 


double *lSigma_bb  = vlldataptr( vldu, bo, 5); 


forallpoints(bo, ijk)                                        
                           FlSigma_bb[ijk] = lSigma_bb[ijk];
}
/* if (nonlin) */



} /* endfor bb */ 

}
/* if (nonlin) */




/* conditional */
if ((MATTRaway || MATTRtouch) && SigmaZeroInOuterBoxes) {



/* conditional */
if (nonlin) {


forallpoints(box, ijk) { 

FSigma[ijk]
=
Sigma[ijk]
;


} /* endfor */ 


} else { /* if (!nonlin) */


forallpoints(box, ijk) { 

FlSigma[ijk]
=
lSigma[ijk]
;


} /* endfor */ 

}
/* if (nonlin) */


}
/* if (nonlin) */



/* end all */ 

} /* end of boxes */


}  /* end of function */

/* set_DNSdata_Sigma_BCs.c */
/* nvars = 124, n* = 602,  n/ = 279,  n+ = 387, n = 1268, O = 1 */
