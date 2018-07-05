/* DNS_CTS.c */
/* Copyright (C) 2005-2008 Wolfgang Tichy, 5.7.2018 */
/* Produced with Mathematica */

#include "sgrid.h"
#include "DNSdata.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Abs(x)     (fabs((double) (x)))
#define Log(x)     (log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))




void DNS_CTS(tVarList *vlFu, tVarList *vlu,       tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs,      int nonlin)
{
int iStart = Set_pdb_iStart_AtPar("DNSdata_rotationstate1");
int VwApprox1 = Getv("DNSdata_rotationstate1","VwApproximation");
int VwApprox2 = Getv("DNSdata_rotationstate2","VwApproximation");
int corot1 = VwApprox1 || Getv("DNSdata_rotationstate1","corotation");
int corot2 = VwApprox2 || Getv("DNSdata_rotationstate2","corotation");
int VwApprox, corot;
int dqFromqg = Getv("DNSdata_q_derivs","dqg");
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
int FakeMatterOutside = Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside");
int LaplaceSigmaOutside = Getv("DNSdata_Sigma_surface_BCs","LaplaceSigmaOutside");
int FakeT0 = Getv("DNSdata_FakeMatterType","rhoEQ-lam");

tGrid *grid = vlu->grid;
int bi;
tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vldu->vlPars;

forallboxes(grid,bi)
{
tBox *box = grid->box[bi];
int ijk;

int MATTRinside = (box->MATTR== INSIDE);

int MATTRtouch  = (box->MATTR== TOUCH);

int isSTAR1     = (box->SIDE == STAR1);



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
int index_DNSdata_Psix = Ind("DNSdata_Psix");
double *dPsi1 = box->v[index_DNSdata_Psix + 0];
double *dPsi2 = box->v[index_DNSdata_Psix + 1];
double *dPsi3 = box->v[index_DNSdata_Psix + 2];
int index_DNSdata_Psixx = Ind("DNSdata_Psixx");
double *ddPsi11 = box->v[index_DNSdata_Psixx + 0];
double *ddPsi12 = box->v[index_DNSdata_Psixx + 1];
double *ddPsi13 = box->v[index_DNSdata_Psixx + 2];
double *ddPsi22 = box->v[index_DNSdata_Psixx + 3];
double *ddPsi23 = box->v[index_DNSdata_Psixx + 4];
double *ddPsi33 = box->v[index_DNSdata_Psixx + 5];
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
int index_DNSdata_Bxxx = Ind("DNSdata_Bxxx");
double *ddB111 = box->v[index_DNSdata_Bxxx + 0];
double *ddB112 = box->v[index_DNSdata_Bxxx + 1];
double *ddB113 = box->v[index_DNSdata_Bxxx + 2];
double *ddB122 = box->v[index_DNSdata_Bxxx + 3];
double *ddB123 = box->v[index_DNSdata_Bxxx + 4];
double *ddB133 = box->v[index_DNSdata_Bxxx + 5];
double *ddB211 = box->v[index_DNSdata_Bxxx + 6];
double *ddB212 = box->v[index_DNSdata_Bxxx + 7];
double *ddB213 = box->v[index_DNSdata_Bxxx + 8];
double *ddB222 = box->v[index_DNSdata_Bxxx + 9];
double *ddB223 = box->v[index_DNSdata_Bxxx + 10];
double *ddB233 = box->v[index_DNSdata_Bxxx + 11];
double *ddB311 = box->v[index_DNSdata_Bxxx + 12];
double *ddB312 = box->v[index_DNSdata_Bxxx + 13];
double *ddB313 = box->v[index_DNSdata_Bxxx + 14];
double *ddB322 = box->v[index_DNSdata_Bxxx + 15];
double *ddB323 = box->v[index_DNSdata_Bxxx + 16];
double *ddB333 = box->v[index_DNSdata_Bxxx + 17];
int index_DNSdata_alphaPx = Ind("DNSdata_alphaPx");
double *dalphaP1 = box->v[index_DNSdata_alphaPx + 0];
double *dalphaP2 = box->v[index_DNSdata_alphaPx + 1];
double *dalphaP3 = box->v[index_DNSdata_alphaPx + 2];
int index_DNSdata_alphaPxx = Ind("DNSdata_alphaPxx");
double *ddalphaP11 = box->v[index_DNSdata_alphaPxx + 0];
double *ddalphaP12 = box->v[index_DNSdata_alphaPxx + 1];
double *ddalphaP13 = box->v[index_DNSdata_alphaPxx + 2];
double *ddalphaP22 = box->v[index_DNSdata_alphaPxx + 3];
double *ddalphaP23 = box->v[index_DNSdata_alphaPxx + 4];
double *ddalphaP33 = box->v[index_DNSdata_alphaPxx + 5];
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
int index_alpha = Ind("alpha");
double *alpha = box->v[index_alpha + 0];
int index_betax = Ind("betax");
double *beta1 = box->v[index_betax + 0];
double *beta2 = box->v[index_betax + 1];
double *beta3 = box->v[index_betax + 2];
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
int index_DNSdata_wBxx = Ind("DNSdata_wBxx");
double *dwB11 = box->v[index_DNSdata_wBxx + 0];
double *dwB12 = box->v[index_DNSdata_wBxx + 1];
double *dwB13 = box->v[index_DNSdata_wBxx + 2];
double *dwB21 = box->v[index_DNSdata_wBxx + 3];
double *dwB22 = box->v[index_DNSdata_wBxx + 4];
double *dwB23 = box->v[index_DNSdata_wBxx + 5];
double *dwB31 = box->v[index_DNSdata_wBxx + 6];
double *dwB32 = box->v[index_DNSdata_wBxx + 7];
double *dwB33 = box->v[index_DNSdata_wBxx + 8];
int index_DNSdata_VRx = Ind("DNSdata_VRx");
double *VR1 = box->v[index_DNSdata_VRx + 0];
double *VR2 = box->v[index_DNSdata_VRx + 1];
double *VR3 = box->v[index_DNSdata_VRx + 2];
int index_X = Ind("X");
double *lam = box->v[index_X + 0];
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
int index_DNSdata_CoordFac = Ind("DNSdata_CoordFac");
double *CoordFac = box->v[index_DNSdata_CoordFac + 0];


double xmax, xC;
double omegMOmeg1, omegMOmeg2, omegMOmeg3;
if(blkinfo!=NULL) if(bi!=blkinfo->bi) continue;
if(nonlin) {
D_and_DD_of_S(box, index_Psi,          Ind("DNSdata_Psix"), Ind("DNSdata_Psixx"));
D_and_DD_of_Sa(box, index_B1,          Ind("DNSdata_Bxx"), Ind("DNSdata_Bxxx"));
D_and_DD_of_S(box, index_alphaP,         Ind("DNSdata_alphaPx"), Ind("DNSdata_alphaPxx"));
D_and_DD_of_S(box, index_Sigma,         Ind("DNSdata_Sigmax"), Ind("DNSdata_Sigmaxx"));
spec_Deriv2(box, 1, Sigma, ddSigmadlam2);
spec_Deriv1(box, 1, ddSigmadlam2, dddSigmadlam3);
spec_Deriv1(box, 1, Sigma, dSigmadlam);
} else {
if(blkinfo!=NULL) {
  if(blkinfo->vari == index_lPsi)                          D_and_DD_of_S(box, index_lPsi,             index_dlPsi1, index_ddlPsi11);
  if(blkinfo->vari == index_lB1)                          D_and_DD_of_S(box, index_lB1,            index_dlB11, index_ddlB111);
  if(blkinfo->vari == index_lB2)                          D_and_DD_of_S(box, index_lB2,            index_dlB21, index_ddlB211);
  if(blkinfo->vari == index_lB3)                          D_and_DD_of_S(box, index_lB3,            index_dlB31, index_ddlB311);
  if(blkinfo->vari == index_lalphaP)                          D_and_DD_of_S(box, index_lalphaP,            index_dlalphaP1, index_ddlalphaP11);
  if(blkinfo->vari == index_lSigma) {
    D_and_DD_of_S(box, index_lSigma,      index_dlSigma1, index_ddlSigma11);
    spec_Deriv2(box, 1, lSigma, ddlSigmadlam2);
    spec_Deriv1(box, 1, ddlSigmadlam2, dddlSigmadlam3);
    spec_Deriv1(box, 1, lSigma, dlSigmadlam);
  } /* end blkinfo->vari == index_lSigma */
} /* end blkinfo!=NULL  */
else {
  D_and_DD_of_S(box, index_lPsi,     index_dlPsi1, index_ddlPsi11);
  D_and_DD_of_Sa(box, index_lB1,     index_dlB11, index_ddlB111);
  D_and_DD_of_S(box, index_lalphaP,    index_dlalphaP1, index_ddlalphaP11);
  D_and_DD_of_S(box, index_lSigma,    index_dlSigma1, index_ddlSigma11);
  spec_Deriv2(box, 1, lSigma, ddlSigmadlam2);
  spec_Deriv1(box, 1, ddlSigmadlam2, dddlSigmadlam3);
  spec_Deriv1(box, 1, lSigma, dlSigmadlam);
} /* end else */
} /* end linear case */
FirstDerivsOf_Sa(box, Ind("DNSdata_wBx"),       Ind("DNSdata_wBxx"));
if(dqFromqg) {
FirstDerivsOf_S(box,  Ind("DNSdata_qg"),          Ind("DNSdata_qx"));
} else {
FirstDerivsOf_S(box,  Ind("DNSdata_q"),          Ind("DNSdata_qx"));
} /* end if */
VwApprox = corot = 0;
if(isSTAR1) {
if(corot1)    corot = 1;
if(VwApprox1) VwApprox = 1;
xmax = xmax1;
omegMOmeg1 = omegax1;
omegMOmeg2 = omegay1;
omegMOmeg3 = omegaz1 - Omega;
} else {
if(corot2)    corot = 1;
if(VwApprox2) VwApprox = 1;
xmax = xmax2;
omegMOmeg1 = omegax2;
omegMOmeg2 = omegay2;
omegMOmeg3 = omegaz2 - Omega;
} /* end if */
xC = xCM + ecc * (xmax - xCM);
SGRID_LEVEL3_Pragma(omp parallel for)
forallpoints(box, ijk) {
double rho0, P, rhoE, drho0dhm1; /* declare them */
/* Jetzt geht's los! */
double alpha2;
double alphaP2;
double alphaP3;
double dalpha1;
double dalpha2;
double dalpha3;
double dbeta11;
double dbeta12;
double dbeta13;
double dbeta21;
double dbeta22;
double dbeta23;
double dbeta31;
double dbeta32;
double dbeta33;
double divbeta;
double divlbeta;
double divlwB;
double divwB;
double dL21;
double dL22;
double dL23;
double dlalpha1;
double dlalpha2;
double dlalpha3;
double dlh1;
double dlh2;
double dlh3;
double dLnalpha1;
double dLnalpha2;
double dLnalpha3;
double dLnalphaoh1;
double dLnalphaoh2;
double dLnalphaoh3;
double dLnalphaP1;
double dLnalphaP2;
double dLnalphaP3;
double dLnalphaPsi21;
double dLnalphaPsi22;
double dLnalphaPsi23;
double dLnalphaPsi2oh1;
double dLnalphaPsi2oh2;
double dLnalphaPsi2oh3;
double dLnalphaPsi6uz1;
double dLnalphaPsi6uz2;
double dLnalphaPsi6uz3;
double dLnalphaPsim61;
double dLnalphaPsim62;
double dLnalphaPsim63;
double dLnh1;
double dLnh2;
double dLnh3;
double dLnPsi1;
double dLnPsi2;
double dLnPsi3;
double dLnuzero1;
double dLnuzero2;
double dLnuzero3;
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
double drho01;
double drho02;
double drho03;
double drho0PLUSrho0dLnalphaoh1;
double drho0PLUSrho0dLnalphaoh2;
double drho0PLUSrho0dLnalphaoh3;
double drho0PLUSrho0dLnalphaPsi2oh1;
double drho0PLUSrho0dLnalphaPsi2oh2;
double drho0PLUSrho0dLnalphaPsi2oh3;
double drho0PLUSrho0dLnalphaPsi6uz1;
double drho0PLUSrho0dLnalphaPsi6uz2;
double drho0PLUSrho0dLnalphaPsi6uz3;
double drhof1;
double drhof2;
double drhof3;
double drhofPLUSrhofdLnalpha1;
double drhofPLUSrhofdLnalpha2;
double drhofPLUSrhofdLnalpha3;
double drhofPLUSrhofdLnalphaPsi21;
double drhofPLUSrhofdLnalphaPsi22;
double drhofPLUSrhofdLnalphaPsi23;
double drhofPLUSrhofdLnalphaPsi6uz1;
double drhofPLUSrhofdLnalphaPsi6uz2;
double drhofPLUSrhofdLnalphaPsi6uz3;
double dSigmaUp1;
double DSigmaUp1;
double dSigmaUp2;
double DSigmaUp2;
double dSigmaUp3;
double DSigmaUp3;
double duzero1;
double duzero2;
double duzero3;
double duzerosqr1;
double duzerosqr2;
double duzerosqr3;
double gdB;
double gdlB;
double h;
double h2;
double hf;
double j1;
double j2;
double j3;
double L2;
double lalpha;
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
double LBLB;
double ldL21;
double ldL22;
double ldL23;
double ldLnalpha1;
double ldLnalpha2;
double ldLnalpha3;
double ldLnalphaoh1;
double ldLnalphaoh2;
double ldLnalphaoh3;
double ldLnalphaP1;
double ldLnalphaP2;
double ldLnalphaP3;
double ldLnalphaPsi21;
double ldLnalphaPsi22;
double ldLnalphaPsi23;
double ldLnalphaPsi2oh1;
double ldLnalphaPsi2oh2;
double ldLnalphaPsi2oh3;
double ldLnalphaPsi6uz1;
double ldLnalphaPsi6uz2;
double ldLnalphaPsi6uz3;
double ldLnalphaPsim61;
double ldLnalphaPsim62;
double ldLnalphaPsim63;
double ldLnh1;
double ldLnh2;
double ldLnh3;
double ldLnPsi1;
double ldLnPsi2;
double ldLnPsi3;
double ldLnuzero1;
double ldLnuzero2;
double ldLnuzero3;
double ldrho01;
double ldrho02;
double ldrho03;
double ldrho0PLUSrho0dLnalphaoh1;
double ldrho0PLUSrho0dLnalphaoh2;
double ldrho0PLUSrho0dLnalphaoh3;
double ldrho0PLUSrho0dLnalphaPsi2oh1;
double ldrho0PLUSrho0dLnalphaPsi2oh2;
double ldrho0PLUSrho0dLnalphaPsi2oh3;
double ldrho0PLUSrho0dLnalphaPsi6uz1;
double ldrho0PLUSrho0dLnalphaPsi6uz2;
double ldrho0PLUSrho0dLnalphaPsi6uz3;
double ldrhofPLUSrhofdLnalpha1;
double ldrhofPLUSrhofdLnalpha2;
double ldrhofPLUSrhofdLnalpha3;
double ldrhofPLUSrhofdLnalphaPsi21;
double ldrhofPLUSrhofdLnalphaPsi22;
double ldrhofPLUSrhofdLnalphaPsi23;
double ldrhofPLUSrhofdLnalphaPsi6uz1;
double ldrhofPLUSrhofdLnalphaPsi6uz2;
double ldrhofPLUSrhofdLnalphaPsi6uz3;
double lduzero1;
double lduzero2;
double lduzero3;
double lduzerosqr1;
double lduzerosqr2;
double lduzerosqr3;
double lh;
double lhfuzeroPsi4;
double lhuzeroPsi4;
double lj1;
double lj2;
double lj3;
double lL2;
double LlB11;
double LlB12;
double LlB13;
double LlB22;
double LlB23;
double LlB33;
double LlBdo11;
double LlBdo12;
double LlBdo13;
double LlBdo21;
double LlBdo22;
double LlBdo23;
double LlBdo31;
double LlBdo32;
double LlBdo33;
double LlBLlB;
double lLnalpha;
double lLnh;
double lLnuzero;
double loouzerosqr;
double lrho;
double lrho0;
double lS;
double luzero;
double luzerosqr;
double lvR1;
double lvR2;
double lvR3;
double lwB1;
double lwB2;
double lwB3;
double lwBDown1;
double lwBDown2;
double lwBDown3;
double OmegaCrossR1;
double OmegaCrossR2;
double OmegaCrossR3;
double oouzerosqr;
double Psi2;
double Psi3;
double Psi4;
double Psi5;
double Psi6;
double Psi7;
double Psim1;
double Psim10;
double Psim2;
double Psim3;
double Psim4;
double Psim5;
double Psim6;
double Psim7;
double Psim8;
double Psim9;
double rdotor;
double rho;
double rhof;
double S;
double uzero;
double uzerosqr;
double vecLapB1;
double vecLapB2;
double vecLapB3;
double vecLaplB1;
double vecLaplB2;
double vecLaplB3;
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
double xMxmax1;
double xMxmax2;
double xMxmax3;
double xrdotor1;
double xrdotor2;
double xrdotor3;



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

LBLB
=
LB11*LBdo11 + LB12*(LBdo12 + LBdo21) + LB22*LBdo22 + 
  LB13*(LBdo13 + LBdo31) + LB23*(LBdo23 + LBdo32) + LB33*LBdo33
;

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

Psi3
=
Psi2*Psi[ijk]
;

Psi4
=
pow2(Psi2)
;

Psi5
=
Psi4*Psi[ijk]
;


DNS_polytrope_EoS_of_hm1(q[ijk],&rho0, &P, &rhoE, &drho0dhm1); 



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

Psim1
=
1/Psi[ijk]
;

Psim2
=
pow2(Psim1)
;

Psim3
=
Psim1*Psim2
;

Psim4
=
pow2(Psim2)
;

Psim8
=
pow2(Psim4)
;

Psim6
=
Psi2*Psim8
;

Psim5
=
Psim6*Psi[ijk]
;

Psim7
=
Psim8*Psi[ijk]
;

Psim9
=
Psim1*Psim8
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

drho01
=
drho0dhm1*dq1[ijk]
;

drho02
=
drho0dhm1*dq2[ijk]
;

drho03
=
drho0dhm1*dq3[ijk]
;

dLnPsi1
=
dPsi1[ijk]/Psi[ijk]
;

dLnPsi2
=
dPsi2[ijk]/Psi[ijk]
;

dLnPsi3
=
dPsi3[ijk]/Psi[ijk]
;

dLnh1
=
dq1[ijk]/h
;

dLnh2
=
dq2[ijk]/h
;

dLnh3
=
dq3[ijk]/h
;

dLnalphaP1
=
dalphaP1[ijk]/alphaP[ijk]
;

dLnalphaP2
=
dalphaP2[ijk]/alphaP[ijk]
;

dLnalphaP3
=
dalphaP3[ijk]/alphaP[ijk]
;

dalpha1
=
-((alphaP[ijk]*dPsi1[ijk])/Psi2) + dalphaP1[ijk]/Psi[ijk]
;

dalpha2
=
-((alphaP[ijk]*dPsi2[ijk])/Psi2) + dalphaP2[ijk]/Psi[ijk]
;

dalpha3
=
-((alphaP[ijk]*dPsi3[ijk])/Psi2) + dalphaP3[ijk]/Psi[ijk]
;

dLnalpha1
=
dalpha1/alpha[ijk]
;

dLnalpha2
=
dalpha2/alpha[ijk]
;

dLnalpha3
=
dalpha3/alpha[ijk]
;

dL21
=
-(dPsi1[ijk]*(4.*Psim5*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
          dSigmaUp3*dSigma3[ijk]) + 
       (8.*Psim9*wBDown1 + 12.*Psim7*dSigma1[ijk])*wB1[ijk] + 
       (8.*Psim9*wBDown2 + 12.*Psim7*dSigma2[ijk])*wB2[ijk] + 
       (8.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk])*wB3[ijk])) + 
  2.*(dLnh1*h2 + Psim4*(dSigmaUp1*ddSigma11[ijk] + 
        dSigmaUp2*ddSigma12[ijk] + dSigmaUp3*ddSigma13[ijk]) + 
     (Psim8*wBDown1 + Psim6*dSigma1[ijk])*dwB11[ijk] + 
     (Psim8*wBDown2 + Psim6*dSigma2[ijk])*dwB21[ijk] + 
     (Psim8*wBDown3 + Psim6*dSigma3[ijk])*dwB31[ijk] + 
     Psim6*(ddSigma11[ijk]*wB1[ijk] + ddSigma12[ijk]*wB2[ijk] + 
        ddSigma13[ijk]*wB3[ijk]))
;

dL22
=
-(dPsi2[ijk]*(4.*Psim5*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
          dSigmaUp3*dSigma3[ijk]) + 
       (8.*Psim9*wBDown1 + 12.*Psim7*dSigma1[ijk])*wB1[ijk] + 
       (8.*Psim9*wBDown2 + 12.*Psim7*dSigma2[ijk])*wB2[ijk] + 
       (8.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk])*wB3[ijk])) + 
  2.*(dLnh2*h2 + Psim4*(dSigmaUp1*ddSigma12[ijk] + 
        dSigmaUp2*ddSigma22[ijk] + dSigmaUp3*ddSigma23[ijk]) + 
     (Psim8*wBDown1 + Psim6*dSigma1[ijk])*dwB12[ijk] + 
     (Psim8*wBDown2 + Psim6*dSigma2[ijk])*dwB22[ijk] + 
     (Psim8*wBDown3 + Psim6*dSigma3[ijk])*dwB32[ijk] + 
     Psim6*(ddSigma12[ijk]*wB1[ijk] + ddSigma22[ijk]*wB2[ijk] + 
        ddSigma23[ijk]*wB3[ijk]))
;

dL23
=
-(dPsi3[ijk]*(4.*Psim5*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
          dSigmaUp3*dSigma3[ijk]) + 
       (8.*Psim9*wBDown1 + 12.*Psim7*dSigma1[ijk])*wB1[ijk] + 
       (8.*Psim9*wBDown2 + 12.*Psim7*dSigma2[ijk])*wB2[ijk] + 
       (8.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk])*wB3[ijk])) + 
  2.*(dLnh3*h2 + Psim4*(dSigmaUp1*ddSigma13[ijk] + 
        dSigmaUp2*ddSigma23[ijk] + dSigmaUp3*ddSigma33[ijk]) + 
     (Psim8*wBDown1 + Psim6*dSigma1[ijk])*dwB13[ijk] + 
     (Psim8*wBDown2 + Psim6*dSigma2[ijk])*dwB23[ijk] + 
     (Psim8*wBDown3 + Psim6*dSigma3[ijk])*dwB33[ijk] + 
     Psim6*(ddSigma13[ijk]*wB1[ijk] + ddSigma23[ijk]*wB2[ijk] + 
        ddSigma33[ijk]*wB3[ijk]))
;

duzerosqr1
=
(dL21 - 2.*L2*(dLnh1 + dalpha1/alpha[ijk]))/(alpha2*h2)
;

duzerosqr2
=
(dL22 - 2.*L2*(dLnh2 + dalpha2/alpha[ijk]))/(alpha2*h2)
;

duzerosqr3
=
(dL23 - 2.*L2*(dLnh3 + dalpha3/alpha[ijk]))/(alpha2*h2)
;

duzero1
=
(0.5*duzerosqr1)/uzero
;

duzero2
=
(0.5*duzerosqr2)/uzero
;

duzero3
=
(0.5*duzerosqr3)/uzero
;

dbeta11
=
rdotor + dB11[ijk]
;

dbeta12
=
-Omega + dB12[ijk]
;

dbeta13
=
dB13[ijk]
;

dbeta21
=
Omega + dB21[ijk]
;

dbeta22
=
rdotor + dB22[ijk]
;

dbeta23
=
dB23[ijk]
;

dbeta31
=
dB31[ijk]
;

dbeta32
=
dB32[ijk]
;

dbeta33
=
rdotor + dB33[ijk]
;

dLnalphaPsi2oh1
=
dLnalphaP1 - dLnh1 + dLnPsi1
;

dLnalphaPsi2oh2
=
dLnalphaP2 - dLnh2 + dLnPsi2
;

dLnalphaPsi2oh3
=
dLnalphaP3 - dLnh3 + dLnPsi3
;

dLnalphaoh1
=
dLnalphaP1 - dLnh1 - dLnPsi1
;

dLnalphaoh2
=
dLnalphaP2 - dLnh2 - dLnPsi2
;

dLnalphaoh3
=
dLnalphaP3 - dLnh3 - dLnPsi3
;

dLnalphaPsi6uz1
=
dLnalphaP1 + 5.*dLnPsi1 + duzero1/uzero
;

dLnalphaPsi6uz2
=
dLnalphaP2 + 5.*dLnPsi2 + duzero2/uzero
;

dLnalphaPsi6uz3
=
dLnalphaP3 + 5.*dLnPsi3 + duzero3/uzero
;

drho0PLUSrho0dLnalphaPsi2oh1
=
drho01 + dLnalphaPsi2oh1*rho0
;

drho0PLUSrho0dLnalphaPsi2oh2
=
drho02 + dLnalphaPsi2oh2*rho0
;

drho0PLUSrho0dLnalphaPsi2oh3
=
drho03 + dLnalphaPsi2oh3*rho0
;

drho0PLUSrho0dLnalphaoh1
=
drho01 + dLnalphaoh1*rho0
;

drho0PLUSrho0dLnalphaoh2
=
drho02 + dLnalphaoh2*rho0
;

drho0PLUSrho0dLnalphaoh3
=
drho03 + dLnalphaoh3*rho0
;

drho0PLUSrho0dLnalphaPsi6uz1
=
drho01 + dLnalphaPsi6uz1*rho0
;

drho0PLUSrho0dLnalphaPsi6uz2
=
drho02 + dLnalphaPsi6uz2*rho0
;

drho0PLUSrho0dLnalphaPsi6uz3
=
drho03 + dLnalphaPsi6uz3*rho0
;

divwB
=
dwB11[ijk] + dwB22[ijk] + dwB33[ijk]
;

divbeta
=
dbeta11 + dbeta22 + dbeta33
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

rho
=
alpha2*rhoE*uzerosqr + P*(-1. + alpha2*uzerosqr)
;

j1
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR1 + beta1[ijk])
;

j2
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR2 + beta2[ijk])
;

j3
=
(P + rhoE)*uzerosqr*alpha[ijk]*(vR3 + beta3[ijk])
;

S
=
3.*P + rho - rhoE
;

dLnalphaPsim61
=
dalphaP1[ijk]/alphaP[ijk] - (7.*dPsi1[ijk])/Psi[ijk]
;

dLnalphaPsim62
=
dalphaP2[ijk]/alphaP[ijk] - (7.*dPsi2[ijk])/Psi[ijk]
;

dLnalphaPsim63
=
dalphaP3[ijk]/alphaP[ijk] - (7.*dPsi3[ijk])/Psi[ijk]
;



/* conditional */
if (nonlin) {

vecLapB1
=
1.3333333333333333333*ddB111[ijk] + ddB122[ijk] + ddB133[ijk] + 
  0.33333333333333333333*(ddB212[ijk] + ddB313[ijk])
;

vecLapB2
=
ddB211[ijk] + 1.3333333333333333333*ddB222[ijk] + ddB233[ijk] + 
  0.33333333333333333333*(ddB112[ijk] + ddB323[ijk])
;

vecLapB3
=
0.33333333333333333333*(ddB113[ijk] + ddB223[ijk]) + ddB311[ijk] + 
  ddB322[ijk] + 1.3333333333333333333*ddB333[ijk]
;

FPsi[ijk]
=
Psi5*((0.03125*LBLB)/alpha2 + 6.283185307179586477*rho) + ddPsi11[ijk] + 
  ddPsi22[ijk] + ddPsi33[ijk]
;

FB1[ijk]
=
-(dLnalphaPsim61*LB11) - dLnalphaPsim62*LB12 - dLnalphaPsim63*LB13 + 
  vecLapB1 - 50.26548245743669182*j1*Psi4*alpha[ijk]
;

FB2[ijk]
=
-(dLnalphaPsim61*LB12) - dLnalphaPsim62*LB22 - dLnalphaPsim63*LB23 + 
  vecLapB2 - 50.26548245743669182*j2*Psi4*alpha[ijk]
;

FB3[ijk]
=
-(dLnalphaPsim61*LB13) - dLnalphaPsim62*LB23 - dLnalphaPsim63*LB33 + 
  vecLapB3 - 50.26548245743669182*j3*Psi4*alpha[ijk]
;

FalphaP[ijk]
=
-(Psi4*((0.21875*LBLB)/alpha2 + 3.1415926535897932385*(2.*rho + 4.*S))*
     alphaP[ijk]) + ddalphaP11[ijk] + ddalphaP22[ijk] + ddalphaP33[ijk]
;



/* conditional */
if (corot) {

FSigma[ijk]
=
Sigma[ijk]
;


} else { /* if (!corot) */


if(MATTRinside) { 

FSigma[ijk]
=
drho0PLUSrho0dLnalphaPsi2oh1*dSigmaUp1 + 
  drho0PLUSrho0dLnalphaPsi2oh2*dSigmaUp2 + 
  drho0PLUSrho0dLnalphaPsi2oh3*dSigmaUp3 - 
  h*Psi4*uzero*(drho0PLUSrho0dLnalphaPsi6uz1*beta1[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz2*beta2[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz3*beta3[ijk]) + 
  rho0*(divwB*Psim2 - divbeta*h*Psi4*uzero + ddSigma11[ijk] + 
     ddSigma22[ijk] + ddSigma33[ijk]) + 
  Psim2*(drho0PLUSrho0dLnalphaoh1*wB1[ijk] + 
     drho0PLUSrho0dLnalphaoh2*wB2[ijk] + drho0PLUSrho0dLnalphaoh3*wB3[ijk])
;


} else if(MATTRtouch) { 


if(FakeMatterOutside) { 

hf
=
1.
;



/* conditional */
if (FakeT0) {

rhof
=
-lam[ijk]
;

drhof1
=
-dlam1[ijk]
;

drhof2
=
-dlam2[ijk]
;

drhof3
=
-dlam3[ijk]
;


} else { /* if (!FakeT0) */

rhof
=
-1.
;

drhof1
=
0
;

drhof2
=
0
;

drhof3
=
0
;

}
/* if (FakeT0) */


dLnalphaPsi21
=
dLnalphaP1 + dLnPsi1
;

dLnalphaPsi22
=
dLnalphaP2 + dLnPsi2
;

dLnalphaPsi23
=
dLnalphaP3 + dLnPsi3
;

dLnalpha1
=
dLnalphaP1 - dLnPsi1
;

dLnalpha2
=
dLnalphaP2 - dLnPsi2
;

dLnalpha3
=
dLnalphaP3 - dLnPsi3
;

dLnalphaPsi6uz1
=
dLnalphaP1 + 5.*dLnPsi1 + duzero1/uzero
;

dLnalphaPsi6uz2
=
dLnalphaP2 + 5.*dLnPsi2 + duzero2/uzero
;

dLnalphaPsi6uz3
=
dLnalphaP3 + 5.*dLnPsi3 + duzero3/uzero
;

drhofPLUSrhofdLnalphaPsi21
=
drhof1 + dLnalphaPsi21*rhof
;

drhofPLUSrhofdLnalphaPsi22
=
drhof2 + dLnalphaPsi22*rhof
;

drhofPLUSrhofdLnalphaPsi23
=
drhof3 + dLnalphaPsi23*rhof
;

drhofPLUSrhofdLnalpha1
=
drhof1 + dLnalpha1*rhof
;

drhofPLUSrhofdLnalpha2
=
drhof2 + dLnalpha2*rhof
;

drhofPLUSrhofdLnalpha3
=
drhof3 + dLnalpha3*rhof
;

drhofPLUSrhofdLnalphaPsi6uz1
=
drhof1 + dLnalphaPsi6uz1*rhof
;

drhofPLUSrhofdLnalphaPsi6uz2
=
drhof2 + dLnalphaPsi6uz2*rhof
;

drhofPLUSrhofdLnalphaPsi6uz3
=
drhof3 + dLnalphaPsi6uz3*rhof
;

FSigma[ijk]
=
drhofPLUSrhofdLnalphaPsi21*dSigmaUp1 + 
  drhofPLUSrhofdLnalphaPsi22*dSigmaUp2 + 
  drhofPLUSrhofdLnalphaPsi23*dSigmaUp3 - 
  hf*Psi4*uzero*(drhofPLUSrhofdLnalphaPsi6uz1*beta1[ijk] + 
     drhofPLUSrhofdLnalphaPsi6uz2*beta2[ijk] + 
     drhofPLUSrhofdLnalphaPsi6uz3*beta3[ijk]) + 
  rhof*(divwB*Psim2 - divbeta*hf*Psi4*uzero + ddSigma11[ijk] + 
     ddSigma22[ijk] + ddSigma33[ijk]) + 
  Psim2*(drhofPLUSrhofdLnalpha1*wB1[ijk] + drhofPLUSrhofdLnalpha2*wB2[ijk] + 
     drhofPLUSrhofdLnalpha3*wB3[ijk])
;


} else if(LaplaceSigmaOutside) { 

FSigma[ijk]
=
ddSigma11[ijk] + ddSigma22[ijk] + ddSigma33[ijk]
;


} else { 

FSigma[ijk]
=
dddSigmadlam3[ijk] + 2.*ddSigmadlam2[ijk] + dSigmadlam[ijk]
;


} /* end Outside cases */ 


} else { 

FSigma[ijk]
=
Sigma[ijk]
;


} 

}
/* if (FakeT0) */


FPsi[ijk]
=
CoordFac[ijk]*FPsi[ijk]
;

FB1[ijk]
=
CoordFac[ijk]*FB1[ijk]
;

FB2[ijk]
=
CoordFac[ijk]*FB2[ijk]
;

FB3[ijk]
=
CoordFac[ijk]*FB3[ijk]
;

FalphaP[ijk]
=
CoordFac[ijk]*FalphaP[ijk]
;

FSigma[ijk]
=
CoordFac[ijk]*FSigma[ijk]
;


} else { /* if (!FakeT0) */

alphaP2
=
pow2(alphaP[ijk])
;

alphaP3
=
alphaP2*alphaP[ijk]
;

Psi6
=
Psi2*Psi4
;

Psi7
=
Psi3*Psi4
;

gdlB
=
dlB11[ijk] + dlB22[ijk] + dlB33[ijk]
;

LlB11
=
-0.66666666666666666667*gdlB + 2.*dlB11[ijk]
;

LlB12
=
dlB12[ijk] + dlB21[ijk]
;

LlB13
=
dlB13[ijk] + dlB31[ijk]
;

LlB22
=
-0.66666666666666666667*gdlB + 2.*dlB22[ijk]
;

LlB23
=
dlB23[ijk] + dlB32[ijk]
;

LlB33
=
-0.66666666666666666667*gdlB + 2.*dlB33[ijk]
;

LlBdo11
=
LlB11
;

LlBdo12
=
LlB12
;

LlBdo13
=
LlB13
;

LlBdo21
=
LlB12
;

LlBdo22
=
LlB22
;

LlBdo23
=
LlB23
;

LlBdo31
=
LlB13
;

LlBdo32
=
LlB23
;

LlBdo33
=
LlB33
;

LlBLlB
=
LlB11*LlBdo11 + LlB12*(LlBdo12 + LlBdo21) + LlB22*LlBdo22 + 
  LlB13*(LlBdo13 + LlBdo31) + LlB23*(LlBdo23 + LlBdo32) + LlB33*LlBdo33
;

vecLaplB1
=
1.3333333333333333333*ddlB111[ijk] + ddlB122[ijk] + ddlB133[ijk] + 
  0.33333333333333333333*(ddlB212[ijk] + ddlB313[ijk])
;

vecLaplB2
=
ddlB211[ijk] + 1.3333333333333333333*ddlB222[ijk] + ddlB233[ijk] + 
  0.33333333333333333333*(ddlB112[ijk] + ddlB323[ijk])
;

vecLaplB3
=
0.33333333333333333333*(ddlB113[ijk] + ddlB223[ijk]) + ddlB311[ijk] + 
  ddlB322[ijk] + 1.3333333333333333333*ddlB333[ijk]
;

lalpha
=
-((alphaP[ijk]*lPsi[ijk])/Psi2) + lalphaP[ijk]/Psi[ijk]
;



/* conditional */
if (corot) {

lvR1
=
0
;

lvR2
=
0
;

lvR3
=
0
;

loouzerosqr
=
2.*(lalpha*alpha[ijk] - Psi4*((vR1 + beta1[ijk])*(lvR1 + lB1[ijk]) + 
        (vR2 + beta2[ijk])*(lvR2 + lB2[ijk]) + 
        (vR3 + beta3[ijk])*(lvR3 + lB3[ijk]))) - 
  Psi3*lPsi[ijk]*(8.*(vR1*beta1[ijk] + vR2*beta2[ijk] + vR3*beta3[ijk]) + 
     4.*(pow2(vR1) + pow2(vR2) + pow2(vR3) + pow2(beta1[ijk]) + 
        pow2(beta2[ijk]) + pow2(beta3[ijk])))
;

luzerosqr
=
-(loouzerosqr*pow2(uzerosqr))
;


} else { /* if (!corot) */

lLnalpha
=
lalpha/alpha[ijk]
;

dlalpha1
=
(2.*alphaP[ijk]*dPsi1[ijk]*lPsi[ijk])/Psi3 - 
  (alphaP[ijk]*dlPsi1[ijk] + dPsi1[ijk]*lalphaP[ijk] + 
     dalphaP1[ijk]*lPsi[ijk])/Psi2 + dlalphaP1[ijk]/Psi[ijk]
;

dlalpha2
=
(2.*alphaP[ijk]*dPsi2[ijk]*lPsi[ijk])/Psi3 - 
  (alphaP[ijk]*dlPsi2[ijk] + dPsi2[ijk]*lalphaP[ijk] + 
     dalphaP2[ijk]*lPsi[ijk])/Psi2 + dlalphaP2[ijk]/Psi[ijk]
;

dlalpha3
=
(2.*alphaP[ijk]*dPsi3[ijk]*lPsi[ijk])/Psi3 - 
  (alphaP[ijk]*dlPsi3[ijk] + dPsi3[ijk]*lalphaP[ijk] + 
     dalphaP3[ijk]*lPsi[ijk])/Psi2 + dlalphaP3[ijk]/Psi[ijk]
;

ldLnalpha1
=
-((dalpha1*lalpha)/alpha2) + dlalpha1/alpha[ijk]
;

ldLnalpha2
=
-((dalpha2*lalpha)/alpha2) + dlalpha2/alpha[ijk]
;

ldLnalpha3
=
-((dalpha3*lalpha)/alpha2) + dlalpha3/alpha[ijk]
;

lh
=
0
;

lLnh
=
0
;

dlh1
=
0
;

dlh2
=
0
;

dlh3
=
0
;

ldLnh1
=
0
;

ldLnh2
=
0
;

ldLnh3
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

lL2
=
-(lPsi[ijk]*(4.*Psim5*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
          dSigmaUp3*dSigma3[ijk]) + 
       (8.*Psim9*wBDown1 + 12.*Psim7*dSigma1[ijk])*wB1[ijk] + 
       (8.*Psim9*wBDown2 + 12.*Psim7*dSigma2[ijk])*wB2[ijk] + 
       (8.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk])*wB3[ijk])) + 
  2.*(h2*lLnh + Psim8*(lwB1*wBDown1 + lwB2*wBDown2 + lwB3*wBDown3) + 
     Psim4*(dSigmaUp1*dlSigma1[ijk] + dSigmaUp2*dlSigma2[ijk] + 
        dSigmaUp3*dlSigma3[ijk]) + 
     Psim6*(lwB1*dSigma1[ijk] + lwB2*dSigma2[ijk] + lwB3*dSigma3[ijk] + 
        dlSigma1[ijk]*wB1[ijk] + dlSigma2[ijk]*wB2[ijk] + 
        dlSigma3[ijk]*wB3[ijk]))
;

luzerosqr
=
(lL2 - 2.*L2*(lLnh + lalpha/alpha[ijk]))/(alpha2*h2)
;

luzero
=
(0.5*luzerosqr)/uzero
;

lwBDown1
=
lwB1
;

lwBDown2
=
lwB2
;

lwBDown3
=
lwB3
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

Psim10
=
Psim1*Psim9
;

ldL21
=
-4.*Psim5*dlPsi1[ijk]*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
     dSigmaUp3*dSigma3[ijk]) - 
  16.*Psim9*((lwB1*wBDown1 + lwB2*wBDown2)*dPsi1[ijk] + 
     wBDown1*dwB11[ijk]*lPsi[ijk]) + 
  (72.*Psim10*wBDown1*dPsi1[ijk]*lPsi[ijk] + 
     dSigma1[ijk]*(-12.*Psim7*dlPsi1[ijk] + 84.*Psim8*dPsi1[ijk]*lPsi[ijk])\
)*wB1[ijk] + (72.*Psim10*wBDown2*dPsi1[ijk]*lPsi[ijk] + 
     dSigma2[ijk]*(-12.*Psim7*dlPsi1[ijk] + 84.*Psim8*dPsi1[ijk]*lPsi[ijk])\
)*wB2[ijk] + (72.*Psim10*wBDown3*dPsi1[ijk]*lPsi[ijk] + 
     dSigma3[ijk]*(-12.*Psim7*dlPsi1[ijk] + 84.*Psim8*dPsi1[ijk]*lPsi[ijk])\
)*wB3[ijk] - 8.*(Psim5*(dSigmaUp1*ddSigma11[ijk] + 
        dSigmaUp2*ddSigma12[ijk])*lPsi[ijk] + 
     Psim9*dlPsi1[ijk]*(wBDown1*wB1[ijk] + wBDown2*wB2[ijk] + 
        wBDown3*wB3[ijk])) + 2.*
   (dlh1*h + dLnh1*h2*lLnh + Psim4*
      (dSigmaUp1*ddlSigma11[ijk] + dSigmaUp2*ddlSigma12[ijk] + 
        dSigmaUp3*ddlSigma13[ijk] + dlSigmaUp1*ddSigma11[ijk] + 
        dlSigmaUp2*ddSigma12[ijk] + dlSigmaUp3*ddSigma13[ijk]) + 
     Psim8*(dlwB11*wBDown1 + dlwB21*wBDown2 + dlwB31*wBDown3 + 
        lwBDown1*dwB11[ijk] + lwBDown2*dwB21[ijk] + lwBDown3*dwB31[ijk]) + 
     Psim6*(lwB1*ddSigma11[ijk] + lwB2*ddSigma12[ijk] + 
        lwB3*ddSigma13[ijk] + dlwB11*dSigma1[ijk] + dlwB21*dSigma2[ijk] + 
        dlwB31*dSigma3[ijk] + dlSigma1[ijk]*dwB11[ijk] + 
        dlSigma2[ijk]*dwB21[ijk] + dlSigma3[ijk]*dwB31[ijk] + 
        ddlSigma11[ijk]*wB1[ijk] + ddlSigma12[ijk]*wB2[ijk] + 
        ddlSigma13[ijk]*wB3[ijk])) + 
  lPsi[ijk]*(dSigmaUp3*(-8.*Psim5*ddSigma13[ijk] + 
        20.*Psim6*dPsi1[ijk]*dSigma3[ijk]) + 
     dSigma2[ijk]*(20.*dSigmaUp2*Psim6*dPsi1[ijk] - 12.*Psim7*dwB21[ijk]) - 
     16.*Psim9*(wBDown2*dwB21[ijk] + wBDown3*dwB31[ijk]) - 
     12.*Psim7*(dSigma1[ijk]*dwB11[ijk] + dSigma3[ijk]*dwB31[ijk] + 
        ddSigma11[ijk]*wB1[ijk] + ddSigma12[ijk]*wB2[ijk] + 
        ddSigma13[ijk]*wB3[ijk])) + 
  dPsi1[ijk]*(-8.*Psim5*(dSigmaUp1*dlSigma1[ijk] + 
        dSigmaUp2*dlSigma2[ijk] + dSigmaUp3*dlSigma3[ijk]) - 
     lwB3*(16.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk]) + 
     20.*dSigmaUp1*Psim6*dSigma1[ijk]*lPsi[ijk] - 
     12.*Psim7*(lwB1*dSigma1[ijk] + lwB2*dSigma2[ijk] + 
        dlSigma1[ijk]*wB1[ijk] + dlSigma2[ijk]*wB2[ijk] + 
        dlSigma3[ijk]*wB3[ijk]))
;

ldL22
=
-4.*Psim5*dlPsi2[ijk]*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
     dSigmaUp3*dSigma3[ijk]) - 
  16.*Psim9*((lwB1*wBDown1 + lwB2*wBDown2)*dPsi2[ijk] + 
     wBDown1*dwB12[ijk]*lPsi[ijk]) + 
  (72.*Psim10*wBDown1*dPsi2[ijk]*lPsi[ijk] + 
     dSigma1[ijk]*(-12.*Psim7*dlPsi2[ijk] + 84.*Psim8*dPsi2[ijk]*lPsi[ijk])\
)*wB1[ijk] + (72.*Psim10*wBDown2*dPsi2[ijk]*lPsi[ijk] + 
     dSigma2[ijk]*(-12.*Psim7*dlPsi2[ijk] + 84.*Psim8*dPsi2[ijk]*lPsi[ijk])\
)*wB2[ijk] + (72.*Psim10*wBDown3*dPsi2[ijk]*lPsi[ijk] + 
     dSigma3[ijk]*(-12.*Psim7*dlPsi2[ijk] + 84.*Psim8*dPsi2[ijk]*lPsi[ijk])\
)*wB3[ijk] - 8.*(Psim5*(dSigmaUp1*ddSigma12[ijk] + 
        dSigmaUp2*ddSigma22[ijk])*lPsi[ijk] + 
     Psim9*dlPsi2[ijk]*(wBDown1*wB1[ijk] + wBDown2*wB2[ijk] + 
        wBDown3*wB3[ijk])) + 2.*
   (dlh2*h + dLnh2*h2*lLnh + Psim4*
      (dSigmaUp1*ddlSigma12[ijk] + dSigmaUp2*ddlSigma22[ijk] + 
        dSigmaUp3*ddlSigma23[ijk] + dlSigmaUp1*ddSigma12[ijk] + 
        dlSigmaUp2*ddSigma22[ijk] + dlSigmaUp3*ddSigma23[ijk]) + 
     Psim8*(dlwB12*wBDown1 + dlwB22*wBDown2 + dlwB32*wBDown3 + 
        lwBDown1*dwB12[ijk] + lwBDown2*dwB22[ijk] + lwBDown3*dwB32[ijk]) + 
     Psim6*(lwB1*ddSigma12[ijk] + lwB2*ddSigma22[ijk] + 
        lwB3*ddSigma23[ijk] + dlwB12*dSigma1[ijk] + dlwB22*dSigma2[ijk] + 
        dlwB32*dSigma3[ijk] + dlSigma1[ijk]*dwB12[ijk] + 
        dlSigma2[ijk]*dwB22[ijk] + dlSigma3[ijk]*dwB32[ijk] + 
        ddlSigma12[ijk]*wB1[ijk] + ddlSigma22[ijk]*wB2[ijk] + 
        ddlSigma23[ijk]*wB3[ijk])) + 
  lPsi[ijk]*(dSigmaUp3*(-8.*Psim5*ddSigma23[ijk] + 
        20.*Psim6*dPsi2[ijk]*dSigma3[ijk]) + 
     dSigma2[ijk]*(20.*dSigmaUp2*Psim6*dPsi2[ijk] - 12.*Psim7*dwB22[ijk]) - 
     16.*Psim9*(wBDown2*dwB22[ijk] + wBDown3*dwB32[ijk]) - 
     12.*Psim7*(dSigma1[ijk]*dwB12[ijk] + dSigma3[ijk]*dwB32[ijk] + 
        ddSigma12[ijk]*wB1[ijk] + ddSigma22[ijk]*wB2[ijk] + 
        ddSigma23[ijk]*wB3[ijk])) + 
  dPsi2[ijk]*(-8.*Psim5*(dSigmaUp1*dlSigma1[ijk] + 
        dSigmaUp2*dlSigma2[ijk] + dSigmaUp3*dlSigma3[ijk]) - 
     lwB3*(16.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk]) + 
     20.*dSigmaUp1*Psim6*dSigma1[ijk]*lPsi[ijk] - 
     12.*Psim7*(lwB1*dSigma1[ijk] + lwB2*dSigma2[ijk] + 
        dlSigma1[ijk]*wB1[ijk] + dlSigma2[ijk]*wB2[ijk] + 
        dlSigma3[ijk]*wB3[ijk]))
;

ldL23
=
-4.*Psim5*dlPsi3[ijk]*(dSigmaUp1*dSigma1[ijk] + dSigmaUp2*dSigma2[ijk] + 
     dSigmaUp3*dSigma3[ijk]) - 
  16.*Psim9*((lwB1*wBDown1 + lwB2*wBDown2)*dPsi3[ijk] + 
     wBDown1*dwB13[ijk]*lPsi[ijk]) + 
  (72.*Psim10*wBDown1*dPsi3[ijk]*lPsi[ijk] + 
     dSigma1[ijk]*(-12.*Psim7*dlPsi3[ijk] + 84.*Psim8*dPsi3[ijk]*lPsi[ijk])\
)*wB1[ijk] + (72.*Psim10*wBDown2*dPsi3[ijk]*lPsi[ijk] + 
     dSigma2[ijk]*(-12.*Psim7*dlPsi3[ijk] + 84.*Psim8*dPsi3[ijk]*lPsi[ijk])\
)*wB2[ijk] + (72.*Psim10*wBDown3*dPsi3[ijk]*lPsi[ijk] + 
     dSigma3[ijk]*(-12.*Psim7*dlPsi3[ijk] + 84.*Psim8*dPsi3[ijk]*lPsi[ijk])\
)*wB3[ijk] - 8.*(Psim5*(dSigmaUp1*ddSigma13[ijk] + 
        dSigmaUp2*ddSigma23[ijk])*lPsi[ijk] + 
     Psim9*dlPsi3[ijk]*(wBDown1*wB1[ijk] + wBDown2*wB2[ijk] + 
        wBDown3*wB3[ijk])) + 2.*
   (dlh3*h + dLnh3*h2*lLnh + Psim4*
      (dSigmaUp1*ddlSigma13[ijk] + dSigmaUp2*ddlSigma23[ijk] + 
        dSigmaUp3*ddlSigma33[ijk] + dlSigmaUp1*ddSigma13[ijk] + 
        dlSigmaUp2*ddSigma23[ijk] + dlSigmaUp3*ddSigma33[ijk]) + 
     Psim8*(dlwB13*wBDown1 + dlwB23*wBDown2 + dlwB33*wBDown3 + 
        lwBDown1*dwB13[ijk] + lwBDown2*dwB23[ijk] + lwBDown3*dwB33[ijk]) + 
     Psim6*(lwB1*ddSigma13[ijk] + lwB2*ddSigma23[ijk] + 
        lwB3*ddSigma33[ijk] + dlwB13*dSigma1[ijk] + dlwB23*dSigma2[ijk] + 
        dlwB33*dSigma3[ijk] + dlSigma1[ijk]*dwB13[ijk] + 
        dlSigma2[ijk]*dwB23[ijk] + dlSigma3[ijk]*dwB33[ijk] + 
        ddlSigma13[ijk]*wB1[ijk] + ddlSigma23[ijk]*wB2[ijk] + 
        ddlSigma33[ijk]*wB3[ijk])) + 
  lPsi[ijk]*(dSigmaUp3*(-8.*Psim5*ddSigma33[ijk] + 
        20.*Psim6*dPsi3[ijk]*dSigma3[ijk]) + 
     dSigma2[ijk]*(20.*dSigmaUp2*Psim6*dPsi3[ijk] - 12.*Psim7*dwB23[ijk]) - 
     16.*Psim9*(wBDown2*dwB23[ijk] + wBDown3*dwB33[ijk]) - 
     12.*Psim7*(dSigma1[ijk]*dwB13[ijk] + dSigma3[ijk]*dwB33[ijk] + 
        ddSigma13[ijk]*wB1[ijk] + ddSigma23[ijk]*wB2[ijk] + 
        ddSigma33[ijk]*wB3[ijk])) + 
  dPsi3[ijk]*(-8.*Psim5*(dSigmaUp1*dlSigma1[ijk] + 
        dSigmaUp2*dlSigma2[ijk] + dSigmaUp3*dlSigma3[ijk]) - 
     lwB3*(16.*Psim9*wBDown3 + 12.*Psim7*dSigma3[ijk]) + 
     20.*dSigmaUp1*Psim6*dSigma1[ijk]*lPsi[ijk] - 
     12.*Psim7*(lwB1*dSigma1[ijk] + lwB2*dSigma2[ijk] + 
        dlSigma1[ijk]*wB1[ijk] + dlSigma2[ijk]*wB2[ijk] + 
        dlSigma3[ijk]*wB3[ijk]))
;

lduzerosqr1
=
(ldL21 - 2.*(L2*(ldLnalpha1 + ldLnh1) + (dLnalpha1 + dLnh1)*lL2))/
   (alpha2*h2) - 2.*duzerosqr1*(lLnalpha + lLnh)
;

lduzerosqr2
=
(ldL22 - 2.*(L2*(ldLnalpha2 + ldLnh2) + (dLnalpha2 + dLnh2)*lL2))/
   (alpha2*h2) - 2.*duzerosqr2*(lLnalpha + lLnh)
;

lduzerosqr3
=
(ldL23 - 2.*(L2*(ldLnalpha3 + ldLnh3) + (dLnalpha3 + dLnh3)*lL2))/
   (alpha2*h2) - 2.*duzerosqr3*(lLnalpha + lLnh)
;

lvR1
=
-lB1[ijk] - (lh*(dSigma1[ijk] + wB1[ijk]))/(h2*uzero) + 
  ((lwB1 + dlSigma1[ijk])/uzero - 
     (luzero*(dSigma1[ijk] + wB1[ijk]))/uzerosqr)/h
;

lvR2
=
-lB2[ijk] - (lh*(dSigma2[ijk] + wB2[ijk]))/(h2*uzero) + 
  ((lwB2 + dlSigma2[ijk])/uzero - 
     (luzero*(dSigma2[ijk] + wB2[ijk]))/uzerosqr)/h
;

lvR3
=
-lB3[ijk] - (lh*(dSigma3[ijk] + wB3[ijk]))/(h2*uzero) + 
  ((lwB3 + dlSigma3[ijk])/uzero - 
     (luzero*(dSigma3[ijk] + wB3[ijk]))/uzerosqr)/h
;

divlwB
=
dlwB11 + dlwB22 + dlwB33
;

divlbeta
=
dlB11[ijk] + dlB22[ijk] + dlB33[ijk]
;

lrho0
=
0
;

ldrho01
=
0
;

ldrho02
=
0
;

ldrho03
=
0
;

ldLnalphaP1
=
dlalphaP1[ijk]/alphaP[ijk] - (dalphaP1[ijk]*lalphaP[ijk])/alphaP2
;

ldLnalphaP2
=
dlalphaP2[ijk]/alphaP[ijk] - (dalphaP2[ijk]*lalphaP[ijk])/alphaP2
;

ldLnalphaP3
=
dlalphaP3[ijk]/alphaP[ijk] - (dalphaP3[ijk]*lalphaP[ijk])/alphaP2
;

ldLnPsi1
=
Psim1*dlPsi1[ijk] - Psim2*dPsi1[ijk]*lPsi[ijk]
;

ldLnPsi2
=
Psim1*dlPsi2[ijk] - Psim2*dPsi2[ijk]*lPsi[ijk]
;

ldLnPsi3
=
Psim1*dlPsi3[ijk] - Psim2*dPsi3[ijk]*lPsi[ijk]
;

}
/* if (corot) */


lrho
=
(P + rhoE)*(alpha2*luzerosqr + 2.*lalpha*uzerosqr*alpha[ijk])
;

lj1
=
luzerosqr*alpha[ijk]*((P + rhoE)*vR1 + P*beta1[ijk]) + 
  rhoE*((lalpha*uzerosqr + luzerosqr*alpha[ijk])*beta1[ijk] + 
     uzerosqr*alpha[ijk]*lB1[ijk]) + 
  uzerosqr*(lalpha*((P + rhoE)*vR1 + P*beta1[ijk]) + 
     alpha[ijk]*(lvR1*(P + rhoE) + P*lB1[ijk]))
;

lj2
=
luzerosqr*alpha[ijk]*((P + rhoE)*vR2 + P*beta2[ijk]) + 
  rhoE*((lalpha*uzerosqr + luzerosqr*alpha[ijk])*beta2[ijk] + 
     uzerosqr*alpha[ijk]*lB2[ijk]) + 
  uzerosqr*(lalpha*((P + rhoE)*vR2 + P*beta2[ijk]) + 
     alpha[ijk]*(lvR2*(P + rhoE) + P*lB2[ijk]))
;

lj3
=
luzerosqr*alpha[ijk]*((P + rhoE)*vR3 + P*beta3[ijk]) + 
  rhoE*((lalpha*uzerosqr + luzerosqr*alpha[ijk])*beta3[ijk] + 
     uzerosqr*alpha[ijk]*lB3[ijk]) + 
  uzerosqr*(lalpha*((P + rhoE)*vR3 + P*beta3[ijk]) + 
     alpha[ijk]*(lvR3*(P + rhoE) + P*lB3[ijk]))
;

lS
=
lrho
;

ldLnalphaPsim61
=
dlalphaP1[ijk]/alphaP[ijk] - (dalphaP1[ijk]*lalphaP[ijk])/alphaP2 + 
  (7.*dPsi1[ijk]*lPsi[ijk])/Psi2 - (7.*dlPsi1[ijk])/Psi[ijk]
;

ldLnalphaPsim62
=
dlalphaP2[ijk]/alphaP[ijk] - (dalphaP2[ijk]*lalphaP[ijk])/alphaP2 + 
  (7.*dPsi2[ijk]*lPsi[ijk])/Psi2 - (7.*dlPsi2[ijk])/Psi[ijk]
;

ldLnalphaPsim63
=
dlalphaP3[ijk]/alphaP[ijk] - (dalphaP3[ijk]*lalphaP[ijk])/alphaP2 + 
  (7.*dPsi3[ijk]*lPsi[ijk])/Psi2 - (7.*dlPsi3[ijk])/Psi[ijk]
;

FlPsi[ijk]
=
((0.0625*(LBdo11*LlB11 + (LBdo12 + LBdo21)*LlB12 + 
          (LBdo13 + LBdo31)*LlB13 + LBdo22*LlB22 + 
          (LBdo23 + LBdo32)*LlB23 + LBdo33*LlB33))/alpha2 + 
     6.283185307179586477*lrho)*Psi5 + ddlPsi11[ijk] + ddlPsi22[ijk] + 
  ddlPsi33[ijk] + 31.415926535897932385*Psi4*rho*lPsi[ijk] + 
  LBLB*((-0.0625*Psi7*lalphaP[ijk])/alphaP3 + 
     (0.21875*Psi6*lPsi[ijk])/alphaP2)
;

FlB1[ijk]
=
-(LB11*ldLnalphaPsim61) - LB12*ldLnalphaPsim62 - LB13*ldLnalphaPsim63 - 
  dLnalphaPsim61*LlB11 - dLnalphaPsim62*LlB12 - dLnalphaPsim63*LlB13 + 
  vecLaplB1 - 3.1415926535897932385*
   (16.*(lj1*Psi4*alpha[ijk] + j1*Psi3*lalphaP[ijk]) + 
     48.*j1*Psi2*alphaP[ijk]*lPsi[ijk])
;

FlB2[ijk]
=
-(LB12*ldLnalphaPsim61) - LB22*ldLnalphaPsim62 - LB23*ldLnalphaPsim63 - 
  dLnalphaPsim61*LlB12 - dLnalphaPsim62*LlB22 - dLnalphaPsim63*LlB23 + 
  vecLaplB2 - 3.1415926535897932385*
   (16.*(lj2*Psi4*alpha[ijk] + j2*Psi3*lalphaP[ijk]) + 
     48.*j2*Psi2*alphaP[ijk]*lPsi[ijk])
;

FlB3[ijk]
=
-(LB13*ldLnalphaPsim61) - LB23*ldLnalphaPsim62 - LB33*ldLnalphaPsim63 - 
  dLnalphaPsim61*LlB13 - dLnalphaPsim62*LlB23 - dLnalphaPsim63*LlB33 + 
  vecLaplB3 - 3.1415926535897932385*
   (16.*(lj3*Psi4*alpha[ijk] + j3*Psi3*lalphaP[ijk]) + 
     48.*j3*Psi2*alphaP[ijk]*lPsi[ijk])
;

FlalphaP[ijk]
=
ddlalphaP11[ijk] + ddlalphaP22[ijk] + ddlalphaP33[ijk] + 
  ((0.21875*LBLB*Psi6)/alphaP2 - 12.566370614359172954*Psi4*S)*
   lalphaP[ijk] - Psi4*(((0.4375*
           (LBdo11*LlB11 + (LBdo12 + LBdo21)*LlB12 + 
             (LBdo13 + LBdo31)*LlB13 + LBdo22*LlB22 + 
             (LBdo23 + LBdo32)*LlB23 + LBdo33*LlB33))/alpha2 + 
        3.1415926535897932385*(2.*lrho + 4.*lS))*alphaP[ijk] + 
     6.283185307179586477*rho*lalphaP[ijk]) - 
  ((1.3125*LBLB*Psi5)/alphaP2 + 3.1415926535897932385*Psi3*(8.*rho + 16.*S))*
   alphaP[ijk]*lPsi[ijk]
;



/* conditional */
if (corot) {

FlSigma[ijk]
=
lSigma[ijk]
;


} else { /* if (!corot) */


if(MATTRinside) { 

lLnuzero
=
luzero/uzero
;

dLnuzero1
=
(0.5*duzerosqr1)/uzerosqr
;

dLnuzero2
=
(0.5*duzerosqr2)/uzerosqr
;

dLnuzero3
=
(0.5*duzerosqr3)/uzerosqr
;

lduzero1
=
(0.5*lduzerosqr1 - (0.25*duzerosqr1*luzerosqr)/uzerosqr)/uzero
;

lduzero2
=
(0.5*lduzerosqr2 - (0.25*duzerosqr2*luzerosqr)/uzerosqr)/uzero
;

lduzero3
=
(0.5*lduzerosqr3 - (0.25*duzerosqr3*luzerosqr)/uzerosqr)/uzero
;

ldLnuzero1
=
-(dLnuzero1*lLnuzero) + lduzero1/uzero
;

ldLnuzero2
=
-(dLnuzero2*lLnuzero) + lduzero2/uzero
;

ldLnuzero3
=
-(dLnuzero3*lLnuzero) + lduzero3/uzero
;

lhuzeroPsi4
=
Psi4*(h*luzero + lh*uzero) + 4.*h*Psi3*uzero*lPsi[ijk]
;

ldLnalphaPsi2oh1
=
ldLnalphaP1 - ldLnh1 + ldLnPsi1
;

ldLnalphaPsi2oh2
=
ldLnalphaP2 - ldLnh2 + ldLnPsi2
;

ldLnalphaPsi2oh3
=
ldLnalphaP3 - ldLnh3 + ldLnPsi3
;

ldLnalphaoh1
=
ldLnalphaP1 - ldLnh1 - ldLnPsi1
;

ldLnalphaoh2
=
ldLnalphaP2 - ldLnh2 - ldLnPsi2
;

ldLnalphaoh3
=
ldLnalphaP3 - ldLnh3 - ldLnPsi3
;

ldLnalphaPsi6uz1
=
ldLnalphaP1 + 5.*ldLnPsi1 + ldLnuzero1
;

ldLnalphaPsi6uz2
=
ldLnalphaP2 + 5.*ldLnPsi2 + ldLnuzero2
;

ldLnalphaPsi6uz3
=
ldLnalphaP3 + 5.*ldLnPsi3 + ldLnuzero3
;

ldrho0PLUSrho0dLnalphaPsi2oh1
=
ldrho01 + dLnalphaPsi2oh1*lrho0 + ldLnalphaPsi2oh1*rho0
;

ldrho0PLUSrho0dLnalphaPsi2oh2
=
ldrho02 + dLnalphaPsi2oh2*lrho0 + ldLnalphaPsi2oh2*rho0
;

ldrho0PLUSrho0dLnalphaPsi2oh3
=
ldrho03 + dLnalphaPsi2oh3*lrho0 + ldLnalphaPsi2oh3*rho0
;

ldrho0PLUSrho0dLnalphaoh1
=
ldrho01 + dLnalphaoh1*lrho0 + ldLnalphaoh1*rho0
;

ldrho0PLUSrho0dLnalphaoh2
=
ldrho02 + dLnalphaoh2*lrho0 + ldLnalphaoh2*rho0
;

ldrho0PLUSrho0dLnalphaoh3
=
ldrho03 + dLnalphaoh3*lrho0 + ldLnalphaoh3*rho0
;

ldrho0PLUSrho0dLnalphaPsi6uz1
=
ldrho01 + dLnalphaPsi6uz1*lrho0 + ldLnalphaPsi6uz1*rho0
;

ldrho0PLUSrho0dLnalphaPsi6uz2
=
ldrho02 + dLnalphaPsi6uz2*lrho0 + ldLnalphaPsi6uz2*rho0
;

ldrho0PLUSrho0dLnalphaPsi6uz3
=
ldrho03 + dLnalphaPsi6uz3*lrho0 + ldLnalphaPsi6uz3*rho0
;

FlSigma[ijk]
=
dlSigmaUp1*drho0PLUSrho0dLnalphaPsi2oh1 + 
  dlSigmaUp2*drho0PLUSrho0dLnalphaPsi2oh2 + 
  dlSigmaUp3*drho0PLUSrho0dLnalphaPsi2oh3 + 
  dSigmaUp1*ldrho0PLUSrho0dLnalphaPsi2oh1 + 
  dSigmaUp2*ldrho0PLUSrho0dLnalphaPsi2oh2 + 
  dSigmaUp3*ldrho0PLUSrho0dLnalphaPsi2oh3 - 
  divbeta*(lhuzeroPsi4*rho0 + h*lrho0*Psi4*uzero) - 
  lhuzeroPsi4*(drho0PLUSrho0dLnalphaPsi6uz1*beta1[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz2*beta2[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz3*beta3[ijk]) + 
  lrho0*(ddSigma11[ijk] + ddSigma22[ijk] + ddSigma33[ijk]) - 
  h*Psi4*uzero*(divlbeta*rho0 + ldrho0PLUSrho0dLnalphaPsi6uz1*beta1[ijk] + 
     ldrho0PLUSrho0dLnalphaPsi6uz2*beta2[ijk] + 
     ldrho0PLUSrho0dLnalphaPsi6uz3*beta3[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz1*lB1[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz2*lB2[ijk] + 
     drho0PLUSrho0dLnalphaPsi6uz3*lB3[ijk]) + 
  rho0*(ddlSigma11[ijk] + ddlSigma22[ijk] + ddlSigma33[ijk] - 
     2.*divwB*Psim3*lPsi[ijk]) - 
  2.*Psim3*lPsi[ijk]*(drho0PLUSrho0dLnalphaoh1*wB1[ijk] + 
     drho0PLUSrho0dLnalphaoh2*wB2[ijk] + drho0PLUSrho0dLnalphaoh3*wB3[ijk]) \
+ Psim2*(divwB*lrho0 + drho0PLUSrho0dLnalphaoh1*lwB1 + 
     drho0PLUSrho0dLnalphaoh2*lwB2 + drho0PLUSrho0dLnalphaoh3*lwB3 + 
     divlwB*rho0 + ldrho0PLUSrho0dLnalphaoh1*wB1[ijk] + 
     ldrho0PLUSrho0dLnalphaoh2*wB2[ijk] + ldrho0PLUSrho0dLnalphaoh3*wB3[ijk])
;


} else if(MATTRtouch) { 


if(FakeMatterOutside) { 

hf
=
1.
;



/* conditional */
if (FakeT0) {

rhof
=
-lam[ijk]
;

drhof1
=
-dlam1[ijk]
;

drhof2
=
-dlam2[ijk]
;

drhof3
=
-dlam3[ijk]
;


} else { /* if (!FakeT0) */

rhof
=
-1.
;

drhof1
=
0
;

drhof2
=
0
;

drhof3
=
0
;

}
/* if (FakeT0) */


dLnalphaPsi21
=
dLnalphaP1 + dLnPsi1
;

dLnalphaPsi22
=
dLnalphaP2 + dLnPsi2
;

dLnalphaPsi23
=
dLnalphaP3 + dLnPsi3
;

dLnalpha1
=
dLnalphaP1 - dLnPsi1
;

dLnalpha2
=
dLnalphaP2 - dLnPsi2
;

dLnalpha3
=
dLnalphaP3 - dLnPsi3
;

dLnalphaPsi6uz1
=
dLnalphaP1 + 5.*dLnPsi1 + duzero1/uzero
;

dLnalphaPsi6uz2
=
dLnalphaP2 + 5.*dLnPsi2 + duzero2/uzero
;

dLnalphaPsi6uz3
=
dLnalphaP3 + 5.*dLnPsi3 + duzero3/uzero
;

drhofPLUSrhofdLnalphaPsi21
=
drhof1 + dLnalphaPsi21*rhof
;

drhofPLUSrhofdLnalphaPsi22
=
drhof2 + dLnalphaPsi22*rhof
;

drhofPLUSrhofdLnalphaPsi23
=
drhof3 + dLnalphaPsi23*rhof
;

drhofPLUSrhofdLnalpha1
=
drhof1 + dLnalpha1*rhof
;

drhofPLUSrhofdLnalpha2
=
drhof2 + dLnalpha2*rhof
;

drhofPLUSrhofdLnalpha3
=
drhof3 + dLnalpha3*rhof
;

drhofPLUSrhofdLnalphaPsi6uz1
=
drhof1 + dLnalphaPsi6uz1*rhof
;

drhofPLUSrhofdLnalphaPsi6uz2
=
drhof2 + dLnalphaPsi6uz2*rhof
;

drhofPLUSrhofdLnalphaPsi6uz3
=
drhof3 + dLnalphaPsi6uz3*rhof
;

lLnuzero
=
luzero/uzero
;

dLnuzero1
=
(0.5*duzerosqr1)/uzerosqr
;

dLnuzero2
=
(0.5*duzerosqr2)/uzerosqr
;

dLnuzero3
=
(0.5*duzerosqr3)/uzerosqr
;

lduzero1
=
(0.5*lduzerosqr1 - (0.25*duzerosqr1*luzerosqr)/uzerosqr)/uzero
;

lduzero2
=
(0.5*lduzerosqr2 - (0.25*duzerosqr2*luzerosqr)/uzerosqr)/uzero
;

lduzero3
=
(0.5*lduzerosqr3 - (0.25*duzerosqr3*luzerosqr)/uzerosqr)/uzero
;

ldLnuzero1
=
-(dLnuzero1*lLnuzero) + lduzero1/uzero
;

ldLnuzero2
=
-(dLnuzero2*lLnuzero) + lduzero2/uzero
;

ldLnuzero3
=
-(dLnuzero3*lLnuzero) + lduzero3/uzero
;

lhfuzeroPsi4
=
hf*(luzero*Psi4 + 4.*Psi3*uzero*lPsi[ijk])
;

ldLnalphaPsi21
=
ldLnalphaP1 + ldLnPsi1
;

ldLnalphaPsi22
=
ldLnalphaP2 + ldLnPsi2
;

ldLnalphaPsi23
=
ldLnalphaP3 + ldLnPsi3
;

ldLnalpha1
=
ldLnalphaP1 - ldLnPsi1
;

ldLnalpha2
=
ldLnalphaP2 - ldLnPsi2
;

ldLnalpha3
=
ldLnalphaP3 - ldLnPsi3
;

ldLnalphaPsi6uz1
=
ldLnalphaP1 + 5.*ldLnPsi1 + ldLnuzero1
;

ldLnalphaPsi6uz2
=
ldLnalphaP2 + 5.*ldLnPsi2 + ldLnuzero2
;

ldLnalphaPsi6uz3
=
ldLnalphaP3 + 5.*ldLnPsi3 + ldLnuzero3
;

ldrhofPLUSrhofdLnalphaPsi21
=
ldLnalphaPsi21*rhof
;

ldrhofPLUSrhofdLnalphaPsi22
=
ldLnalphaPsi22*rhof
;

ldrhofPLUSrhofdLnalphaPsi23
=
ldLnalphaPsi23*rhof
;

ldrhofPLUSrhofdLnalpha1
=
ldLnalpha1*rhof
;

ldrhofPLUSrhofdLnalpha2
=
ldLnalpha2*rhof
;

ldrhofPLUSrhofdLnalpha3
=
ldLnalpha3*rhof
;

ldrhofPLUSrhofdLnalphaPsi6uz1
=
ldLnalphaPsi6uz1*rhof
;

ldrhofPLUSrhofdLnalphaPsi6uz2
=
ldLnalphaPsi6uz2*rhof
;

ldrhofPLUSrhofdLnalphaPsi6uz3
=
ldLnalphaPsi6uz3*rhof
;

FlSigma[ijk]
=
dlSigmaUp1*drhofPLUSrhofdLnalphaPsi21 + 
  dlSigmaUp2*drhofPLUSrhofdLnalphaPsi22 + 
  dlSigmaUp3*drhofPLUSrhofdLnalphaPsi23 + 
  dSigmaUp1*ldrhofPLUSrhofdLnalphaPsi21 + 
  dSigmaUp2*ldrhofPLUSrhofdLnalphaPsi22 + 
  dSigmaUp3*ldrhofPLUSrhofdLnalphaPsi23 - 
  rhof*(divbeta*lhfuzeroPsi4 + divlbeta*hf*Psi4*uzero) - 
  (drhofPLUSrhofdLnalphaPsi6uz1*lhfuzeroPsi4 + 
     hf*ldrhofPLUSrhofdLnalphaPsi6uz1*Psi4*uzero)*beta1[ijk] - 
  (drhofPLUSrhofdLnalphaPsi6uz2*lhfuzeroPsi4 + 
     hf*ldrhofPLUSrhofdLnalphaPsi6uz2*Psi4*uzero)*beta2[ijk] - 
  (drhofPLUSrhofdLnalphaPsi6uz3*lhfuzeroPsi4 + 
     hf*ldrhofPLUSrhofdLnalphaPsi6uz3*Psi4*uzero)*beta3[ijk] - 
  hf*Psi4*uzero*(drhofPLUSrhofdLnalphaPsi6uz1*lB1[ijk] + 
     drhofPLUSrhofdLnalphaPsi6uz2*lB2[ijk] + 
     drhofPLUSrhofdLnalphaPsi6uz3*lB3[ijk]) + 
  rhof*(ddlSigma11[ijk] + ddlSigma22[ijk] + ddlSigma33[ijk] - 
     2.*divwB*Psim3*lPsi[ijk]) - 
  2.*Psim3*lPsi[ijk]*(drhofPLUSrhofdLnalpha1*wB1[ijk] + 
     drhofPLUSrhofdLnalpha2*wB2[ijk] + drhofPLUSrhofdLnalpha3*wB3[ijk]) + 
  Psim2*(drhofPLUSrhofdLnalpha1*lwB1 + drhofPLUSrhofdLnalpha2*lwB2 + 
     drhofPLUSrhofdLnalpha3*lwB3 + divlwB*rhof + 
     ldrhofPLUSrhofdLnalpha1*wB1[ijk] + ldrhofPLUSrhofdLnalpha2*wB2[ijk] + 
     ldrhofPLUSrhofdLnalpha3*wB3[ijk])
;


} else if(LaplaceSigmaOutside) { 

FlSigma[ijk]
=
ddlSigma11[ijk] + ddlSigma22[ijk] + ddlSigma33[ijk]
;


} else { 

FlSigma[ijk]
=
dddlSigmadlam3[ijk] + 2.*ddlSigmadlam2[ijk] + dlSigmadlam[ijk]
;


} /* end Outside cases */ 


} else { 

FlSigma[ijk]
=
lSigma[ijk]
;


} 

}
/* if (FakeT0) */


FlPsi[ijk]
=
CoordFac[ijk]*FlPsi[ijk]
;

FlB1[ijk]
=
CoordFac[ijk]*FlB1[ijk]
;

FlB2[ijk]
=
CoordFac[ijk]*FlB2[ijk]
;

FlB3[ijk]
=
CoordFac[ijk]*FlB3[ijk]
;

FlalphaP[ijk]
=
CoordFac[ijk]*FlalphaP[ijk]
;

FlSigma[ijk]
=
CoordFac[ijk]*FlSigma[ijk]
;

}
/* if (FakeT0) */



} /* end of points loop */ 

} /* end of boxes */


}  /* end of function */

/* DNS_CTS.c */
/* nvars = 181, n* = 1454,  n/ = 219,  n+ = 1188, n = 2861, O = 1 */
