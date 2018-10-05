/* DNS_compute_chi.c */
/* dtim 2/2014, WT 10/2018 */

/* WT: we compute this:

   chieq   = dh/dx = dq/dx evaluated on equator at point lam=1, A=B=0
   chipole = dh/dz = dq/dz evaluated at pole at point lam=1, A=B=0

   ==> chi = chieq/chipole = (dq/dx)/(dq/dz)

   Note: h = q + 1,
   d\ln h = dh/h = dq/h
   On star surface h=1 => d\ln h = dq

   The mass shedding indicator is defined in Eq. (104) of
   https://journals.aps.org/prd/pdf/10.1103/PhysRevD.63.064029
   where H = \ln h = \ln(q+1).
   
   chi goes to zero if mass shedding happens */

#include "sgrid.h"
#include "DNSdata.h"


int DNS_compute_chi(tGrid *grid)
{
  int bi;
  int index_DNSdata_q     = Ind("DNSdata_q");
  int index_DNSdata_qx    = Ind("DNSdata_qx");
  int index_DNSdata_temp1 = Ind("DNSdata_temp1");
  double chieq1=1e20, chipole1=1e10, chieq2=1e20, chipole2=1e10;
  double chi1, chi2;

  /* do nothing if BNSdata_Interpolate_pointsfile exists */
  if(GetsLax("BNSdata_Interpolate_pointsfile")!=0) return 0;

  /* loop over boxes that have star parts we need */
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
    double *qx = box->v[index_DNSdata_qx];
    double *qz = box->v[index_DNSdata_qx+2];
    double *c  = box->v[index_DNSdata_temp1];
    int hasSSURF    = (box->BOUND == SSURF);
    int MATTRinside = (box->MATTR == INSIDE);
    int isxinDom    = (box->CI->dom == box->SIDE - STAR1); // has xmin1/2
    int iszmaxDom   = (box->CI->dom == 5); // has northpole of star
    int star        = (box->SIDE); 

    if(MATTRinside && hasSSURF)
    {
      double chieq, chipole;
      if(isxinDom)
      {
        FirstDerivsOf_S(box, index_DNSdata_q, index_DNSdata_qx);
        spec_Coeffs(box, qx, c);
        chieq = spec_interpolate(box, c, 1.,0.,0.); /* interp. to lam=1 */
        if(star==STAR1) chieq1 = -chieq;
        if(star==STAR2) chieq2 = +chieq;
      }
      if(iszmaxDom)
      {
        FirstDerivsOf_S(box, index_DNSdata_q, index_DNSdata_qx);
        spec_Coeffs(box, qz, c);
        chipole = spec_interpolate(box, c, 1.,0.,0.); /* interp. to lam=1 */
        if(star==STAR1) chipole1 = chipole;
        if(star==STAR2) chipole2 = chipole;
      }
    }
  }

  /* compute and set chi */
  chi1 = chieq1/chipole1;
  chi2 = chieq2/chipole2;
  Setd("DNSdata_mass_shedding1", chi1);
  Setd("DNSdata_mass_shedding2", chi2);

  return 0;
}
