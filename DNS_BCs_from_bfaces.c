
////////////////////////////////////////////////////////////////////////////
// This funcs in this file should later go into src/utility/boundary/
////////////////////////////////////////////////////////////////////////////

/* BCs_from_bfaces.c */
/* Wolfgang Tichy 2/2018 */

#include "sgrid.h"



// the next 2 are for Psi(box) = Psi(obox)

// we will also need n*dPsi(box) = n*dPsi(obox)
// ??? so much work ... :(


/* interpolate onto points in box on bface fi 
   this uses 1d interplolation is done in other box, in direction idir
   FPsi = Psi(box) - Psi_interp(obox) */
void FPsi_1Dinterp_for_bface(int iFPsi, int iPsi, tBox *box, int fi, int idir)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  tBface *bface = box->bface[fi];
  int dir = 1 + bface->f/2;
  int ob  = bface->ob;
  int ofi = bface->ofi;
  int oXi = bface->oXi;
  int oYi = bface->oYi;
  int oZi = bface->oZi;
  tBox *obox = NULL;
  int oCi;
  tBface *obface;
  int of, odir, op, plN;
  double *P, *Pcoeffs, Pinterp;
  int pi, ind;

  /* check if there is another box and another bface */
  if(ob<0)  errorexit("can only interpolate if there is another box");
  if(ofi<0) errorexit("can only interpolate if there is another bface");

  /* set other box props */
  obox = grid->box[ob];
  obface = obox->bface[ofi];
  of = obface->f;
  odir = 1+of/2;

  P = obox->v[iPsi];              /* values in box ob */
  Pcoeffs = dmalloc(obox->nnodes); /* get mem for coeffs of P */

  /* if interp-dir is in oither face plane */
  if(odir!=idir)
  {
    plN = odir;
    if(odir==1)      op = ( (obox->n1-1) )*(of%2);
    else if(odir==2) op = ( (obox->n2-1) )*(of%2);
    else             op = ( (obox->n3-1) )*(of%2);
  }
  else errorexit("implement odir==idir case");

  /* set oCi to oXi, oYi, or oZi */
  if(idir==1)      oCi = oXi;
  else if(idir==2) oCi = oYi;
  else if(idir==3) oCi = oZi;
  else errorexit("1<=idir<=3 is required");
  if(oCi<=0) errorexit("oCi <= 0");

  /* gwt coeffs in dir idir and plane plN,op */
  spec_analysis1_inplaneN(obox, idir, plN, op, P, Pcoeffs);

  forPointList_inbox(bface->fpts, box, pi, ind)
  {
    int i1, i2;
    int ok = kOfInd_n1n2(ind, n1,n2);
    int oj = jOfInd_n1n2_k(ind, n1,n2, ok);
    int oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);
    double oC = box->v[oCi][ind]; /* this can be oX, oY, or oZ */

    errorexit("oi, oj, ok are only correct if coords in box and obox are "
              "aligned and the number of points agree n1=on1, ... ");

    if(odir==1)      oi = op;
    else if(odir==2) oj = op;
    else             ok = op;

    if(idir==1)      { i1 = oj; i2 = ok; }
    else if(idir==2) { i1 = oi; i2 = ok; }
    else             { i1 = oi; i2 = oj; }

    Pinterp = spec_interpolate_in_dir_at_i1_i2(obox, idir, 
                                               i1,i2, Pcoeffs, oC);
    FPsi[ind] = Psi[ind] - Pinterp;
    if(!finite(Pinterp))
    {
      printf("Pinterp=%g  i1=%d oC=%.13g i2=%d  ind=%d\n",
             Pinterp, i1,oC,i2, ind);
      printbface(bface);
      printbface(obface);
      grid->time  = 42;
      write_grid(grid);
      errorexit("Pinterp is not finite!");
    }
  }
  free(Pcoeffs);
}

/* interpolate onto points in box on bface fi 
   this uses 2d interplolation is done in other box, in plane plN
   FPsi = Psi(box) - Psi_interp(obox) */
void FPsi_2Dinterp_for_bface(int iFPsi, int iPsi, tBox *box, int fi, int plN)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  tBface *bface = box->bface[fi];
  int dir = 1 + bface->f/2;
  int ob  = bface->ob;
  int ofi = bface->ofi;
  int oXi = bface->oXi;
  int oYi = bface->oYi;
  int oZi = bface->oZi;
  tBox *obox = NULL;
  int oCi1, oCi2;
  tBface *obface;
  int of, op;
  double *P, *Pcoeffs, Pinterp;
  int pi, ind;

  /* check if there is another box and another bface */
  if(ob<0)  errorexit("can only interpolate if there is another box");
  if(ofi<0) errorexit("can only interpolate if there is another bface");

  /* set other box props */
  obox = grid->box[ob];
  obface = obox->bface[ofi];
  of = obface->f;
  //odir = 1+of/2;

  P = obox->v[iPsi];              /* values in box ob */
  Pcoeffs = dmalloc(obox->nnodes); /* get mem for coeffs of P */

  /* find plane index and set oCi1/2 to oXi, oYi, or oZi */
  if(plN==1)
  {
    op = ( (obox->n1-1) )*(of%2);
    oCi1 = oYi;
    oCi2 = oZi;
  }
  else if(plN==2)
  {
    op = ( (obox->n2-1) )*(of%2);
    oCi1 = oXi;
    oCi2 = oZi;
  }
  else if(plN==3)
  {
    op = ( (obox->n3-1) )*(of%2);
    oCi1 = oXi;
    oCi2 = oYi;
  }
  else errorexit("1<=plN<=3 is required");
  if(oCi1<=0 || oCi2<=0) errorexit("oCi1/2 > 0 is required");

  /* get values by 2d interpolation from box ob */
  spec_Coeffs_inplaneN(obox, plN,op, P, Pcoeffs);

  forPointList_inbox(bface->fpts, box, pi, ind)
  {
    double X1 = box->v[oCi1][ind];
    double X2 = box->v[oCi2][ind];

    Pinterp = spec_interpolate_inplaneN(obox, plN, op,
                                        Pcoeffs, X1,X2); 
    FPsi[ind] = Psi[ind] - Pinterp;
    if(!finite(Pinterp))
    {
      printf("Pinterp=%g  X1=%.13g X2=%.13g  ind=%d\n",
              Pinterp, X1,X2, ind);
      printbface(bface);
      printbface(obface);
      grid->time  = 42;
      write_grid(grid);
      errorexit("Pinterp is not finite!");
    }
  }
  free(Pcoeffs);
}


/* interpolate onto points in box on bface fi 
   this uses 3d interplolation is done in other box
   FPsi = Psi(box) - Psi_interp(obox) */
void FPsi_3Dinterp_for_bface(int iFPsi, int iPsi, tBox *box, int fi)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  tBface *bface = box->bface[fi];
  tBox *obox = NULL;
  double *P, *Pcoeffs, Pinterp;
  int ob  = bface->ob;
  int oXi = bface->oXi;
  int oYi = bface->oYi;
  int oZi = bface->oZi;
  int pi, ind;

  /* check if there is another box */
  if(ob<0)  errorexit("can only interpolate if there is another box");
  if(oXi<=0 || oYi <=0 || oZi<=0)
    errorexit("oXi, oYi, oZi are not all > 0");

  /* set other box props */
  obox = grid->box[ob];
  P = obox->v[iPsi];              /* values in box ob */
  Pcoeffs = dmalloc(obox->nnodes); /* get mem for coeffs of P */

    /* get values by 3d interpolation from box ob */
  spec_Coeffs(obox, P, Pcoeffs);
  forPointList_inbox(bface->fpts, box, pi, ind)
  {
    double X = box->v[oXi][ind];
    double Y = box->v[oYi][ind];
    double Z = box->v[oZi][ind];
    Pinterp = spec_interpolate(obox, Pcoeffs, X,Y,Z);
    if(!finite(Pinterp))
    {
      printf("Pinterp=%g  X=%.13g Y=%.13g Z=%.13g  ind=%d\n",
             Pinterp, X,Y,Z, ind);
      NumberChecker_CheckIfFinite(grid, "DNSdata_temp1");
      printbox(obox);
      grid->time  = 42;
      write_grid(grid);
      errorexit("Pinterp is not finite!");
    }
    FPsi[ind] = Psi[ind] - Pinterp;
  }
  free(Pcoeffs);
}
