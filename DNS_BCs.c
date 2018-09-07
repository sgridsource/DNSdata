/* set BCs */
/* Wolfgang Tichy 2010 */


#include "sgrid.h"
#include "DNSdata.h"

#define Power pow


/* functions in this file */
void set_Sigma_C0_in1ATlamA0B0_BC(tBox *box, int iFSigma, int iSigma,
                                  int domOfStar1, double lam);



/* set BC's between boxes and at outerbound */
void DNS_set_interbox_and_outer_BCs(tBox *box, int iFPsi, int iPsi,
                                    int iPsix, int iPsiy, int iPsiz,
                                    double PsiFarLimit, int setFarLimit,
                                    intList *skip_f)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  int idPsi[4];
  int fi;
  idPsi[1] = iPsix;
  idPsi[2] = iPsiy;
  idPsi[3] = iPsiz;

  /* loop over bfaces */
  forallbfaces(box, fi)
  {
    tBface *bface = box->bface[fi];
    int ob  = bface->ob;
    int pi, ind;

    /* do nothing if bface->f is in skip_f */
    if(in_intList(skip_f, bface->f)) continue;

    /* check if there is other box */
    if(ob>=0)
    {
      /* set BCs for cases where there is another box */
      set_interbox_BCs_for_bface(iFPsi, bface, iPsi, idPsi);
    }
    else  /* there is no box */
    {
      /* set far limit BC */
      if(bface->outerbound && setFarLimit)
        forPointList_inbox(bface->fpts, box, pi, ind)
          FPsi[ind] = Psi[ind] - PsiFarLimit;
    }

  } /* end of of forallbfaces */
}



/* standard BCs for all fields */
void set_DNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs,
                     int nonlin)
{
  int FakeMatterOutside = Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside");
  tGrid *grid = vlu->grid;
  int vind;
  int vindDerivs=0;
  tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vlu->vlPars;

  for(vind=0; vind<vlu->n; vind++)
  {
    int b;
    int iFPsi = vlFu->index[vind];
    int iPsi  = vlu->index[vind];
    int iPsix = vluDerivs->index[vindDerivs];
    int iPsiy = vluDerivs->index[vindDerivs+1];
    int iPsiz = vluDerivs->index[vindDerivs+2];
    int ncomp = VarNComponents(iPsi);
    double PsiFarLimit = VarFarLimit(iPsi)*nonlin;
    char *varname = VarName(vlu->index[vind]);
    int is_Sigma = 0;
    intList *skip_f;

    /* do nothing and goto end of loop if var with vind is not the one
       of the current block */
    if(blkinfo!=NULL) if(vlu->index[vind] != blkinfo->vari)
                        goto Incr_vindDerivs;

    if(strstr(varname, "DNSdata_Sigma"))
    {
      is_Sigma = 1;
      /* printf("varname=%s\n", varname); */
    }

    /* box loop */
    skip_f = alloc_intList(); /* list of face we want to skip */
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[iFPsi];
      double *Psi  = box->v[iPsi];
      double *Psix = box->v[iPsix];
      double *Psiy = box->v[iPsiy];
      double *Psiz = box->v[iPsiz];

      /* do nothing and continue if current block is not in box b */
      if(blkinfo!=NULL) if(b!=blkinfo->bi) continue;

      /* list of faces we want to skip */
      clear_intList(skip_f);

      /* DNSdata_Sigma? */
      if(is_Sigma)
      {
        /* do nothing else for DNSdata_Sigma away from stars */
        if(box->MATTR == AWAY) continue;

        /* we need BCs for Sigma in touching box at lam=0, but not elsewhere */
        if(box->MATTR == TOUCH)
        {
          int f;
          /* make a list of faces to omit in DNS_set_interbox_and_outer_BCs */
          push_intList(skip_f, 1);
          //if(!FakeMatterOutside) for(f=2; f<6; f++) push_intList(skip_f, f);

          /* set normal deriv equal at lam=0 of TOUCH box */
          if(0 && FakeMatterOutside)
          {
            int idPsi[4];
            idPsi[1] = iPsix;
            idPsi[2] = iPsiy;
            idPsi[3] = iPsiz;
            push_intList(skip_f, 0);
            /* set normal derivs the same for face0 */
            set_interbox_BC_onface(box, iFPsi, 0, iPsi, idPsi, 1);
            ///*  if we are in box2 or box16 set C0 in Psi=Sigma */
            //set_Sigma_C0_in1ATlamA0B0_BC(box, iFPsi, iPsi, 1, 0.0);
          }
        }
        /* we may need a special BC at lam=1 of the inside boxes */
        if( (box->MATTR == INSIDE) && (box->BOUND == SSURF) )
        {
          /* set Sigma at lam=1 of INSIDE box at one point */
          if(0 && FakeMatterOutside)
          {
            /* make list of faces to omit in DNS_set_interbox_and_outer_BCs */
            push_intList(skip_f, 1);

            /* if we are in box1 or box15 set C0 in Psi=Sigma */
            set_Sigma_C0_in1ATlamA0B0_BC(box, iFPsi, iPsi, 0, 1.0);
          }
        }
      }

      /* set some BCs for each box */
      DNS_set_interbox_and_outer_BCs(box, iFPsi, iPsi, iPsix,iPsiy,iPsiz,
                                     PsiFarLimit, 1, skip_f);
    } /* end forallboxes */
    free_intList(skip_f);

    Incr_vindDerivs:
      /* increase index for derivs */
      vindDerivs += 3;
      if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}


/* impose: DNSdata_Sigma[i] = Omega*(xc1-xCM) * y
   or:     DNSdata_Sigma[i] = Omega*(xc2-xCM) * y
   on outside (lam=1) of TOUCH boxes */
void set_Sigma_Omega_r_y_BC(tVarList *vlFu, tVarList *vlu,
                            tVarList *vluDerivs, int nonlin)
{
  double Omega = Getd("DNSdata_Omega");
  double xCM = Getd("DNSdata_x_CM");
  double xmax1 = Getd("DNSdata_actual_xmax1");
  double xmax2 = Getd("DNSdata_actual_xmax2");
  tGrid *grid = vlu->grid;
  int vind;
  int vindDerivs=0;
  tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vlu->vlPars;

  for(vind=0; vind<vlu->n; vind++)
  {
    int b;
    int ix = Ind("x");
    int iX = Ind("X");
    int iFSigma = vlFu->index[vind];
    int iSigma  = vlu->index[vind];
    int ncomp = VarNComponents(iSigma);
    char *varname = VarName(vlu->index[vind]);

    /* do nothing and goto end of loop if var with vind is not the one
       of the current block */
    if(blkinfo!=NULL) if(vlu->index[vind] != blkinfo->vari)
                        goto Incr_vindDerivs2;

    /* do nothing if var is not Sigma */
    if(!strstr(varname, "DNSdata_Sigma")) goto Incr_vindDerivs2;

    /* box loop */
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int i,j,k;
      double *FSigma = box->v[iFSigma];
      double *Sigma  = box->v[iSigma];
      double *y = box->v[ix+1];
      double Omega_r, xc;

      /* do nothing and continue if current block is not in box b */
      if(blkinfo!=NULL) if(b!=blkinfo->bi) continue;

      /* do nothing else for DNSdata_Sigma inside and away from stars */
      if(box->MATTR == AWAY || box->MATTR == INSIDE) continue;

      if(box->x_of_X[1]==NULL) /* if Cartesian box */
        y = box->v[iX+1];

      /* figure out xc and Omega_r */
      if (box->SIDE==STAR1) xc = xmax1;
      else                  xc = xmax2;
      Omega_r = Omega*(xc-xCM);

      /* we need BCs for Sigma in touching box at lam=1, but not elsewhere */
      forplane1(i,j,k, n1,n2,n3, n1-1)
      {
        int ijk=Index(i,j,k);
        /* FSigma[ijk] = Sigma[ijk] - (Omega*(xc1-xCM) * y[ijk]) * nonlin; */
        FSigma[ijk] = Sigma[ijk] - (Omega_r * y[ijk]) * nonlin;
      }
    } /* end forallboxes */

    Incr_vindDerivs2:
      /* increase index for derivs */
      vindDerivs += 3;
      if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}

/* impose: DNSdata_Sigma[i] = 0 for each star
   at ONE point on outside (lam=1,A=B=0) of one TOUCH box,
   set this BC at the point with index (i,j,k) = (n1-1, n2/2, n3/2) */
void set_Sigma_0_in1TOUCHATlam1A0B0_BC(tVarList *vlFu, tVarList *vlu,
                                       tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int vind;
  int vindDerivs=0;
  tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vlu->vlPars;

  for(vind=0; vind<vlu->n; vind++)
  {
    int b;
    int iFSigma = vlFu->index[vind];
    int iSigma  = vlu->index[vind];
    int ncomp = VarNComponents(iSigma);
    char *varname = VarName(vlu->index[vind]);

    /* do nothing and goto end of loop if var with vind is not the one
       of the current block */
    if(blkinfo!=NULL) if(vlu->index[vind] != blkinfo->vari)
                        goto Incr_vindDerivs2;

    /* do nothing if var is not Sigma */
    if(!strstr(varname, "DNSdata_Sigma")) goto Incr_vindDerivs2;

    /* box loop */
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int ijk;
      double *FSigma = box->v[iFSigma];
      double *Sigma  = box->v[iSigma];
      double *co;
      double Sig_100;
      int isXinDom    = (box->CI->dom == box->SIDE - STAR1);

      /* do nothing and continue if current block is not in box b */
      if(blkinfo!=NULL) if(b!=blkinfo->bi) continue;

      /* do nothing else for DNSdata_Sigma inside and away from stars */
      if(box->MATTR == AWAY || box->MATTR == INSIDE) continue;

      /* do nothing if we are not in dom0 for STAR1 or dom1 for STAR2 */
      if(!isXinDom) continue;

      /* find Sigma at (lam,A,B)=(1,0,0) by 2D interpolation */
      co = dmalloc(box->nnodes); /* get mem for coeffs in co */
      spec_Coeffs_inplaneN(box, 1,n1-1, Sigma, co);
      Sig_100 = spec_interpolate_inplaneN(box, 1,n1-1, co, 0.,0.);
      free(co);

      /* we need BCs for Sigma in touching box at lam=1, A=B=0 */
      ijk = Index(n1-1, n2/2, n3/2);
      FSigma[ijk] = Sig_100; // - 0. * nonlin;
    } /* end forallboxes */

    Incr_vindDerivs2:
      /* increase index for derivs */
      vindDerivs += 3;
      if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}


/* impose: DNSdata_Sigma[i]_INSIDE = DNSdata_Sigma[i]_TOUCH for each star
   at ONE point at (lam=1/0,A=B=0) of one INSIDE/TOUCH box,
   set this BC at the point with index (i,j,k) = ((n1-1)*lam, n2/2, n3/2) */
void set_Sigma_C0_in1ATlamA0B0_BC(tBox *box, int iFSigma, int iSigma,
                                  int domOfStar1, double lam)
{
  tGrid *grid = box->grid;
  tBox *obox;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int b = box->b;
  int ijk, i, oi;
  double *FSigma = box->v[iFSigma];
  double *Sigma  = box->v[iSigma];
  double *co, *oSigma;
  double Sig_100, oSig_000;
  int BCdom;

  /* find dom in which we apply BC */
  if(box->SIDE == STAR1)
    BCdom = domOfStar1;
  else
    switch(domOfStar1)
    {
    case 0:
    case 2:
      BCdom = domOfStar1 + 1;
      break;
    case 1:
    case 3:
      BCdom = domOfStar1 - 1;
      break;
    default:
      BCdom = domOfStar1;
    }

  /* do nothing if we are not in dom0 for STAR1 or dom1 for STAR2, ... */
  if(box->CI->dom != BCdom) return;

  /* get mem for coeffs in co */
  co = dmalloc(box->nnodes);

  /* find Sigma at (lam,A,B)=(1,0,0) by 2D interpolation */
  i = (n1-1)*lam;
  spec_Coeffs_inplaneN(box, 1,i, Sigma, co);
  Sig_100 = spec_interpolate_inplaneN(box, 1,i, co, 0.,0.);

  /* find oSigma at (lam,A,B)=(0,0,0) by 2D interpolation */
  obox = grid->box[b+6]; // for our setup the TOUCH box is 6 boxes further
  oi = (obox->n1-1)*(1.0-lam);
  oSigma = obox->v[iSigma];
  spec_Coeffs_inplaneN(obox, 1,oi, oSigma, co);
  oSig_000 = spec_interpolate_inplaneN(obox, 1,oi, co, 0.,0.);

  free(co);

  /* we need BCs for Sigma in touching box at lam=1, A=B=0 */
  ijk = Index(i, n2/2, n3/2);
  FSigma[ijk] = Sig_100 - oSig_000; // - 0. * nonlin;
}
