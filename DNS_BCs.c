/* set BCs */
/* Wolfgang Tichy 2010 */


#include "sgrid.h"
#include "DNSdata.h"

#define Power pow


/* functions in this file */
/* ... */




/* find approx normal vector of box faces */
void DNS_approx_normal(tBface *bface, int ijk, double *nx, double *ny, double *nz)
{
  tGrid *grid = bface->grid;
  int b = bface->b;
  //int f = bface->f;
  tBox *box = grid->box[b];

  if(box->COORD == CART)
    errorexit("implement normal for Cartesian boxes");
  else /* assume AnsorgNS */
  {
    double s = 1.0 - 2.0*(box->SIDE == STAR2);
    double A = box->v[Ind("X")][ijk];

    if(dequal(A,1.0))
    {
      *nx = s;
      *ny = *nz = 0.0;
    }
    else
    {
      double B   = box->v[Ind("Y")][ijk];
      double phi = box->v[Ind("Z")][ijk];
      /* We use an approximate normal vec */
      *nx = s*cos(PI*B);
      *ny = sin(PI*B)*cos(phi);
      *nz = sin(PI*B)*sin(phi);
    }
  }
}


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


/* set BC's between boxes and at outerbound */
void set_interbox_and_FarLimit_BCs(tBox *box, int iFPsi, int iPsi,
                                   int iPsix, int iPsiy, int iPsiz,
                                   double PsiFarLimit, int setFarLimit,
                                   intList *skip_f)
{
  tGrid *grid = box->grid;
  double *FPsi = box->v[iFPsi];
  double *Psi  = box->v[iPsi];
  double *Psix = box->v[iPsix];
  double *Psiy = box->v[iPsiy];
  double *Psiz = box->v[iPsiz];
  int fi;

  /* loop over bfaces */
  for(fi=0; fi<box->nbfaces; fi++)
  {
    tBface *bface = box->bface[fi];
    int ob  = bface->ob;
    int ofi = bface->ofi;
    tBox *obox = NULL;
    int oXi = bface->oXi;
    int oYi = bface->oYi;
    int oZi = bface->oZi;
    int pi, ind;

    /* do nothing if bface->f is in skip_f */
    if(in_intList(skip_f, bface->f)) continue;

    /* check if there is another box */
    if(ob>=0)
    {
      double *P;
      double *dP[4];

      obox = grid->box[ob];

      if(bface->touch) /* bface touches one other face */
      {
        if(bface->sameX && bface->sameY && bface->sameZ)
        {
          if(bface->setnormalderiv) /* field derivs are equal */
          {
            dP[1] = obox->v[iPsix]; /* derivs in other box */
            dP[2] = obox->v[iPsiy];
            dP[3] = obox->v[iPsiz];
            forPointList_inbox(bface->fpts, box, pi, ind)
            {
              double nx, ny, nz;
              DNS_approx_normal(bface, ind, &nx, &ny, &nz);
              FPsi[ind] = nx * (Psix[ind] - dP[1][ind]) +
                          ny * (Psiy[ind] - dP[2][ind]) +
                          nz * (Psiz[ind] - dP[3][ind]);
            }
          }
          else /* fields are equal */
          {
            P = obox->v[iPsi]; /* values in box ob */
            forPointList_inbox(bface->fpts, box, pi, ind)
              FPsi[ind] = Psi[ind] - P[ind];
          }
        }
        else if(bface->sameX && bface->sameZ)
        {
          int n1 = box->n1;
          int n2 = box->n2;
          int n3 = box->n3;
          int dir = 1 + bface->f/2;
          tBface *obface = obox->bface[ofi];
          int of = obface->f;
          int odir = 1+of/2;
          int op;

          if(odir==1)      op = ( (obox->n1-1) )*(of%2);
          else if(odir==2) op = ( (obox->n2-1) )*(of%2);
          else if(odir==3) op = ( (obox->n3-1) )*(of%2);

          if(dir==1) /* X-dir => interpolate in Y-dir */
          {
            if(bface->setnormalderiv) /* field derivs are equal */
            {
              /* get deriv values by 1d interpolation from box ob */
              double dPinterp[4];
              double *dcoeffs[4];
              dcoeffs[1] = obox->v[Ind("DNSdata_temp1")];
              dcoeffs[2] = obox->v[Ind("DNSdata_temp2")];
              dcoeffs[3] = obox->v[Ind("DNSdata_temp3")];
              dP[1] = obox->v[iPsix]; /* derivs in other box */
              dP[2] = obox->v[iPsiy];
              dP[3] = obox->v[iPsiz];
              spec_analysis1_inplaneN(obox, 2, odir, op, dP[1], dcoeffs[1]);
              spec_analysis1_inplaneN(obox, 2, odir, op, dP[2], dcoeffs[2]);
              spec_analysis1_inplaneN(obox, 2, odir, op, dP[3], dcoeffs[3]);
////spec_Coeffs(obox, dP[1], dcoeffs[1]);
////spec_Coeffs(obox, dP[2], dcoeffs[2]);
////spec_Coeffs(obox, dP[3], dcoeffs[3]);

              forPointList_inbox(bface->fpts, box, pi, ind)
              {
                int ok = kOfInd_n1n2(ind, n1,n2);
                int oj = jOfInd_n1n2_k(ind, n1,n2, ok);
                int oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);
                //double X = box->v[oXi][ind];
                double Y = box->v[oYi][ind];
                //double Z = box->v[oZi][ind];
                double nx, ny, nz;

                if(odir==1)      oi = op;
                else if(odir==2) oj = op;
                else if(odir==3) ok = op;

                dPinterp[1] = spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok,
                                                               dcoeffs[1], Y);
                dPinterp[2] = spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok,
                                                               dcoeffs[2], Y);
                dPinterp[3] = spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok,
                                                               dcoeffs[3], Y);
//double X = box->v[oXi][ind];
//double Z = box->v[oZi][ind];
////dPinterp[1] = spec_interpolate(obox, dcoeffs[1], X,Y,Z);
////dPinterp[2] = spec_interpolate(obox, dcoeffs[2], X,Y,Z);
////dPinterp[3] = spec_interpolate(obox, dcoeffs[3], X,Y,Z);
//printf("b%d ob%d %d %d %d X=%g Y=%g Z=%g : %g %g %g %s\n",
//box->b,ob,oi,oj,ok, X,Y,Z, dPinterp[1],dPinterp[2],dPinterp[3], VarName(iPsi));

                DNS_approx_normal(bface, ind, &nx, &ny, &nz);
                FPsi[ind] = nx * (Psix[ind] - dPinterp[1]) +
                            ny * (Psiy[ind] - dPinterp[2]) +
                            nz * (Psiz[ind] - dPinterp[3]);
                if(!finite(dPinterp[1]) ||
                   !finite(dPinterp[2]) || !finite(dPinterp[3]))
                {
                  NumberChecker_CheckIfFinite(grid, "DNSdata_temp1");
                  NumberChecker_CheckIfFinite(grid, "DNSdata_temp2");
                  NumberChecker_CheckIfFinite(grid, "DNSdata_temp3");
                  printf("dPinterp[1]=%g  oi=%d Y=%.13g ok=%d  ind=%d\n",
                         dPinterp[1], oi,Y,ok, ind);
                  printf("dPinterp[2]=%g\n", dPinterp[2]);
                  printf("dPinterp[3]=%g\n", dPinterp[3]);
                  printbface(bface);
                  printbface(obface);
                  grid->time  = 42;
                  write_grid(grid);
                  errorexit("dPinterp[] is not finite!");
                }
              }
            }
            else /* fields are equal */
            {
              /* get values by 1d interpolation from box ob */
              double Pinterp;
              double *Pcoeffs = obox->v[Ind("DNSdata_temp1")];
              P = obox->v[iPsi]; /* values in box ob */
              spec_analysis1_inplaneN(obox, 2, odir, op, P, Pcoeffs);

////spec_Coeffs(obox, P, Pcoeffs);

              forPointList_inbox(bface->fpts, box, pi, ind)
              {
                int ok = kOfInd_n1n2(ind, n1,n2);
                int oj = jOfInd_n1n2_k(ind, n1,n2, ok);
                int oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);
                //double X = box->v[oXi][ind];
                double Y = box->v[oYi][ind];
                //double Z = box->v[oZi][ind];

                if(odir==1)      oi = op;
                else if(odir==2) oj = op;
                else if(odir==3) ok = op;

                Pinterp =
                  spec_interpolate_in_dir_at_i1_i2(obox, 2, oi,ok, Pcoeffs, Y);

//double X = box->v[oXi][ind];
//double Z = box->v[oZi][ind];
////Pinterp = spec_interpolate(obox, Pcoeffs, X,Y,Z);
//printf("b%d ob%d %d %d %d X=%g Y=%g Z=%g : %g %s\n",
//box->b,ob, oi,oj,ok, X,Y,Z, Pinterp, VarName(iPsi));

                FPsi[ind] = Psi[ind] - Pinterp;
                if(!finite(Pinterp))
                {
                  NumberChecker_CheckIfFinite(grid, "DNSdata_temp1");
                  printf("Pinterp=%g  oi=%d Y=%.13g ok=%d  ind=%d\n",
                         Pinterp, oi,Y,ok, ind);
                  printbface(bface);
                  printbface(obface);
                  grid->time  = 42;
                  write_grid(grid);
                  errorexit("Pinterp is not finite!");
                }
              }
            }
          }
          else
            errorexit("implement more dir cases!");
        }
        else
          errorexit("implement more bface->same? cases!");
      }
      else /* bface is not associated with just one other face */
      {
          /* get values by 3d interpolation from box ob */
          double *P = obox->v[iPsi]; /* values in other box ob */
          double *Pcoeffs = obox->v[Ind("DNSdata_temp1")];
          spec_Coeffs(obox, P, Pcoeffs);

          forPointList_inbox(bface->fpts, box, pi, ind)
          {
            double X = box->v[oXi][ind];
            double Y = box->v[oYi][ind];
            double Z = box->v[oZi][ind];
            double Pinterp = spec_interpolate(obox, Pcoeffs, X,Y,Z);
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
      }
    }
    else /* there is no other box */
    {
      /* set far limit BC */
      if(bface->outerbound && setFarLimit)
        forPointList_inbox(bface->fpts, box, pi, ind)
          FPsi[ind] = Psi[ind] - PsiFarLimit;
    }
  }
}


/* new main BC routine, replaces set_BNSdata_BCs__old */
void general_DNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs,
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
        if( (box->MATTR == TOUCH) && (!FakeMatterOutside) )
        {
          int f;
          /* make a list of faces to omit */
          for(f=1; f<6; f++) push_intList(skip_f, f);
        }
      }

      /* set some BCs for each box */
      //set_interbox_and_FarLimit_BCs(box, iFPsi, iPsi, iPsix,iPsiy,iPsiz,
      //                              PsiFarLimit, 1, skip_f);
      DNS_set_interbox_and_outer_BCs(box, iFPsi, iPsi, iPsix,iPsiy,iPsiz,
                                     PsiFarLimit, 1, skip_f);
// replace above by DNS_set_interbox_and_outer_BCs

    } /* end forallboxes */
    free_intList(skip_f);

    Incr_vindDerivs:
      /* increase index for derivs */
      vindDerivs += 3;
      if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}

/* standard BCs for all fields */
/* set BCs for a varlist */
void set_DNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs,
                     int nonlin)
{
    general_DNSdata_BCs(vlFu, vlu, vluDerivs, nonlin);
}


/* impose: DNSdata_Sigma[i] = Omega*(xc1-xCM) * y
   or:     DNSdata_Sigma[i] = Omega*(xc2-xCM) * y
   on outside (lam=1) of TOUCH boxes */
void set_Sigma_Omega_r_y_BCs(tVarList *vlFu, tVarList *vlu,
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
    if(!strstr(varname, "DNSdata_Sigma")) continue;

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
