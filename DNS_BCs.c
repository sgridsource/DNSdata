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



/* new main BC routine, replaces set_interbox_and_FarLimit_BCs */
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
