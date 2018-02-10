/* set BCs */
/* Wolfgang Tichy 2010 */


#include "sgrid.h"
#include "DNSdata.h"

#define Power pow


/* functions in this file */
/* ... */



/* old intricate version of BCs */
void set_DNSdata_BCs__old(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin)
{
  tGrid *grid = vlu->grid;
  int b;
  int vind;
  int vindDerivs=0;
  tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vlu->vlPars;

  if( grid->box[0]->n1 != grid->box[1]->n1 ||
      grid->box[3]->n1 != grid->box[2]->n1 ||
      grid->box[1]->n1 != grid->box[2]->n1 ) 
    errorexit("all n1 in boxes0-3 must be the same because we currently use "
              "lines like:\n"
              "FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];\n"
              "where Psi = grid->box[0]->v[vi], P = grid->box[1]->v[vi]");
                      
  for(vind=0; vind<vlu->n; vind++)
  {
    int ncomp = VarNComponents(vlu->index[vind]);
    double PsiFarLimit = VarFarLimit(vlu->index[vind])*nonlin;
    int BCs_atInf = 1;
    double FarLimitBC_R = Getd("DNSdata_FarLimitBC_Radius");
    int BCs_box1_2 = 1;
    int BCs_AxisAtOuterInterfaces = 1;
    int BCs_box0_A0 = 1; /* set BCs at A=0 in box0 */
    int BCs_box3_A0 = 1; /* set BCs at A=0 in box3 */
    //int BCs_box0_1 = 1;
    //int BCs_box3_2 = 1;
    char *varname = VarName(vlu->index[vind]);

    /* do nothing and goto end of loop if var with vind is not the one 
       of the current block */ 
    if(blkinfo!=NULL) if(vlu->index[vind] != blkinfo->vari) 
                        goto Inc_vindDerivs;

    if(strstr(varname, "DNSdata_Sigma"))
    { 
      BCs_atInf = 0;
      FarLimitBC_R = -1.0;
      BCs_box1_2 = 0;
      BCs_AxisAtOuterInterfaces = 0;
      BCs_box0_A0 = 0;
      BCs_box3_A0 = 0;
      /* printf("varname=%s\n", varname); */
    }

    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[vlFu->index[vind]];
      double *Psi  = box->v[vlu->index[vind]];
      double *Psix = box->v[vluDerivs->index[vindDerivs]];
      double *Psiy = box->v[vluDerivs->index[vindDerivs+1]];
      double *Psiz = box->v[vluDerivs->index[vindDerivs+2]];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;
      int i,j,k;

      /* do nothing and continue if current block is not in box b */
      if(blkinfo!=NULL) if(b!=blkinfo->bi) continue;

      /* BCs */
      if(Getv("DNSdata_grid", "SphericalDF"))
      {
        forplane1(i,j,k, n1,n2,n3, 0)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 1.0*(vind+1)*nonlin;

        forplane1(i,j,k, n1,n2,n3, n1-1)
          FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - 0.5*(vind+1)*nonlin;
      }
      else if (Getv("DNSdata_grid", "AnsorgNS"))
      {
        double *P;
        double *dP[4];
        double *BM;

        /* special rho=0 case??? */
        if(b==0 || b==1 || b==2 || b==3)
        {
          int pl;
          char str[1000];
          snprintf(str, 999, "box%d_basis2", b);
          if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
          {
            double *Psi_phi_phi = box->v[Ind("DNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("DNSdata_temp2")];
            double *temp3 = box->v[Ind("DNSdata_temp3")];
            double *temp4 = box->v[Ind("DNSdata_temp4")];

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* loop over rho=0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B=0 and B=1 */
              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                  FPsi[Index(i,j,k)] = Psi_phi_phi[Index(i,j,k)];
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                {
                  double Psi_rho = Psiy[Index(i,j,k)];
                  double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
                  FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
                }
              }
          }
          /* same as before, but also interpolate to rho=0 */
          else if(Getv("DNSdata_regularization", "regularity_on_axis"))
          {
            double *Psi_phi_phi = box->v[Ind("DNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("DNSdata_temp2")];
            double *temp3 = box->v[Ind("DNSdata_temp3")];
            double *temp4 = box->v[Ind("DNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              forplane2(i,j,k, n1,n2,n3, pl)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
          /* same as before, but do it only in box0/3 at A,B=1,0 and A,B=1,1 */
          else if(Getv("DNSdata_regularization", 
                       "regularity_on_axis_at_center") && (b==0 || b==3) )
          {
            double *Psi_phi_phi = box->v[Ind("DNSdata_temp1")];
            double *Psi_y_phi_phi = box->v[Ind("DNSdata_temp2")];
            double *temp3 = box->v[Ind("DNSdata_temp3")];
            double *temp4 = box->v[Ind("DNSdata_temp4")];
            double *line = (double *) calloc(n2, sizeof(double));
            double *BM[2];
            BM[0] = (double *) calloc(n2, sizeof(double));
            BM[1] = (double *) calloc(n2, sizeof(double));

            /* get u_phi_phi */
            spec_Deriv2(box, 3, Psi, Psi_phi_phi);
            
            /* get u_rho_phi_phi at phi=0 */
            /* d/drho = dx^i/drho d/dx^i, 
               dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
               ==> d/drho u = d/dy u  at phi=0 */           
            /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
            cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

            /* obtain BM vectors for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

            /* loop over rho~0 boundary */
            for(pl=0; pl<n2; pl=pl+n2-1)  /* <-- B~0 and B~1 */
            {
              int l;
              double U0, V0;

              i=n1-1;   /* do it only at A=1 */
              j=pl;
              for(k=0; k<n3; k++)
              {
                if(k>0) /* phi>0: impose u_phi_phi=0 */
                {
                  /* find value Psi_phi_phi at B=0 or 1 */
                  get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
                  FPsi[Index(i,j,k)] = U0;
                }
                else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
                { /* Psi_rho = Psiy  
                     Psi_rho_phi_phi = Psi_y_phi_phi */
                  /* find value Psi_rho at B=0 or 1 */
                  get_memline(Psiy, line, 2, i,k, n1,n2,n3);
                  for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

                  /* find value Psi_rho_phi_phi at B=0 or 1 */
                  get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
                  for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

                  FPsi[Index(i,j,k)] = U0 + V0;
                }
              }
            }
            free(BM[0]);
            free(BM[1]);
            free(line);
          }
        } /* end: special rho=0 case??? */

        if(b==0)  /* in box0 */
        {
          BM = (double *) calloc(n1, sizeof(double));

          /* values at A=0 are equal in box0 and box1 */
          P = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          spec_Basis_times_CoeffMatrix_direc(box, 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            int l;
            double *line = (double *) calloc(n1, sizeof(double));
            double U0;

            /* find values U0 in box0 at A=0*/
            get_memline(Psi, line, 1, j,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
            free(line);
          }
          free(BM);
        }
        else if(b==3)  /* in box3 */
        {
          BM = (double *) calloc(n1, sizeof(double));

          /* values at A=0 are equal in box3 and box2 */
          P = grid->box[2]->v[vlu->index[vind]]; /* values in box2 */
          spec_Basis_times_CoeffMatrix_direc(box, 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            int l;
            double *line = (double *) calloc(n1, sizeof(double));
            double U0;

            /* find values U0 in box3 at A=0*/
            get_memline(Psi, line, 1, j,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n1; l++)  U0 += BM[l]*line[l];
            FPsi[Index(i,j,k)] = U0 - P[Index(i,j,k)];
            free(line);
          }
          free(BM);
        }
        else if(b==1)
        {
          BM = (double *) calloc(max3(grid->box[0]->n1, box->n1, box->n2), 
                                 sizeof(double));

          /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
          dP[1] = grid->box[2]->v[vluDerivs->index[vindDerivs]];
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];

          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vindDerivs+2]];
          spec_Basis_times_CoeffMatrix_direc(grid->box[0], 1, BM, 0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];
            double DP[4];
            int m,l;
            /* find derivs of Psi at A=0 in box0 and store them in DP[m] */
            for(m=1; m<=3; m++)
            {
              int n1 = grid->box[0]->n1;
              double *line = (double *) calloc(n1, sizeof(double));

              DP[m] = 0.0;
              get_memline(dP[m], line, 1, j,k, n1,n2,n3);
              for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
              free(line);
            }
            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          }

          /* Psi=0 at infinity */
          if(Getv("box1_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
          }

          free(BM);
        }
        else if(b==2)
        {
          BM = (double *) calloc(max3(grid->box[3]->n1, box->n1, box->n2), 
                                 sizeof(double));

          /* values at A=1 are equal in box1 and box2 */
          P  = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
          /* Below we use the approximate normal vec 
             ( -cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[3]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vindDerivs+2]];
          spec_Basis_times_CoeffMatrix_direc(grid->box[3], 1, BM, 0.0);
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];
            double DP[4];
            int m,l;
            /* find derivs of Psi at A=0 in box3 and store them in DP[m] */
            for(m=1; m<=3; m++)
            {
              int n1 = grid->box[3]->n1;
              double *line = (double *) calloc(n1, sizeof(double));

              DP[m] = 0.0;
              get_memline(dP[m], line, 1, j,k, n1,n2,n3);
              for(l=0; l<n1; l++)  DP[m] += BM[l]*line[l];
              free(line);
            }
            FPsi[Index(i,j,k)] = 
              -cos(PI*B)         * (Psix[Index(i,j,k)] - DP[1]) +
               sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - DP[2]) +
               sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - DP[3]);
          }

          /* Psi=0 at infinity */
          if(Getv("box2_basis2", "ChebExtrema"))
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
          else // B=0 is not on grid for ChebZeros!!!
          {
            int l;
            double U0;
            double *line = (double *) calloc(n2, sizeof(double));

            /* obtain BM vector for interpolation along B */
            spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
            for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
            {
              /* find value of Psi at A=1, B=0 */
              get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
              for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
              FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
            }
            free(line);
          }

          free(BM);
        }
        else errorexiti("b=%d should be impossible!", b);
      } /* end: else if (Getv("DNSdata_grid", "AnsorgNS")) */
      else if (Getv("DNSdata_grid", "4ABphi_2xyz"))
      {
        double *P;
        double *dP[4];
        double *X, *Y, *Z,  *xp, *yp, *zp;
        double *Pcoeffs;
        double Pinterp;
        double x,y,z;

        /* special rho=0 case??? */
        DNSdata_RegularityConditions_for_Var_at_rho_eq_0(box, FPsi,
                                                         Psi, Psix,Psiy,Psiz);
        /* cases for each box */
        if(b==0)  /* in box0 */
        {
          /* values at A=0 are equal in box0 and box1 */
          if(BCs_box0_A0)
          {
            P = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
            forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
              FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          }
          /* values at A=Amax are interpolated from box5 */
          xp = box->v[Ind("x")];
          yp = box->v[Ind("y")];
          zp = box->v[Ind("z")];
          P = grid->box[5]->v[vlu->index[vind]]; /* values in box5 */
          Pcoeffs = grid->box[5]->v[Ind("DNSdata_temp1")];
          spec_Coeffs(grid->box[5], P, Pcoeffs);
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amax */
          {
            int ind=Index(i,j,k);
            x = xp[ind]; 
            y = yp[ind]; 
            z = zp[ind];
            Pinterp = spec_interpolate(grid->box[5], Pcoeffs, x,y,z);
            if(!finite(Pinterp))
            {
              printf("Pinterp=%g  x=%.13g y=%.13g z=%.13g  ind=%d\n", 
                     Pinterp, x,y,z, ind);
              NumberChecker_CheckIfFinite(grid, "DNSdata_temp1");
              printbox(grid->box[5]);
              grid->time  = 42;
              write_grid(grid);
              errorexit("Pinterp is not finite!");
            }
            FPsi[ind] = Psi[ind] - Pinterp;
          }
        }
        else if(b==3)  /* in box3 */
        {
          /* values at A=0 are equal in box3 and box2 */
          if(BCs_box3_A0)
          {
            P = grid->box[2]->v[vlu->index[vind]]; /* values in box2 */
            forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
              FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          }
          /* values at A=Amax are interpolated from box4 */
          xp = box->v[Ind("x")];
          yp = box->v[Ind("y")];
          zp = box->v[Ind("z")];
          P = grid->box[4]->v[vlu->index[vind]]; /* values in box4 */
          Pcoeffs = grid->box[4]->v[Ind("DNSdata_temp1")];
          spec_Coeffs(grid->box[4], P, Pcoeffs);
          forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=Amax */
          {
            int ind=Index(i,j,k);
            x = xp[ind]; 
            y = yp[ind]; 
            z = zp[ind]; 
            Pinterp = spec_interpolate(grid->box[4], Pcoeffs, x,y,z);
            FPsi[ind] = Psi[ind] - Pinterp;
          }
        }
        else if(b==1)  /* in box1 */
        {
          /* normal derivs (d/dx) at A=1 are equal in box1 and box2 */
          if(BCs_box1_2)
          {
            dP[1] = grid->box[2]->v[vluDerivs->index[vindDerivs]];
            forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
              FPsi[Index(i,j,k)] = Psix[Index(i,j,k)] - dP[1][Index(i,j,k)];
          }
          /* normal derivs (~d/dA) at A=0 are equal in box1 and box0 */
          /* Below we use the approximate normal vec 
             ( cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[0]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[0]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[0]->v[vluDerivs->index[vindDerivs+2]];
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];

            FPsi[Index(i,j,k)] = 
             cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
             sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
             sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          }
          
          /* Psi=0 at infinity */
          if(BCs_atInf)
          {
            if(Getv("box1_basis2", "ChebExtrema"))
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
                FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
            else // B=0 is not on grid for ChebZeros!!!
            {
              int l;
              double U0;
              double *BM = (double *) calloc(n2, sizeof(double));
              double *line = (double *) calloc(n2, sizeof(double));

              /* obtain BM vector for interpolation along B */
              spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
              {
                /* find value of Psi at A=1, B=0 */
                get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
                FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
              }
              free(line);
              free(BM);
            }
          }

          /* Psi=FarLimit if r>FarLimitBC_R */
          if(FarLimitBC_R>0.0)
            DNSdata_FarLimit_Beyond_R(box, FPsi, Psi,Psix,Psiy,Psiz,
                                      PsiFarLimit, FarLimitBC_R);
        }
        else if(b==2)  /* in box2 */
        {
          /* values at A=1 are equal in box1 and box2 */
          if(BCs_box1_2)
          {
            P  = grid->box[1]->v[vlu->index[vind]]; /* values in box1 */
            forplane1(i,j,k, n1,n2,n3, n1-1) /* <-- A=1 */
              FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];
          }
          /* normal derivs (d/d?) at A=0 are equal in box2 and box3 */
          /* Below we use the approximate normal vec 
             ( -cos(PI*B), sin(PI*B)*cos(phi), sin(PI*B)*sin(phi) ) */
          dP[1] = grid->box[3]->v[vluDerivs->index[vindDerivs]];
          dP[2] = grid->box[3]->v[vluDerivs->index[vindDerivs+1]];
          dP[3] = grid->box[3]->v[vluDerivs->index[vindDerivs+2]];
          forplane1(i,j,k, n1,n2,n3, 0) /* <-- A=0 */
          {
            double B   = box->v[Ind("Y")][Index(i,j,k)];
            double phi = box->v[Ind("Z")][Index(i,j,k)];

            /* NOTE: The minus sign in front of cos(PI*B) below is correct.
                     It used to be +cos(PI*B) !!! But this leads to a vector
                     that can be tangent at some B in (0,0.5), which
                     would be bad for a normal vector! */ 
            FPsi[Index(i,j,k)] = 
              -cos(PI*B)         * (Psix[Index(i,j,k)] - dP[1][Index(i,j,k)]) +
               sin(PI*B)*cos(phi)* (Psiy[Index(i,j,k)] - dP[2][Index(i,j,k)]) +
               sin(PI*B)*sin(phi)* (Psiz[Index(i,j,k)] - dP[3][Index(i,j,k)]);
          }

          /* Psi=0 at infinity */
          if(BCs_atInf)
          {
            /* Psi=0 at infinity */
            if(Getv("box2_basis2", "ChebExtrema"))
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)=(1,0) */
                FPsi[Index(n1-1,0,k)] = Psi[Index(n1-1,0,k)]-PsiFarLimit;
            else // B=0 is not on grid for ChebZeros!!!
            {
              int l;
              double U0;
              double *BM = (double *) calloc(n2, sizeof(double));
              double *line = (double *) calloc(n2, sizeof(double));

              /* obtain BM vector for interpolation along B */
              spec_Basis_times_CoeffMatrix_direc(box, 2, BM, 0.0);
              for(k=0;k<n3;k++)  /* <--loop over all phi with (A,B)~(1,0) */
              {
                /* find value of Psi at A=1, B=0 */
                get_memline(Psi, line, 2, n1-1,k, n1,n2,n3);
                for(U0=0.0, l=0; l<n2; l++)  U0 += BM[l]*line[l];
                FPsi[Index(n1-1,0,k)] = U0-PsiFarLimit;
              }
              free(line);
              free(BM);
            }
          }

          /* Psi=FarLimit if r>FarLimitBC_R */
          if(FarLimitBC_R>0.0)
            DNSdata_FarLimit_Beyond_R(box, FPsi, Psi,Psix,Psiy,Psiz,
                                      PsiFarLimit, FarLimitBC_R);
        }
        else if(b==5)  /* in box5 */
        {
          /* values at border are interpolated from box0 */
          double A,B,phi;
          int pl; //, k_phi;
          double *pA = box->v[Ind("DNSdata_A")];
          double *pB = box->v[Ind("DNSdata_B")];
          double *pphi = box->v[Ind("DNSdata_phi")];
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[0]->v[vlu->index[vind]]; /* values in box0 */
          Pcoeffs = grid->box[0]->v[Ind("DNSdata_temp1")];
          spec_Coeffs(grid->box[0], P, Pcoeffs);
          for(pl=0; pl<n1; pl=pl+n1-1)
          {
            int ind=Index(pl,0,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n2; pl=pl+n2-1)
          {
            int ind=Index(0,pl,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];
                             
              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n3; pl=pl+n3-1)
          {
            int ind=Index(0,0,pl);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[0]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[0], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[0], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[0], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else if(b==4)  /* in box4 */
        {
          /* values at border are interpolated from box3 */
          double A,B,phi;
          int pl; //, k_phi;
          double *pA = box->v[Ind("DNSdata_A")];
          double *pB = box->v[Ind("DNSdata_B")];
          double *pphi = box->v[Ind("DNSdata_phi")];
          X = box->v[Ind("X")];
          Y = box->v[Ind("Y")];
          Z = box->v[Ind("Z")];
          P = grid->box[3]->v[vlu->index[vind]]; /* values in box3 */
          Pcoeffs = grid->box[3]->v[Ind("DNSdata_temp1")];
          spec_Coeffs(grid->box[3], P, Pcoeffs);
          for(pl=0; pl<n1; pl=pl+n1-1)
          {
            int ind=Index(pl,0,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane1_nojump(i,j,k, n1,n2,n3, pl) /* <-- x=xmin and xmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n2; pl=pl+n2-1)
          {
            int ind=Index(0,pl,0);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane2_nojump(i,j,k, n1,n2,n3, pl) /* <-- y=ymin and ymax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
          for(pl=0; pl<n3; pl=pl+n3-1)
          {
            int ind=Index(0,0,pl);
            //int i0;
            //phi   = Arg(Y[ind],Z[ind]);   if(phi<0) phi = 2.0*PI+phi;
            //k_phi = grid->box[3]->n3 * phi/(2.0*PI);
            //nearestXYZ_of_xyz_inplane(grid->box[3], &i0, &A,&B,&phi,
            //                          X[ind],Y[ind],Z[ind], 3, k_phi);
            forplane3_nojump(i,j,k, n1,n2,n3, pl) /* <-- z=zmin and zmax */
            {
              ind=Index(i,j,k);
              //compute_ABphi_from_xyz(grid->box[3], &A,&B,&phi, X[ind],Y[ind],Z[ind]);
              A=pA[ind];  B=pB[ind];  phi=pphi[ind];

              Pinterp = spec_interpolate(grid->box[3], Pcoeffs, A,B,phi);
              FPsi[ind] = Psi[ind] - Pinterp;
            }
          }
        }
        else errorexiti("b=%d should be impossible!", b);

        /* special rho=0 case again at A=0 and A=1 ???  */
        if( ( (b==0 || b==3) || 
              (b==1 || b==2) && BCs_AxisAtOuterInterfaces ) && 
            Getv("DNSdata_regularization",
                 "regularity_on_axis_at_interfaces") )
        {
          DNSdata_RegularityConditions_for_Var_at_rho_eq_0(box, FPsi,
                                                           Psi, Psix,Psiy,Psiz);
        } /* end: special rho=0 case again ??? */

      } /* end: else if (Getv("DNSdata_grid", "4ABphi_2xyz")) */

    } /* end forallboxes */

    Inc_vindDerivs:
    /* increase index for derivs */
    vindDerivs += 3;
    if(VarComponent(vlu->index[vind])==ncomp-1) vindDerivs += 6*ncomp;
  } /* end loop over vars */
}


/* treat rho=0 case: */
void DNSdata_RegularityConditions_for_Var_at_rho_eq_0(tBox *box, double *FPsi,
                        double *Psi, double *Psix, double *Psiy, double *Psiz)
{
  int b = box->b;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int i,j,k;

  /* special rho=0 case??? */
  if(box->COORD != CART)
  {
    int pl, pl0, pllim;
    char str[1000];

    /* if B=0 or B=1 are not in box adjust plane loop for rho=0 boundary */
    if( dequal(box->bbox[2], 0.0) )  pl0 = 0;
    else                             pl0 = n2-1;
    if( dequal(box->bbox[3], 1.0) )  pllim = n2;
    else                             pllim = 1;

    snprintf(str, 999, "box%d_basis2", b);
    if(Getv(str, "ChebExtrema"))  /* treat rho=0 case */
    {
      double *Psi_phi_phi = box->v[Ind("DNSdata_temp1")];
      double *Psi_y_phi_phi = box->v[Ind("DNSdata_temp2")];
      double *temp3 = box->v[Ind("DNSdata_temp3")];
      double *temp4 = box->v[Ind("DNSdata_temp4")];

      /* get u_phi_phi */
      spec_Deriv2(box, 3, Psi, Psi_phi_phi);
      
      /* get u_rho_phi_phi at phi=0 */
      /* d/drho = dx^i/drho d/dx^i, 
         dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
         ==> d/drho u = d/dy u  at phi=0 */           
      /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
      cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

      /* loop over rho=0 boundary */
      for(pl=pl0; pl<pllim; pl=pl+n2-1)  /* <-- B=0 and B=1 */
        forplane2(i,j,k, n1,n2,n3, pl)
        {
          if(k>0) /* phi>0: impose u_ijk = u_ij0 (not u_phi_phi=0) */
            FPsi[Index(i,j,k)] = Psi[Index(i,j,k)]-Psi[Index(i,j,0)];
          else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
          {
            double Psi_rho = Psiy[Index(i,j,k)];
            double Psi_rho_phi_phi = Psi_y_phi_phi[Index(i,j,k)];
            FPsi[Index(i,j,k)] = Psi_rho + Psi_rho_phi_phi;
          }
        }
    }
    /* same as before, but also interpolate to rho=0 */
    else if(Getv("DNSdata_regularization", "regularity_on_axis"))
    {
      double *Psi_phi_phi = box->v[Ind("DNSdata_temp1")];
      double *Psi_y_phi_phi = box->v[Ind("DNSdata_temp2")];
      double *temp3 = box->v[Ind("DNSdata_temp3")];
      double *temp4 = box->v[Ind("DNSdata_temp4")];
      double *line = (double *) calloc(n2, sizeof(double));
      double *BM[2];
      BM[0] = (double *) calloc(n2, sizeof(double));
      BM[1] = (double *) calloc(n2, sizeof(double));

      /* get u_phi_phi */
      spec_Deriv2(box, 3, Psi, Psi_phi_phi);
      
      /* get u_rho_phi_phi at phi=0 */
      /* d/drho = dx^i/drho d/dx^i, 
         dx/drho=0, dy/drho=cos(phi), dz/drho=sin(phi)
         ==> d/drho u = d/dy u  at phi=0 */           
      /* get u_rho_phi_phi at phi=0: u_rho_phi_phi = d/dy u_phi_phi */
      cart_partials(box, Psi_phi_phi, temp3, Psi_y_phi_phi, temp4);

      /* obtain BM vectors for interpolation along B */
      spec_Basis_times_CoeffMatrix_direc(box, 2, BM[0], 0);
      spec_Basis_times_CoeffMatrix_direc(box, 2, BM[1], 1);

      /* loop over rho~0 boundary */
      for(pl=pl0; pl<pllim; pl=pl+n2-1)  /* <-- B~0 and B~1 */
      {
        int l;
        double U0, V0;

        forplane2(i,j,k, n1,n2,n3, pl)
        {
          if(k>0) /* phi>0: impose u_phi_phi=0 */
          {
            /* find value Psi_phi_phi at B=0 or 1 */
            get_memline(Psi_phi_phi, line, 2, i,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];
            FPsi[Index(i,j,k)] = U0;
          }
          else /* phi=0: impose u_rho + u_rho_phi_phi=0 */
          { /* Psi_rho = Psiy  
               Psi_rho_phi_phi = Psi_y_phi_phi */
            /* find value Psi_rho at B=0 or 1 */
            get_memline(Psiy, line, 2, i,k, n1,n2,n3);
            for(U0=0.0, l=0; l<n2; l++)  U0 += BM[j>0][l]*line[l];

            /* find value Psi_rho_phi_phi at B=0 or 1 */
            get_memline(Psi_y_phi_phi, line, 2, i,k, n1,n2,n3);
            for(V0=0.0, l=0; l<n2; l++)  V0 += BM[j>0][l]*line[l];

            FPsi[Index(i,j,k)] = U0 + V0;
          }
        }
      }
      free(BM[0]);
      free(BM[1]);
      free(line);
    }
  } /* end: special rho=0 case */
}

/* impose:  Psi = PsiFarLimit  for  r >= R */
void DNSdata_FarLimit_Beyond_R(tBox *box,
              double *FPsi,
              double *Psi, double *Psix, double *Psiy, double *Psiz,
              double PsiFarLimit, double R)
{
  double *xp = box->v[Ind("x")];
  double *yp = box->v[Ind("y")];
  double *zp = box->v[Ind("z")];
  double R2 = R*R;
  int i;

  forallpoints(box, i)
  {
    double r2 = xp[i]*xp[i] + yp[i]*yp[i] + zp[i]*zp[i];
    /* do nothing if r<R */
    if(r2<R2) continue;
    /* otherwise set: */ 
    FPsi[i] = Psi[i]-PsiFarLimit;
  }
}


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


/* new main BC routine, replaces set_DNSdata_BCs__old */
void general_DNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs,
                         int nonlin)
{
  tGrid *grid = vlu->grid;
  int vind;
  int vindDerivs=0;
  tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vlu->vlPars;

  if( grid->box[0]->n1 != grid->box[1]->n1 ||
      grid->box[3]->n1 != grid->box[2]->n1 ||
      grid->box[1]->n1 != grid->box[2]->n1 )
    errorexit("all n1 in boxes0-3 must be the same because we currently use "
              "lines like:\n"
              "FPsi[Index(i,j,k)] = Psi[Index(i,j,k)] - P[Index(i,j,k)];\n"
              "where Psi = grid->box[0]->v[vi], P = grid->box[1]->v[vi]");
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
    double FarLimitBC_R = Getd("DNSdata_FarLimitBC_Radius");
    char *varname = VarName(vlu->index[vind]);
    int is_Sigma = 0;

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
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *FPsi = box->v[iFPsi];
      double *Psi  = box->v[iPsi];
      double *Psix = box->v[iPsix];
      double *Psiy = box->v[iPsiy];
      double *Psiz = box->v[iPsiz];
      intList *skip_f = alloc_intList(); /* list of face we want to skip */

      /* do nothing and continue if current block is not in box b */
      if(blkinfo!=NULL) if(b!=blkinfo->bi) continue;

      /* special rho=0 case??? */
      DNSdata_RegularityConditions_for_Var_at_rho_eq_0(box, FPsi,
                                                       Psi, Psix,Psiy,Psiz);
      /* DNSdata_Sigma? */
      if(is_Sigma)
      {
        /* do nothing else for DNSdata_Sigma away from stars */
        if(box->MATTR == AWAY) continue;

        /* we need BCs for Sigma in box1/2 at A=0, but not at A=1 !!! */
        /* impose BC only at A=1 for boxes (1 and 2) that touch star */
        if(box->MATTR == TOUCH)
        {
          int f;
          /* make a list of faces to omit */
          for(f=1; f<6; f++) push_intList(skip_f, f);
        }
      }

      /* set some BCs for each box */
      set_interbox_and_FarLimit_BCs(box, iFPsi, iPsi, iPsix,iPsiy,iPsiz,
                                    PsiFarLimit, 1, skip_f);
      /* additional special BCs */
      if (Getv("DNSdata_grid", "4ABphi_2xyz"))
      {
        if(b==1 || b==2)  /* in box1/2 */
        {
          /* Psi=FarLimit if r>FarLimitBC_R */
          if(FarLimitBC_R>0.0)
            DNSdata_FarLimit_Beyond_R(box, FPsi, Psi,Psix,Psiy,Psiz,
                                      PsiFarLimit, FarLimitBC_R);
        }
      } /* end: if (Getv("DNSdata_grid", "4ABphi_2xyz")) */

      /* special rho=0 case again at A=0 and A=1 ???  */
      if( box->COORD != CART && Getv("DNSdata_regularization",
                                     "regularity_on_axis_at_interfaces") )
        DNSdata_RegularityConditions_for_Var_at_rho_eq_0(box, FPsi,
                                                         Psi, Psix,Psiy,Psiz);
      free_intList(skip_f);
    } /* end forallboxes */

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
  /* Note: general_DNSdata_BCs does not have all the options of
     set_DNSdata_BCs__old when it comes to ChebZeros and such. But for
     for the standard 4ABphi_2xyz, they both give the same results up
     to 14 digits. */
  if(Getv("DNSdata_grid", "6ABphi_2xyz") ||
     Getv("DNSdata_grid", "general_DNSdata_BCs"))
    general_DNSdata_BCs(vlFu, vlu, vluDerivs, nonlin);
  else
    set_DNSdata_BCs__old(vlFu, vlu, vluDerivs, nonlin);
}
