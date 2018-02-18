/* compute_chi.c */
/* dtim 2/2014 */

/* WT: as far as I can tell, Tim computes this for each star:
   h = q+1
   chieq   = (dh/dx)/h evaluated at (A,B,phi)=(0,1,0)
   chipole = (dh/dz)/h evaluated at (A,B,phi)=(0,Bm,pi/2)
              where Bm is chosen such that z is max
   ==> chi = |chieq/chipole|

   the max in z is called "zmin" and is found only approximately
   by looking at m=1000 points                                    */

#include "sgrid.h"
#include "DNSdata.h"
#define Log(x)     log((double) (x)))


int DNS_compute_chi(tGrid *grid)
{
	 
     int bi,i;     
     int m = 1000;
     double chieqx, chipolez, chi, Bstep, z, zmin, Bstepmin;
     double *chieq, *chipole, *temp4;
     int index_DNSdata_q       = Ind("DNSdata_q");
     int index_DNSdata_qx      = Ind("DNSdata_qx");
     int index_DNSdata_temp1   = Ind("DNSdata_temp1");
     int index_DNSdata_temp2   = Ind("DNSdata_temp2");
     int index_DNSdata_temp3   = Ind("DNSdata_temp3");
     int index_DNSdata_temp4   = Ind("DNSdata_temp4");

     /* do nothing if DNSdata_Interpolate_pointsfile exists */
     if(GetsLax("DNSdata_Interpolate_pointsfile")!=0) return 0;

printf("DNS_compute_chi needs to be updated to not use box0/3 and AnsorgNS!!!\n");
printf("*** Skipping DNS_compute_chi ***\n");
return 77;
 
     forallboxes(grid,bi) {
         
         tBox *box = grid->box[bi];
         int ijk;
         double *q   = box->v[index_DNSdata_q + 0];
         double *dq1 = box->v[index_DNSdata_qx + 0];
         double *dq2 = box->v[index_DNSdata_qx + 1];
         double *dq3 = box->v[index_DNSdata_qx + 2];
         double *h = box->v[index_DNSdata_temp1 + 0];
         double *chieq = box->v[index_DNSdata_temp2 + 0];
         double *chipole = box->v[index_DNSdata_temp3 + 0];

         FirstDerivsOf_S(box, Ind("DNSdata_q"), Ind("DNSdata_qx")); 

         forallpoints(box, ijk) { 

            h[ijk]       = q[ijk] + 1.;
            chieq[ijk]   = dq1[ijk]/h[ijk];
            chipole[ijk] = dq3[ijk]/h[ijk];

          } /* end of points loop */ 

      } /* end of boxes */


    /* Computation of chi*/
    /* for box 3*/
errorexit("DNS_compute_chi should not be done in box3");
     chieq   = grid-> box[3]->v[index_DNSdata_temp2 + 0];
     chipole = grid-> box[3]->v[index_DNSdata_temp3 + 0];
     temp4   = grid-> box[3]->v[index_DNSdata_temp4 + 0];
  
     zmin = 0;
     Bstepmin =0; 

     for (i = 0; i<m; i++ ){
            Bstep = (double) (i) / (double) (m);
            z     = z_of_NAnsorgNS3(grid-> box[3], -1, 0, Bstep, 0.5*PI);
            if(z > zmin) {
    	          zmin     = z;
    	          Bstepmin = Bstep;}
      }
             
     spec_Coeffs(grid-> box[3], chieq, temp4); 
     chieqx   = spec_interpolate(grid-> box[3], temp4, 0, 1 , 0); 
     spec_Coeffs(grid-> box[3], chipole, temp4); 
     chipolez = spec_interpolate(grid-> box[3], temp4, 0, Bstepmin , 0.5*PI); 
     chi      = chieqx/chipolez;
    	  
     printf(" Mass shedding: chi1=%e \n", chi);
     Setd("DNSdata_mass_shedding1",chi);
  
     /*for box 0*/
errorexit("DNS_compute_chi should not be done in box0");
printf("  WARNING: DNS_compute_chi results are WRONG!!!");
      chieq   = grid-> box[0]->v[index_DNSdata_temp2 + 0];
      chipole = grid-> box[0]->v[index_DNSdata_temp3 + 0];
      temp4   = grid-> box[0]->v[index_DNSdata_temp4 + 0];
 
      zmin     = 0;
      Bstepmin = 0; 
 
      for (i = 0; i<m; i++ ){
             Bstep = (double) (i) / (double) (m);
             z     = z_of_NAnsorgNS3(grid-> box[0], -1, 0, Bstep, 0.5*PI);
             if(z > zmin) {
    	            zmin     = z;
    	            Bstepmin = Bstep;}
        }
 
       spec_Coeffs(grid-> box[0], chieq, temp4); 
       chieqx   = spec_interpolate(grid-> box[0], temp4, 0, 1 , 0); 
       spec_Coeffs(grid-> box[0], chipole, temp4); 
       chipolez = spec_interpolate(grid-> box[0], temp4, 0, Bstepmin , 0.5*PI); 
       chi      = - chieqx/chipolez;  	  
    	  
       printf(" Mass shedding: chi2=%e \n", chi);
       Setd("DNSdata_mass_shedding2",chi);

    return 0;  
}  /* end of function */

