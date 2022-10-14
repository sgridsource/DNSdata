/* DNS_Interpolate_ADMvars.c */
/* Wolfgang Tichy 10/2009 */
/* interpolate onto points which we read from a file */

#include "sgrid.h"
#include "DNSdata.h"



/* position filepointer below the label str */
int position_fileptr_below_label(FILE *in, char *str)
{
  char line[1000];
  
  while(fgets(line, 999, in)!=NULL)
  {
    if(strstr(line, str)!=NULL) return 1; //break;
  }
  return EOF;
}

/* read xyz of next point */
int read_next_xyz_from_pointsfile(FILE *in, double *x, double *y, double *z)
{
  double xyz[3];
  int n;

  n=fread(xyz, sizeof(double), 3, in);
  if(n!=3) n=EOF;
  *x=xyz[0];
  *y=xyz[1];
  *z=xyz[2];
  return n;
}


/* interpolate ADM initial data onto points listed in pointsfile */
/* Note: here we check extra pars that are not introduced in sgrid_DNSdata.c :
   BNSdata_Interpolate_pointsfile = some_filename
   BNSdata_Interpolate_output = some_output_filename
   BNSdata_Interpolate_verbose = no # yes
   BNSdata_Interpolate_max_xyz_diff = 1e-8
   The last 2 are just for debugging. Their values are set e.g. by BAM.
   We do nothing if BNSdata_Interpolate_pointsfile = *** NONE ***        */
int DNS_Interpolate_ADMvars(tGrid *grid)
{
  int pr = GetvLax("BNSdata_Interpolate_verbose", "yes");
  FILE *in, *out;
  char *pointsfile;
  char *outfile;
  int nelem = 3456000; /* number of elements here 5*5*5*3*3*3*1024 */
  double *inbuf;
  double *outbuf;
  int insize, outsize, bpos, ncount, wcount;
  tVarList *vlu;
  tVarList *vlc;
  int corot1 = Getv("DNSdata_rotationstate1","corotation");
  int corot2 = Getv("DNSdata_rotationstate2","corotation");
  int b, j;

  if(GetsLax("BNSdata_Interpolate_pointsfile")==0) return 0;
  if(strcmp(Gets("BNSdata_Interpolate_pointsfile"), "*** NONE ***")==0)
    return 0;
  prdivider(0);
  printf("DNS_Interpolate_ADMvars:\n");
  prTimeIn_s("WallTime: ");

  /* allocate varlists */
  vlu = vlalloc(grid);
  vlc = vlalloc(grid);

  /* add all vars to vlu */
  vlpush(vlu, Ind("alpha"));
  vlpush(vlu, Ind("DNSdata_Bx"));
  vlpush(vlu, Ind("gxx"));
  vlpush(vlu, Ind("Kxx"));
  vlpush(vlu, Ind("DNSdata_q"));
  vlpush(vlu, Ind("DNSdata_VRx"));
  /* Those we do not need:
    vlpush(vlu, Ind("DNSdata_Sigmax"));
    vlpush(vlu, Ind("DNSdata_wBx"));
    vlpush(vlu, Ind("rho"));
    vlpush(vlu, Ind("jx"));
    vlpush(vlu, Ind("Sxx"));
  */

  /* now duplicate vlu to get vlc */  
  vlc = AddDuplicateEnable(vlu, "_c");

  /* set buffer sizes */
  /* size of outbuf depends on number of vars */
  insize  = 3 * sizeof(double);
  outsize = vlc->n * sizeof(double);
  inbuf  = calloc(nelem,  insize);
  outbuf = calloc(nelem, outsize);
  if(inbuf==NULL)  errorexit("couldn't allocate inbuf");
  if(outbuf==NULL) errorexit("couldn't allocate outbuf");

  /* write coeffs of vlu in all boxes into vlc */
  forallboxes(grid, b)
    spec_Coeffs_varlist(grid->box[b], vlu, vlc);

  /* filenames */
  pointsfile = Gets("BNSdata_Interpolate_pointsfile");
  outfile    = Gets("BNSdata_Interpolate_output");
  printf(" BNSdata_Interpolate_pointsfile = %s\n", pointsfile);
  printf(" BNSdata_Interpolate_output = %s\n", outfile);
  prTimeIn_s(" WallTime: ");

  /* open both files */
  in = fopen(pointsfile, "rb");
  if(!in) errorexits("failed opening %s", pointsfile);
  out = fopen(outfile, "wb");
  if(!out) errorexits("failed opening %s", outfile);


  /* write header info */  
  fprintf(out, "%s", "#");
  for(j=0; j<vlu->n; j++)
    fprintf(out, " %s", VarName(vlu->index[j]));
  fprintf(out, "%s", "\n");
  fprintf(out, "%s\n", "$BEGIN_data:");
  
  j=position_fileptr_below_label(in, "$BEGIN_data:");
  if(j==EOF) errorexits("could not find $BEGIN_data: in %s", pointsfile);

  do /* read buffer as long as there is data */
  {
    ncount = fread(inbuf, insize, nelem, in);
    SGRID_TOPLEVEL_Pragma(omp parallel)
    {
      tGrid *grid_p = make_empty_grid(grid->nvariables, 0);
      copy_grid(grid, grid_p, 0);
      intList *bl = alloc_intList(); /* list that contains boxes */

      SGRID_TOPLEVEL_Pragma(omp for)
      for(bpos=0; bpos<ncount; bpos++)
      {
        double x,y,z, val;
        double gxx=-1e300;
        double X,Y,Z;
        int ind, j, b;
        int star;
        int inind  =      3 * bpos;
        int outind = vlc->n * bpos;

        /* get x,y,z */
        x = inbuf[inind];
        y = inbuf[inind+1];
        z = inbuf[inind+2];
        if(pr) printf("(x,y,z)=(%g,%g,%g)\n", x,y,z);

        /* initial guess for X,Y,Z, b: */
        if(x>=0.0) star=STAR1;
        else       star=STAR2;
        clear_intList(bl);
        bladd_ifAttrib(grid, iSIDE, star, bl); 
        bladd_ifAttrib(grid, iSIDE, ZERO, bl); /* boxes with x>0 and x<0 */
        /*
        errorexit("bl is missing boxes that do not belong to just one star");
        nearest_b_XYZ_of_xyz_inboxlist(grid, bl->e,bl->n,
                                       &b, &ind, &X,&Y,&Z, x,y,z);
        if(grid->box[b]->COORD==CART)
        {
          X=x;
          Y=y;
          Z=z;
        }
        if(pr) printf("guess:  b=%d (X,Y,Z)=(%g,%g,%g)  nearest ind=%d\n", b, X,Y,Z, ind);
        */
        /* get X,Y,Z, b of x,y,z */
        b=DNSgrid_Get_BoxAndCoords_of_xyz(grid_p, &X,&Y,&Z, bl, x,y,z);
        if(pr) printf("actual: b=%d (X,Y,Z)=(%g,%g,%g)\n", b, X,Y,Z);
        if(b<0)
        {
          printf("point: (x,y,z)=(%g,%g,%g)\n", x,y,z);
          printf("error: b=%d (X,Y,Z)=(%g,%g,%g)\n", b, X,Y,Z); 
          errorexit("could not find point");
        }
        if( GetsLax("BNSdata_Interpolate_max_xyz_diff")!=0 )
        {
          double tol= Getd("BNSdata_Interpolate_max_xyz_diff");
          tBox *box = grid_p->box[b];
          double xs,ys,zs, diff;

          if(box->x_of_X[1]==NULL) { xs=X;  ys=Y;  zs=Z; }
          else
          {
            xs = box->x_of_X[1]((void *) box, -1, X,Y,Z);
            ys = box->x_of_X[2]((void *) box, -1, X,Y,Z);
            zs = box->x_of_X[3]((void *) box, -1, X,Y,Z);
          }
          diff = sqrt( (x-xs)*(x-xs) + (y-ys)*(y-ys) + (z-zs)*(z-zs) );
          if(diff>tol)
          {
            printf("error: b=%d {X,Y,Z}={%.15g,%.15g,%.15g}\n", b, X,Y,Z);
            printf("=> {xs,ys,zs}={x(X,Y,Z), y(X,Y,Z), z(X,Y,Z)} with:\n");
            printf(" {xs,ys,zs}={%.15g,%.15g,%.15g}\n", xs,ys,zs);
            printf("BUT {x,y,z}={%.15g,%.15g,%.15g}\n", x,y,z);
            printf(" diff=%.15g\n", diff);
            fflush(stdout);
            errorexit("x(X,Y,Z), y(X,Y,Z), z(X,Y,Z) are wrong!");
          }
        }

        /* interpolate vlu (using coeffs in vlc) to X,Y,Z in box b */
        for(j=0; j<vlc->n; j++)
        {
          tBox *box = grid_p->box[b];
          double *c = box->v[vlc->index[j]];

          /* HACK: don't interpolate for some vars that are zero
             or already known */
          if( strcmp(VarName(vlu->index[j]),"gxx")==0 ) 
            gxx = val = spec_interpolate(box, c, X,Y,Z);
          else if( strcmp(VarName(vlu->index[j]),"gyy")==0 ||
                   strcmp(VarName(vlu->index[j]),"gzz")==0   )
            val = gxx;
          else if( strcmp(VarName(vlu->index[j]),"gxy")==0 ||
                   strcmp(VarName(vlu->index[j]),"gxz")==0 ||
                   strcmp(VarName(vlu->index[j]),"gyz")==0   )
            val=0.0;
          else if( strcmp(VarName(vlu->index[j]),"DNSdata_q")==0 )
          {
            if(box->MATTR!=INSIDE) val=0.0;
            else             val=spec_interpolate(box, c, X,Y,Z);
          }
          else if( strcmp(VarName(vlu->index[j]),"DNSdata_VRx")==0 ||
                   strcmp(VarName(vlu->index[j]),"DNSdata_VRy")==0 ||
                   strcmp(VarName(vlu->index[j]),"DNSdata_VRz")==0   )
          {
            if(box->MATTR!=INSIDE) val=0.0;
            else if( (box->SIDE==STAR1) && corot1 ||
                     (box->SIDE==STAR2) && corot2 )  val=0.0;
            else val=spec_interpolate(box, c, X,Y,Z);
          }
          else val = spec_interpolate(box, c, X,Y,Z);
          /* if we always interpolate we need:
          val = spec_interpolate(box, c, X,Y,Z); */
          if(!finit(val))
          {
            printf("point:  (x,y,z)=(%g,%g,%g)\n", x,y,z);
            printf("NAN at: b=%d (X,Y,Z)=(%g,%g,%g)\n", b, X,Y,Z); 
            errorexit("spec_interpolate returned NAN, probably (X,Y,Z) was bad!");
          }
          if(pr) printf("%s=%g\n", VarName(vlu->index[j]), val);
          /* put val into outbuf */
          outbuf[outind + j] = val;
        }
      }
      /* free local copies */
      free_intList(bl);
      free_grid(grid_p);
    }
    wcount = fwrite(outbuf, outsize, ncount, out);
    if(wcount!=ncount)  errorexit("write count is wrong");
  } while(ncount==nelem); /* if we could fill buffer once, try again */

  /* close files */
  fclose(out);
  fclose(in);
  
  /* free var and buffers */
  vlfree(vlu);
  VLDisableFree(vlc);
  free(outbuf);
  free(inbuf);

  printf("DNS_Interpolate_ADMvars: finished interpolations "
         "and wrote all output.\n");
  prdivider(0);
  prTimeIn_s("WallTime: ");
     
  return 0;
}
