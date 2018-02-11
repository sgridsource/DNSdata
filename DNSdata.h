/* DNSdata.h */
/* Wolfgang Tichy, 2007 */

/* DNSdata box attribute indices */
enum
{
  iZERO,
  iSIDE,    /* attrib. val. can be STAR1, STAR2 */
  iMATTR,   /* attrib. val. can be INSIDE, TOUCH, AWAY */
  iBOUND,   /* attrib. val. can be SSURF, ZERO */
  iCOORD    /* attrib. val. can be CART, CUBSPH, SCUBSH */
};

/* for convenience: DNSdata box attribute names */
#define SIDE   Attrib[iSIDE]  /* can be STAR1, STAR2 */
#define MATTR  Attrib[iMATTR] /* can be INSIDE, TOUCH, AWAY */
#define BOUND  Attrib[iBOUND] /* can be SSURF, ZERO */
#define COORD  Attrib[iCOORD] /* can be CART, CUBSPH, SCUBSH */

/* possible DNSdata box attribute values */
enum
{
  ZERO,     /* attrib value is zero */
  STAR1,    /* side of star1 (needs to be here in pos2) */
  STAR2,    /* side of star2 (needs to be here in pos3) */
  INSIDE,   /* star matter is inside this box */
  TOUCH,    /* star matter touches box */
  AWAY,     /* star matter is away from box */
  SSURF,    /* if BOUND==SSURF box boundary is star surface */
  CART,     /* Cartesian coords */
  CUBSPH,   /* one of the cubed spheres */
  SCUBSH,   /* stretched cubed shell that can go far out */
  NDNSATTRIBUTES
};


/* main functions */
int DNSdata_setup_boxsizes(tGrid *grid);
int DNSdata_startup(tGrid *grid);
void   DNS_compute_new_centered_q(tGrid *grid);
double DNS_compute_new_centered_q_atXYZ(tGrid *grid, int bi,
                                        double X, double Y, double Z);
int DNSdata_verify_solution(tGrid *grid);
int DNSdata_analyze(tGrid *grid);
int DNSdata_solve(tGrid *grid);
int setDNSdata(tGrid *grid);
void setADMvars(tGrid *grid);


/* more from  DNSdata.c */
double GetInnerRestMass(tGrid *grid, int bi);


/* funtions from mathematica */
void DNS_CTS(tVarList *vlFu, tVarList *vlu, tVarList *vlJdu, 
             tVarList *vldu, tVarList *vlduDerivs, int nonlin);
void DNS_compute_new_q(tGrid *grid, int iq);
double DNS_compute_new_q_atXYZ(tGrid *grid, int bi,
                               double X, double Y, double Z);
void DNS_set_restmassintegrand(tGrid *grid, int iInteg);
void DNS_set_J_ADM_VolInt_integrand(tGrid *grid, int iIntegx, int iIntegy, int iIntegz);
void DNS_set_M_ADM_VolInt_integrand(tGrid *grid, int iInteg);
void DNS_set_P_ADM_VolInt_integrand_Om_xcm(tGrid *grid,
                                           int iIntegx, int iIntegy, int iIntegz,
                                           double Om, double xcm);
void DNS_set_P_ADM_VolInt_integrand(tGrid *grid, int iIntegx, int iIntegy, int iIntegz);
void set_DNSdata_Sigma_BC(tVarList *vlFu, tVarList *vlu,
                          tVarList *vlJdu, tVarList *vldu,
                          tVarList *vlduDerivs, int nonlin);
void DNS_set_dlnIntegEuler(tGrid *grid, int ilnIntegEuler,
                           int idlnIntegEuler, double Om, double xcm);

/* for solving all ell. eqns together */
void F_DNSdata(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2);
void J_DNSdata(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu);

/* for solving the ell. eqns sequentially */
void F_oneComp(tVarList *vlFw, tVarList *vlw,
               tVarList *vlwDerivs, tVarList *vlc2);
void J_oneComp(tVarList *vlJdw, tVarList *vldw,
               tVarList *vldwDerivs, tVarList *vlw);

/* funcs from TOV */
int TOV_init(double Pc, int pr, double *rf_surf,
             double *m, double *Phi_c, double *Psi_c, double *m0);
int TOV_m_P_Phi_Psi_m0_OF_rf(double rf, double rf_surf,
                             double Pc, double Phic, double Psic,
                             double *m, double *P, double *Phi, double *Psi,
                             double *m0);

/* from DNSgrid.c */
int set_DNS_boxsizes(tGrid *grid);
int DNSdata_setup_boxes(tGrid *grid);
int set_DNS_box_attribs(tGrid *grid);
int set_sigma_pm_vars(tGrid *grid);
void reset_Coordinates_CubedSphere_sigma01(tGrid *grid, tGrid *gridnew,
                                           int star);
double InnerVolumeIntegral(tGrid *grid, int star, int vind);
double GridVolumeIntegral(tGrid *grid, int vind);
double VolumeIntegral_inDNSgridBox(tGrid *grid, int b, int vind);
tGrid *make_grid_with_sigma_pm(tGrid *grid, int nAB, int nphi, int nxyz);
int DNSgrid_Get_BoxAndCoords_of_xyz(tGrid *grid1,
                                    double *X1, double *Y1, double *Z1,
                                    tBox *box, int ind, double x, double y, double z);
int DNSgrid_set_bfaces(tGrid *grid, int set_fpts, int pr);
void Interpolate_Var_From_Grid1_To_Grid2(tGrid *grid1, tGrid *grid2, int vind);
void Interp_Var_From_Grid1_To_Grid2_star(tGrid *grid1, tGrid *grid2, int vind,
                                         int star);
void copy_Var_Box1ATlam1_to_Box2ATlam0(tGrid *grid, int vind, int b1, int b2);
void DNSgrid_init_Coords_for_star(tGrid *grid, int star);
void DNSgrid_init_Coords(tGrid *grid);
void DNSgrid_scale_Coordinates_CubSph_sigma(tGrid *grid, double fac, int star);
//void DNSgrid_copy_DomainShape(tGrid *grid, int ibd);
//void DNSgrid_set_Var_equalmasses_sym(tGrid *grid, int ibd, int iv, int sym);
void DNS_set_wB(tGrid *grid, int star, double xc,double yc,double zc);
void DNSgrid_load_initial_guess_from_checkpoint(tGrid *grid, char *filename);
void set_Var_to_Val_if_below_limit_or_outside(tGrid *grid, int vi, 
                                              double Val, double lim);
void set_Var_to_Val_atSurface(tGrid *grid, int vi, double Val);
int set_DNSdata_CoordFac(tGrid *grid);
//void debug_Coordinates_AnsorgNS_sigma_pm_B01(tGrid *grid, int dom, char *label);
//void debug_Coordinates_AnsorgNS_sigma_pm_B01_dom03(tGrid *grid, char *label);

/* from DNS_Interpolate_ADMvars.c */
int DNS_Interpolate_ADMvars(tGrid *grid);

/* from DNS_BCs.c */
void set_DNSdata_BCs(tVarList *vlFu, tVarList *vlu, tVarList *vluDerivs, int nonlin);

/* DNS_compute_chi */
int DNS_compute_chi(tGrid *grid);

/* from DNS_EoS.c */
int EoS_pwp; /* switch to pwp*/

void DNS_polytrope_EoS_of_hm1(double hm1,
                              double *rho0, double *P, double *rhoE,
                              double *drho0dhm1);
double DNS_polytrope_rho0_of_hm1(double hm1);
double DNS_polytrope_P_of_hm1(double hm1);
double DNS_polytrope_hm1_of_P(double P);
void DNS_polytrope_rho0_rhoE_of_P(double P, double *rho0, double *rhoE);
int DNS_pwp_init_file();
int DNS_pwp_init_parameter();
int DNS_poly_init();
void find_n_kappa(double P, double *n, double *kappa);
