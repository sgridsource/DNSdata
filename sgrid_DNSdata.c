/* sgrid_DNSdata.c */
/* Wolfgang Tichy, Jan 2018 */

#include "sgrid.h"
#include "DNSdata.h"


int sgrid_DNSdata() 
{
  if (!Getv("physics", "DNSdata")) return 0;
  printf("Adding DNSdata\n");

  /* functions */
  AddFun(PRE_GRID, set_DNS_boxsizes, "setup initial box sizes");
  AddFun(PRE_COORDINATES, DNSdata_setup_boxes, "setup boxes e.g. place Cub. Sph.");
  AddFun(PRE_INITIALDATA, set_DNS_box_attribs, "set box Attribs for DNSdata");
  AddFun(PRE_INITIALDATA, DNSdata_startup, "initialize DNSdata");
  AddFun(INITIALDATA, setDNSdata, "set the DNS data");
  AddFun(ANALYZE, DNSdata_analyze, "compute properties of DNS data");
  AddFun(ANALYZE, DNS_compute_chi, "Compute chi");
  AddFun(OUTPUT, DNS_Interpolate_ADMvars,
         "interpolate ADM initial data using DNSdata_Interpolate_pointsfile");

  /* variables */
  AddVar("DNSdata_Psi",     "",     "new conf. factor");
  AddVar("DNSdata_Psi",     "i",    "1st deriv of Psi");
  AddVar("DNSdata_Psi",     "(ij)", "2nd deriv of Psi");
  AddVar("DNSdata_B",       "I",    "shift in inertial coords");
  AddVar("DNSdata_B",       "Ij",   "1st deriv of B^i");
  AddVar("DNSdata_B",       "I(jk)","2nd deriv of B^i");
  AddVar("DNSdata_alphaP",  "",     "lapse times Psi");
  AddVar("DNSdata_alphaP",  "i",    "1st deriv of alphaP");
  AddVar("DNSdata_alphaP",  "(ij)", "2nd deriv of alphaP");
  AddVar("DNSdata_Sigma",   "",     "Sigma is potential for irrot. part of h*u^i");
  AddVar("DNSdata_Sigma",   "i",    "1st deriv of Sigma");
  AddVar("DNSdata_Sigma",   "(ij)", "2nd deriv of Sigma");
  AddVar("DNSdata_SigmaX",    "",     "1st A-deriv of Sigma");
  AddVar("DNSdata_SigmaXX",   "",     "2nd A-deriv of Sigma");
  AddVar("DNSdata_SigmaXXX",  "",     "3rd A-deriv of Sigma");
  AddVar("DNSdata_lSigmaX",   "",     "1st A-deriv of Sigma");
  AddVar("DNSdata_lSigmaXX",  "",     "2nd A-deriv of linearized Sigma");
  AddVar("DNSdata_lSigmaXXX", "",     "3rd A-deriv of linearized Sigma");

  AddVar("DNSdata_VR",  "I", "velocity in rotating frame: "
                             "V^i = u^i/u^0 - xi^i,  xi^i = KV");
  AddVar("DNSdata_wB",  "I", "conformal rotational part of h*u^i: "
                             "h u^i = w^i + D^i Sigma - h u^0 beta^i, "
                             "w^i = Psi^{-6} wB^i, D^i Sigma = g^{ij} d_i Sigma. "
                             "Note: v^i = u^i/u^0, i.e. it is not u^i/(alpha u^0)");
  AddVar("DNSdata_q",   "",  "q := h-1");
  AddVar("DNSdata_wB",  "Ij","1st deriv of wB");
  AddVar("DNSdata_q",   "i", "1st deriv of q");

  /* q we get from DNS_compute_new_q, with neg. values kept in boxes around stars */
  AddVar("DNSdata_qg", "", "global smooth q with neg. values kept");
  /* sometimes we save the old vars before the ell. solve */
  AddVar("DNSdata_Psiold",    "",  "old Psi");
  AddVar("DNSdata_Bold",      "I", "old B");
  AddVar("DNSdata_alphaPold", "",  "old alphaP");
  AddVar("DNSdata_Sigmaold",  "",  "old Sigma");
  AddVar("DNSdata_qgold",     "",  "old qg");
  AddVar("DNSdata_qnocent",   "",  "q we get without centering");

  AddVar("DNSdata_temp1", "", "temporary variable(e.g. to store derivs)");
  AddVar("DNSdata_temp2", "", "temporary variable(e.g. to store derivs)");
  AddVar("DNSdata_temp3", "", "temporary variable(e.g. to store derivs)");
  AddVar("DNSdata_temp4", "", "temporary variable(e.g. to store derivs)");

  AddVar("DNSdata_CoordFac", "", "scale PDEs by a coord. dependent factor");

  /* parameters */
  AddPar("DNSdata_rotationstate1", "corotation",
         "rotation state of NS1 [corotation,VwApproximation,rotation");
  AddPar("DNSdata_rotationstate2", "corotation",
         "rotation state of NS2 [corotation,VwApproximation,rotation");
  AddPar("DNSdata_wB_factor", "Psi6",
         "wB = DNSdata_omega cross (r-r_c) * factor [1,Psi6,h,1/alpha]");
  AddPar("DNSdata_omegax1", "0", "x-comp of angular velocity of NS1");
  AddPar("DNSdata_omegay1", "0", "y-comp of angular velocity of NS1");
  AddPar("DNSdata_omegaz1", "0", "z-comp of angular velocity of NS1");
  AddPar("DNSdata_omegax2", "0", "x-comp of angular velocity of NS2");
  AddPar("DNSdata_omegay2", "0", "y-comp of angular velocity of NS2");
  AddPar("DNSdata_omegaz2", "0", "z-comp of angular velocity of NS2");
  AddPar("DNSdata_m01",   "0.141202", "rest mass of NS1");
  AddPar("DNSdata_m02",   "0.141202", "rest mass of NS2");
  AddPar("DNSdata_iterate_m0", "no", "whether we iterate rest masses [no,yes]");
  if(Getv("DNSdata_iterate_m0", "yes"))
  {
    AddPar("DNSdata_init_m01", "-1", "initial rest mass of NS1, if positive "
           "DNSdata_m01 will be changed during the iterations until we reach "
           "DNSdata_desired_m01");
    AddPar("DNSdata_init_m02", "-1", "initial rest mass of NS2, if positive "
           "DNSdata_m02 will be changed during the iterations until we reach "
           "DNSdata_desired_m02");
    AddPar("DNSdata_m0change", "0.1", "amount by which change DNSdata_m01/2 "
           "during iterations");
    AddPar("DNSdata_desired_m01", "-1", "desired rest mass1 if we iterate it");
    AddPar("DNSdata_desired_m02", "-1", "desired rest mass2 if we iterate it");
  }
  AddPar("DNSdata_qm1",   "-1", "max q of NS1 we want");
  AddPar("DNSdata_qm2",   "-1", "max q of NS2 we want");
  AddPar("DNSdata_qmax1", "0", "max q of NS1 along x-axis");
  AddPar("DNSdata_qmax2", "0", "max q of NS2 along x-axis");
  AddPar("DNSdata_xmax1", "0", "pos. of max q of NS1 along x-axis");
  AddPar("DNSdata_xmax2", "0", "pos. of max q of NS2 along x-axis");
  AddPar("DNSdata_find_position_of_qmax", "find_XYZ",
         "how we set X,Y,Z of qmax in this func [find_XYZ,initial_XYZ_guess]");
  AddPar("DNSdata_find_xyz_of_qmax", "globalmax", "how we set x,y,x of qmax "
         "in this func [globalmax,xmax_along_axis,max_along_axis]");
  AddPar("DNSdata_actual_xmax1", "0", "x-pos. of actual global max of q in NS1");
  AddPar("DNSdata_actual_ymax1", "0", "y-pos. of actual global max of q in NS1");
  AddPar("DNSdata_actual_zmax1", "0", "z-pos. of actual global max of q in NS1");
  AddPar("DNSdata_actual_xmax2", "0", "x-pos. of actual global max of q in NS2");
  AddPar("DNSdata_actual_ymax2", "0", "y-pos. of actual global max of q in NS2");
  AddPar("DNSdata_actual_zmax2", "0", "z-pos. of actual global max of q in NS2");
  AddPar("DNSdata_q_derivs", "dq", "how we compute the derivs of q [dq,dqg]");
  AddPar("DNSdata_drho0_inBC", "dq", "what we use for drho0 in BC [dq,dA]");
  AddPar("DNSdata_new_q", "FromFields", "how we compute the new q "
         "[FromFields,Fromqgold]");
  AddPar("DNSdata_center_new_q", "no",
         "if and how we center new q on (DNSdata_xmax1/2,0,0) "
         "[no,center_yz,center_xyz,adjust_domainshapes,"
          "set_DNSdata_actual_xyzmax]");
  AddPar("DNSdata_center_new_q_timebin", "before_ell_solve", 
         "when we center q "
         "[before_ell_solve,after_adjusting_Omega_xCM]");
  AddPar("DNSdata_center_new_q_first_at", "0", "first iteration when we "
         "center q [#,-1]. -1 means never");
  AddPar("DNSdata_center_new_q_flag", "no",
         "set automatically in main iteration [no,yes]");
  AddPar("DNSdata_center_new_q_fac", "0.1",
         "by how much we try to move the qmax location");
  AddPar("DNSdata_center_new_q_maxdx", "0",
         "max dx allowed before we try to center");
  AddPar("DNSdata_center_new_q_maxdy", "0",
         "max dy allowed before we try to center");
  AddPar("DNSdata_center_new_q_maxdz", "0",
         "max dz allowed before we try to center");
//  AddPar("DNSdata_center_fields", "no",
//         "if and how we center fields on (DNSdata_xmax1/2,0,0) "
//         "[no,center_yz,center_xyz,reset_q,adjust_domainshapes]");
// there were more center fields options in BNSdata ...
  AddPar("DNSdata_Omega_init", "DNSdata_Omega", "ini. orbital angular velocity "
         "[#,estimate,estimate_from_desired_m0,DNSdata_Omega]");
  AddPar("DNSdata_x_CM_init", "DNSdata_x_CM",
         "initial center of mass in x-direction [#,estimate,DNSdata_x_CM]");
  AddPar("DNSdata_Omega", "estimate",  "orbital angular velocity "
         "[#,estimate,estimate_from_desired_m0]");
  AddPar("DNSdata_b",     "1",  "separation parameter (distance=2b)");
  AddPar("DNSdata_n",     "1",  "polytropic index n, Gamma = 1 + 1/n");
  AddPar("DNSdata_kappa", "1",  "kappa in EOS: P = kappa rho0^Gamma");
  AddPar("DNSdata_x_CM",  "estimate",  "center of mass in x-direction");
  AddPar("DNSdata_C1",    "-1", "C1 in q = (C1/F-1)/(n+1) "
         "[needs to be adjusted so that m01 stays the constant]");
  AddPar("DNSdata_C2",    "-1", "C2 in q = (C2/F-1)/(n+1) "
         "[needs to be adjusted so that m02 stays the constant]");
  AddPar("DNSdata_guess", "TOV", "init. guess for Psi, alphaP and "
         "q [initialize_from_checkpoint,TOV,TOVaverage,TOVproduct] "
         "for shift [TaniguchiShift]");
  AddPar("DNSdata_initfile", "", "name of first initialization file to read");
  AddPar("DNSdata_init_q_fromfields", "no", "init q from other fields [no,yes]");
  if(Getd("DNSdata_m02")==0.0)
  {
    AddPar("DNSdata_yshift1", "0", "shift NS1 in y-direction for testing");
    AddPar("DNSdata_adjustdomain01", "yes", "if we adjust domainshapes "
           "after shift [yes,no]");
  }
  AddPar("DNSdata_Sigma_surface_BCs", "AddInnerVolIntToBC",
         "BCs for Sigma on star surfaces "
         "[AddInnerVolIntToBC,AddInnerSumToBC,FakeMatterOutside,"
         "ZeroAtPoint,InnerVolIntZero,InnerSumZero,"
         "none,ZeroInOuterBoxes,EllEqn]");
  AddPar("DNSdata_InnerToOuterSigmaTransition", "C2", 
         "smoothness (at A=0) as we go from inner to outer box [C2,C1,C0]");
  AddPar("DNSdata_set_desired_VolAvSigmas", "no",
         "whether we set desired_VolAvSigma1/2 before ell. solves [no,yes]");
  AddPar("DNSdata_desired_VolAvSigma1", "0", "desired value of VolAvSigma "
         "(i.e. InnerVolInt or InnerSum) in Sigma surface BC for star1");
  AddPar("DNSdata_desired_VolAvSigma2", "0", "desired value of VolAvSigma "
         "(i.e. InnerVolInt or InnerSum) in Sigma surface BC for star2");
  AddPar("DNSdata_KeepInnerSigma", "no", "keep Sigma in inner boxes [no,yes]");
  AddPar("DNSdata_itmax", "10", "max. number of iterations in DNSdata_solve");
  AddPar("DNSdata_break_if_err_below_tol", "at_iterationend after_ell_solve",
         "which if clauses we use to break out of main iteration loop "
         "[at_iterationend,after_ell_solve]");
  AddPar("DNSdata_tol",   "1e-6", "tolerance for DNSdata_solve");
  AddPar("DNSdata_Newton_tolFac", "0.1", "tol for Newton is "
         "Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac)");
  AddPar("DNSdata_esw",   "1", "ell. solve weight: after ell. solve new "
         "and old values are averaged according to: "
         "new = esw*new + (1-esw)*old");
  AddPar("DNSdata_esw1",  "1", "second weight esw1 is used if better");
  AddPar("DNSdata_allow_esw1_first_at", "-1", "first iteration when esw=esw1 "
         "will be tried and allowed if better [#,-1]. -1 means never");
  AddPar("DNSdata_Sigma_esw", "0.2", "ell. solve weight for Sigma");
  AddPar("DNSdata_Sigma_esw1",  "1", "second Sigma weight is used if better");
  AddPar("DNSdata_allow_Sigma_esw1_first_at", "-1",
         "first iteration when Sigma_esw=Sigma_esw1 "
         "will be tried and allowed if better [#,-1]. -1 means never");
  AddPar("DNSdata_reset_qmax_xmax_pars_at", "-1", "list of its when we "
         "reset DNSdata_qmax1/2, DNSdata_xmax1/2");
  AddPar("DNSdata_analyze_xmax", "set_DNSdata_xmax",
         "what we do with DNSdata_xmax1/2 inside DNSdata_analyze "
         "[set_DNSdata_xmax,print_xmax]");
  AddPar("DNSdata_adjust_C1C2", "refineguess", "how to adjust C1/2 "
         "[refineguess,noguess,no]");
  AddPar("DNSdata_adjust", "nothing", "what we adjust (apart from C1/2) "
         "after ell. solve. E.g. \"forcebalance always\" adjusts Omega "
         "and x_CM to keep xmax1 and xmax2 in place. "
         "[nothing,forcebalance [always],"
         "Py0_forcebalance"
         "]");
  AddPar("DNSdata_adjust_first_at", "0", 
         "first iteration when we use DNSdata_adjust. -1 means never");
  AddPar("DNSdata_adjust_mintol", "1e-10", "always use tol>=mintol in adjust");
  AddPar("DNSdata_dOmega_fac", "0.1", "dOmega = Omega*dOmega_Fac");
  AddPar("DNSdata_dx_CM_fac",  "0.1", "dx_CM = b*dx_CM_fac");
  AddPar("DNSdata_set_Surface_q", "no", "what to do with q at surface [no,zero]");
  AddPar("DNSdata_set_negative_q", "no", "what to do with q<0 [no,zero]");
  AddPar("DNSdata_EllSolver_method", "DNS_Eqn_Iterator",
         "how we solve for Psi,B^i,alphaP,Sigma "
         "[allatonce, DNS_Eqn_Iterator, DNS_ordered_Eqn_Iterator,"
         " DNS_ordered_Var_Eqn_Iterator]");
  AddPar("DNSdata_CTS_Eqs_Iteration_order", 
         "DNSdata_Psi DNSdata_Bx DNSdata_By DNSdata_Bz DNSdata_alphaP",
         "Order we use in function DNS_ordered_Eqn_Iterator. "
         "This par has to contain all we solve for!");
  AddPar("DNSdata_extraSigmaSolve_fac", "0",
         "solve for Sigma if res_old/res <= fac");
  AddPar("DNSdata_extraSigmaSolve_every", "10",
         "how often we do the extra Sigma solve");
  AddPar("DNSdata_SigmaSolve_tolFac", "0",
         "use solver for Sigma if norm(Sigma_Err) >= norm(Rest_Err) * tolFac");
  AddPar("DNSdata_FinalEllSolveVars", "", "vars for which do another elliptic "
         "solve after all other iterations are done. E.g. DNSdata_Psi");
  AddPar("DNSdata_linSolver", "UMFPACK", 
         "linear solver used [LAPACK,templates_GMRES,bicgstab,UMFPACK]");
  AddPar("DNSdata_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("DNSdata_linSolver_tolFac","0.1", "tol for linSolver is " 
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");
  AddPar("DNSdata_linSolver_tol","0", "tol for linSolver is "
         "max2((*normres)*linSolv_tolFac, linSolv_tol)");

  /* pick root finders and their options */
  AddPar("DNSdata_sigp_Bphi_FINDER_reset_Coordinates_AnsorgNS_sigma_pm",
         "newton_linesrch_itsP", "[newton_linesrch_itsP,zbrent_itsP]");
  AddPar("DNSdata_m0_error_VectorFuncP_grid0", "grid",
         "grid from which we interpolate vars when domains change "
         "[grid,grid_bak]");

  AddPar("DNSdata_grid", "36CS_2xyz", "what grid we use [36CS_2xyz]");

  /* some par that contains values of variables controlling the 
     "main iteration loop", all variables are saved in this on par */
  AddPar("DNSdata_Main_Iteration_Loop_State","", "set by sgrid, not by user!");

  /* Tim's chi indicators */
  AddPar("DNSdata_mass_shedding1", "1", "mass shedding indicator chi for star 1");
  AddPar("DNSdata_mass_shedding2", "1", "mass shedding indicator chi for star 2");

  /* par for ecc. and inspiral */
  AddPar("DNSdata_ecc", "0", "Eccentricity parameter e of orbits");
  AddPar("DNSdata_rdot", "0", "time deriv. of relative coordinate r=|r1-r2|");
  AddPar("DNSdata_ADMshift", "B+phidotphi^i+rdotor0r^i", "what we put in ADM "
         "var beta at the very end in setADMvars "
         "[B^i+xi^i,B^i+phidotphi^i+rdotor0r^i]");

  /* pars for piecewise polytropes */
  AddPar("DNSdata_EoS_type", "poly", "poly. or piecewise poly. [poly,pwp]"); 
  if(Getv("DNSdata_EoS_type", "pwp"))
  {
     AddPar("DNSdata_pwp_rho0", "0.00023695302 0.00081212262 0.00162039766",
            "rho0 were we use transition between pieces");
  }

  /* par to improve condition number of PDE matrix */
  AddPar("DNSdata_CoordFacPower", "0", "scale PDEs by a coord. "
         "dependent factor raised to this Power [-2,-1,0,1,2]");

  return 0;
}
