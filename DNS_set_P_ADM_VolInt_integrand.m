(* DNS_set_P_ADM_VolInt_integrand.m 
   Wolfgang Tichy  11/2008       *)

(* set integrand for inner volume integrals of J_ADM components *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dSigma[a], q, wB[a], x,y,z,
            Integx,Integy,Integz}


constvariables = {OmegaCrossR[a], OmCrossR[a], xrdotor[a], xrd[a], 
                  omegMOmeg[a], xMxmax[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_S(box, Ind(\"DNSdata_Sigma\"),
					Ind(\"DNSdata_Sigmax\"));",

  (* which star are we considering? *)
  Cinstruction == "VwApprox = corot = 0;",
  Cif == (isSTAR1),
    Cinstruction == "if(corot1)    corot = 1;",
    Cinstruction == "if(VwApprox1) VwApprox = 1;",
    xmax == xmax1,
    omegMOmeg1 == omegax1,
    omegMOmeg2 == omegay1,
    omegMOmeg3 == omegaz1 - Omega,

  Cif == else,
    Cinstruction == "if(corot2)    corot = 1;",
    Cinstruction == "if(VwApprox2) VwApprox = 1;",
    xmax == xmax2,
    omegMOmeg1 == omegax2,
    omegMOmeg2 == omegay2,
    omegMOmeg3 == omegaz2 - Omega,
  Cif == end,

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",
  Cinstruction == "double rho0, P, rhoE, drho0dhm1; /* declare them */",

  (* set lapse and Psi4 *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,

  (* center of circle *)
  xC == xCM + ecc (xmax - xCM),
  
  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xC),
  OmegaCrossR3 == 0,

  (* Om \times r term, depends on Om,xcm args of this function *)
  xc == xcm + ecc (xmax - xcm),
  OmCrossR1 == - Om y,
  OmCrossR2 == + Om (x-xc),
  OmCrossR3 == 0,

  (* radial x^i rdot / r term *)
  rdotor == rdot / (xmax1-xmax2),
  xrdotor1 == (x - xCM) rdotor,
  xrdotor2 == y rdotor,
  xrdotor3 == z rdotor,

  (* radial x^i rdot / r term with xcm arg of this function *)
  rdotor == rdot / (xmax1-xmax2),
  xrd1 == (x - xcm) rdotor,
  xrd2 == y rdotor,
  xrd3 == z rdotor,

  (* the shift in the inertial frame is B^i *)
  (* the beta[a] used here is equal to \beta^i + \xi^i *)
  (* \beta^i + \xi^i has the same value in both inertial and comoving frames *)
  (* in an inertial frame xi[a] = OmegaCrossR[a] + xrdotor[a] *)
  (* in a comoving frame where \xi^i = 0, beta[a] is the shift in this frame *)
  beta[a] == B[a] + OmegaCrossR[a] + xrdotor[a],
  bet[a]  == B[a] + OmCrossR[a] + xrd[a],

  (**************)
  (* corotation *)
  (**************)
  Cif == ( corot ),

    (* check if we approximate V^i and w^i and set h from them and xi^i *)
    Cif == (VwApprox),
      (* V^i = (omega-Omega) /times r *)
      xMxmax1 == x - xmax,
      xMxmax2 == y,
      xMxmax3 == z,
      vR[a] == epsmatrix3d[a,b,c] omegMOmeg[b] xMxmax[c],
    (* pure corotation with V^i=0 and w^i=0 *)
    Cif == else,
      (* vR[a] is zero for corotation *)
      vR[a] == 0,
    Cif == end,

    (* compute square of u^0 in rotating frame *)
    oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
    uzerosqr == 1.0/oouzerosqr,

  (****************)
  (* general case *)
  (****************)
  Cif == else,
    Psim2 == 1/Psi2,
    Psim4 == Psim2*Psim2,
    Psim6 == Psim4*Psim2,
    DSigmaUp[a] == Psim4 dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],
    h == q + 1,
    h2 == h*h,
    uzerosqr == (1 + (wDown[a] + dSigma[a]) (w[a] + DSigmaUp[a])/h2)/alpha2,
    uzero == sqrt[uzerosqr],
    vR[a] == (w[a] + DSigmaUp[a])/(uzero*h) - beta[a],
  Cif == end,

  (* rest mass density, pressure, and total energy density *)
  Cinstruction == "EoS->vars_from_hm1(q[ijk], &rho0, &P, &rhoE, &drho0dhm1);",

  (* if q=0 all matter vars are zero, which can be enforced by uzerosqr=0 *)
  Cif == (q==0),
    uzerosqr == 0,
  Cif == end,

  (* get fluid var j in 3+1 *)
  jup[a] == alpha (rhoE + P) uzerosqr (vR[a]+bet[a]),

 (* Integrand for rest mass with 3-volume element factor Psi^6 *)
  Integx == Psi4 (  jup1 ) * Psi4*Psi2,
  Integy == Psi4 (  jup2 ) * Psi4*Psi2,
  Integz == Psi4 (  jup3 ) * Psi4*Psi2,

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "DNS_set_P_ADM_VolInt_integrand.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"DNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n"];
  pr["extern tEoS EoS[1];"];
  pr["\n\n\n"];

  pr["void DNS_set_P_ADM_VolInt_integrand_Om_xcm(tGrid *grid,
                                                 int iIntegx, int iIntegy, int iIntegz,
                                                 double Om, double xcm)\n"];
  pr["{\n"];

  pr["int VwApprox1 = Getv(\"DNSdata_rotationstate1\",\"VwApproximation\");\n"];
  pr["int VwApprox2 = Getv(\"DNSdata_rotationstate2\",\"VwApproximation\");\n"];
  pr["int corot1 = VwApprox1 || Getv(\"DNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = VwApprox2 || Getv(\"DNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int VwApprox, corot;\n"];
  pr["double Omega = Getd(\"DNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"DNSdata_x_CM\");\n"];
  pr["double ecc = Getd(\"DNSdata_ecc\");\n"];  
  pr["double rdot = Getd(\"DNSdata_rdot\");\n"];  
  pr["double xmax1 = Getd(\"DNSdata_actual_xmax1\");\n"];
  pr["double xmax2 = Getd(\"DNSdata_actual_xmax2\");\n"];
  pr["double omegax1 = Getd(\"DNSdata_omegax1\");\n"];
  pr["double omegay1 = Getd(\"DNSdata_omegay1\");\n"];
  pr["double omegaz1 = Getd(\"DNSdata_omegaz1\");\n"];
  pr["double omegax2 = Getd(\"DNSdata_omegax2\");\n"];
  pr["double omegay2 = Getd(\"DNSdata_omegay2\");\n"];
  pr["double omegaz2 = Getd(\"DNSdata_omegaz2\");\n"];
  pr["\n"];

  pr["int bi;\n"];
  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["int isSTAR1 = (box->SIDE == STAR1);\n\n"];
  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];
  prdecvarname[{z},      "z"];

  prdecvarname[{Psi},     "DNSdata_Psi"];
  prdecvarname[{B[a]},    "DNSdata_Bx"];
  prdecvarname[{alphaP},  "DNSdata_alphaP"];
  prdecvarname[{Sigma},   "DNSdata_Sigma"];
  prdecvarname[{q},      "DNSdata_q"];
  prdecvarname[{wB[a]}, "DNSdata_wBx"];
  prdecvarname[{dSigma[a]},     "DNSdata_Sigmax"];

  prdecvar[{Integx}, "iIntegx"];
  prdecvar[{Integy}, "iIntegy"];
  prdecvar[{Integz}, "iIntegz"];

  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["} /* end of boxes */\n"];
  pr["\n\n"];
  pr["}  /* end of function */\n\n"];

  (* add a second function *)
  pr["void DNS_set_P_ADM_VolInt_integrand(tGrid *grid, int iIntegx, int iIntegy, int iIntegz)\n"];
  pr["{\n"];
  pr["double Omega = Getd(\"DNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"DNSdata_x_CM\");\n"];
  pr["DNS_set_P_ADM_VolInt_integrand_Om_xcm(grid, iIntegx, iIntegy, iIntegz ,Omega, xCM);\n"];
  pr["}  /* end of 2nd function */\n\n"];

];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"
