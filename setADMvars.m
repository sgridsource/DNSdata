(* setADMvars.m 
   Wolfgang Tichy  12/2007       *)

(* compute ADMvars from the DNSdata vars *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dB[a,b], dSigma[a],
	     psi, g[a,b], alpha, beta[a], K[a,b], rho, jdo[a], Sdo[a,b],
             q, wB[a], VR[a], x, y, z, rhobar}

constvariables = {OmegaCrossR[a], xrdotor[a], omegMOmeg[a], xMxmax[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_Sa(box, Ind(\"DNSdata_Bx\"), 
					 Ind(\"DNSdata_Bxx\"));",
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

  (* radial x^i rdot / r term *)
  rdotor == rdot / (xmax1-xmax2),
  xrdotor1 == (x - xCM) rdotor,
  xrdotor2 == y rdotor,
  xrdotor3 == z rdotor,

  (* the shift in the inertial frame is B^i *)
  (* the beta[a] used here is equal to \beta^i + \xi^i *)
  (* \beta^i + \xi^i has the same value in both inertial and comoving frames *)
  (* in an inertial frame xi[a] = OmegaCrossR[a] + xrdotor[a] *)
  (* in a comoving frame where \xi^i = 0, beta[a] is the shift in this frame *)
  beta[a] == B[a] + OmegaCrossR[a] + xrdotor[a],

  (* get 1st derivs of B *)
  (* Note: if B^i = beta^i - (Omega \times r)^i 
	=> vecLapB = vecLapbeta , LB = Lbeta, 
	since the L of any Killingvec is zero *)
  gdB == delta[a,b] dB[a,b],
  LB[a,b] == dB[a,b] + dB[b,a] -(2/3) delta[a,b] gdB,
  LBdo[a,b] == delta[a,c] delta[b,d] LB[c,d], 

  (* set psi, g_ij and K_ij *)
  psi == 1, (* set ADMvars psi to one *)
  g[a,b] == Psi4 delta[a,b],
  K[a,b] == Psi4 LBdo[a,b] / (2 alpha),

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

  (* set VR on grid equal to local vR, so when setADMvars is called in the
     end DNSdata_VR is set correctly! *)
  VR[a] == vR[a],

  (* rest mass density, pressure, and total energy density *)
  Cinstruction == "DNS_polytrope_EoS_of_hm1(q[ijk], &rho0, &P, &rhoE, &drho0dhm1);",

  (* if q=0 all matter vars are zero, which can be enforced by uzerosqr=0 *)
  Cif == (q==0),
    uzerosqr == 0,
  Cif == end,

  (* set fluid vars in 3+1 *)
  rho  == alpha2 (rhoE + P) uzerosqr - P,
  jup[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),

  (* set ADMvars *)
  jdo[a] == g[a,b] jup[b],
  (* Sup[a,b] == (rhoE + P) uzerosqr (vR[a]+beta[a]) (vR[b]+beta[b]) +
                  P ginv[a,b], *)
  vRplusbetado[a] == g[a,b] (vR[b]+beta[b]),
  Sdo[a,b] == (rhoE + P) uzerosqr vRplusbetado[a] vRplusbetado[b] + P g[a,b],   

  (* do we reset shift beta? (because so far beta = \beta^i + \xi^i,
     which is discontinuous!) *)
  Cif == shiftver1,
    (* new Omega \times r term *)
    OmegaCrossR1 == - (1-ecc) Omega y,
    OmegaCrossR2 == + (1-ecc) Omega (x-xCM),
    OmegaCrossR3 == 0,
    (* use new Omega \times r in ADM shift *)
    beta[a] == B[a] + OmegaCrossR[a] + xrdotor[a],
  Cif == end,

  (* set rescaled rhobar *)
  rhobar == rho Psi4 Psi4,

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_]   :=  g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_]   :=  K[b,a] /; !OrderedQ[{a,b}]
Sdo[a_,b_] :=  Sdo[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Kup[a_,b_]  := Kup[b,a]  /; !OrderedQ[{a,b}]

ddPsi[a_,b_]     := ddPsi[b,a]    /; !OrderedQ[{a,b}]
ddB[a_,b_,c_]    := ddB[a,c,b]    /; !OrderedQ[{b,c}]
LB[a_,b_]        := LB[b,a]       /; !OrderedQ[{a,b}]
ddalphaP[a_,b_]  := ddalphaP[b,a] /; !OrderedQ[{a,b}]
ddSigma[a_,b_]   := ddSigma[b,a]  /; !OrderedQ[{a,b}]

ddlPsi[a_,b_]     := ddlPsi[b,a]    /; !OrderedQ[{a,b}]
ddlB[a_,b_,c_]    := ddlB[a,c,b]    /; !OrderedQ[{b,c}]
LlB[a_,b_]        := LlB[b,a]       /; !OrderedQ[{a,b}]
ddlalphaP[a_,b_]  := ddlalphaP[b,a] /; !OrderedQ[{a,b}]
ddlSigma[a_,b_]   := ddlSigma[b,a]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "setADMvars.c"

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

  pr["\n\n\n"];

  pr["void setADMvars(tGrid *grid)\n"];
  pr["{\n"];

  pr["int VwApprox1 = Getv(\"DNSdata_rotationstate1\",\"VwApproximation\");\n"];
  pr["int VwApprox2 = Getv(\"DNSdata_rotationstate2\",\"VwApproximation\");\n"];
  pr["int corot1 = VwApprox1 || Getv(\"DNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = VwApprox2 || Getv(\"DNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int VwApprox, corot;\n"];
  pr["int shiftver1 = Getv(\"DNSdata_ADMshift\",\"B^i+phidotphi^i+rdotor0r^i\");\n"];
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

  prdecvarname[{Psi},     "DNSdata_Psi"];
  prdecvarname[{B[a]},    "DNSdata_Bx"];
  prdecvarname[{alphaP},  "DNSdata_alphaP"];
  prdecvarname[{Sigma},   "DNSdata_Sigma"];
  prdecvarname[{q},       "DNSdata_q"];
  prdecvarname[{wB[a]},   "DNSdata_wBx"];
  prdecvarname[{VR[a]},   "DNSdata_VRx"];
  prdecvarname[{rhobar},  "DNSdata_rhobar"];

  prdecvarname[{dB[a,b]}, 	"DNSdata_Bxx"];
  prdecvarname[{dSigma[a]},     "DNSdata_Sigmax"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];
  prdecvarname[{z},      "z"];

  prdecvarname[{g[a,b]}, "gxx"];
  prdecvarname[{psi}, "psi"];
  prdecvarname[{K[a,b]}, "Kxx"];
  prdecvarname[{alpha},  "alpha"];
  prdecvarname[{beta[a]},"betax"];
  prdecvarname[{rho},    "rho"];
  prdecvarname[{jdo[a]},   "jx"];
  prdecvarname[{Sdo[a,b]}, "Sxx"];

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
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"
