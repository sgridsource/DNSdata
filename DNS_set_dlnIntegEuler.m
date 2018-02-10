(* DNS_set _dlnIntegEuler.m 
   Wolfgang Tichy  11/2011    *)

(* compute the updated q *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dSigma[a], q, wB[a], x, y, z,
             lnIntegEuler, dlnIntegEuler[a]}

constvariables = {OmegaCrossR[a], OmCrossR[a], xrdotor[a], xrd[a],
                 omegMOmeg[a], xMxmax[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_S(box,index_DNSdata_Sigma,
                                   Ind(\"DNSdata_Sigmax\"));",
  Cif == (bi==1),
    Cinstruction == "\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmax\"), 0,1);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmay\"), 0,1);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmaz\"), 0,1);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBx\"), 0,1);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBy\"), 0,1);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBz\"), 0,1);",
  Cif == end,
  Cif == (bi==2),
    Cinstruction == "\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmax\"), 3,2);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmay\"), 3,2);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmaz\"), 3,2);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBx\"), 3,2);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBy\"), 3,2);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBz\"), 3,2);",
  Cif == end,

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

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,
  Psim4 == 1/Psi4,
  Psim2 == Psim4 Psi2,

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

      (* compute u^0 in rotating frame *)
      betVRDown[a] == Psi4 delta[a,b] (bet[b] + vR[b]),
      oouzerosqr == alpha2 - betVRDown[a] (bet[a] + vR[a]),
      Cif == (oouzerosqr==0),
        oouzerosqr == -1,
      Cif == end,
      Cif == (oouzerosqr<0),
        uzero == -1, (* -Sqrt[-1/oouzerosqr], *)
      Cif == else,
        uzero == Sqrt[1/oouzerosqr],
      Cif == end,

      (* more terms *)
      h == q + 1,
      wBDown[a] == wB[a],
      wDown[a] == Psim2 wBDown[a],

      (* Integ. Euler is:
         h (1/uzero + uzero vR[a] betVRDown[a]) - vR[a] wDown[a] = -C *)
      (* set log of integrated Euler *)
      lnIntegEuler == 
              log[ h/uzero + h uzero vR[a] betVRDown[a] - vR[a] wDown[a] ],

    (* pure corotation with V^i=0 and w^i=0 *)
    Cif == else,

      (* vR[a] is zero for corotation *)
      (* vR[a] == 0,  ==>  Integ. Euler is h/uzero = - C  *)

      (* compute u^0 in rotating frame *)
      oouzerosqr == alpha2 - Psi4 delta[b,c] (bet[b]) (bet[c]),
      Cif == (oouzerosqr==0),
        oouzerosqr == 1,
      Cif == end,

      (* set log of integrated Euler *)
      lnIntegEuler == log[oouzerosqr],

    Cif == end,

  (****************)
  (* general case *)
  (****************)
  Cif == else,

    Psim6 == Psim4 Psim2,
    DSigmaUp[a] == Psim4 dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],

    h == q + 1,
    h2 == h h, 
    L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
    uzerosqr == L2/(alpha2 h2),
    uzero == sqrt[uzerosqr],

    (* killing vec xi^i in rotating frame, xi^a = (xi^0, xi^i) = (1,0,0,0) *)
    (* xi[a] == 0, *)
    (* below all xi[a] have been removed because beta[a] is already
       \beta^i + \xi^i in comoving frame! *)

    (* U and U0 from Gorgoulhon, PRD 63, 064029 *)
    (* here I compute all those using Omega, xCM *)
    U0[a] == ( beta[a] + w[a]/(h uzero) )/alpha,
    U0Down[a] == Psi4 U0[a],
    U[a] == DSigmaUp[a]/(alpha h uzero),
    (* Gamma factors *)
    Gamman == alpha uzero,
    Gamma0 == 1/sqrt[1 - U0[a] U0Down[a]],
    Gamma  == Gamman Gamma0 * 
              ( 1 - U0Down[a] U[a] - wDown[a] w[a]/(alpha2 h2 uzerosqr) ),

    (* set log of integrated Euler *)
    (* Om, xcm enter only the next 2 terms *)
    betw[a] == bet[a] + w[a]/(h uzero),
    betwDown[a] == Psi4 betw[a],
    lnIntegEuler == log[alpha2 - betw[a] betwDown[a] ] + 2 log[Gamma],

  Cif == end,

  Cinstruction == "} /* end of points loop */\n",

  (* compute cart derivs of lnIntegEuler *)
  (* CAUTION: this will write into the vars
     idlnIntegEuler, idlnIntegEuler+1, idlnIntegEuler+2 *)
  Cinstruction == "FirstDerivsOf_S(box, ilnIntegEuler, idlnIntegEuler);"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "DNS_set_dlnIntegEuler.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"DNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define Abs(x)     (fabs((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void DNS_set_dlnIntegEuler(tGrid *grid, int ilnIntegEuler,
                                 int idlnIntegEuler, double Om, double xcm)\n"];
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

  prdecvar[{lnIntegEuler},  "ilnIntegEuler"];
  prdecvar[{dlnIntegEuler}, "idlnIntegEuler"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];
  prdecvarname[{z},      "z"];

  (* prdecvarname[{g[a,b]},    "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{q},         "DNSdata_q"];
  prdecvarname[{Psi},       "DNSdata_Psi"];
  prdecvarname[{alphaP},    "DNSdata_alphaP"];
  prdecvarname[{B[a]},      "DNSdata_Bx"];
  prdecvarname[{wB[a]},     "DNSdata_wBx"];
  prdecvarname[{Sigma},     "DNSdata_Sigma"];
  prdecvarname[{dSigma[a]}, "DNSdata_Sigmax"];


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
