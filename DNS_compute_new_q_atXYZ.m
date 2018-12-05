(* DNS_compute _new _q _atXYZ.m 
   Wolfgang Tichy  2/2008       *)

(* compute the updated q *)


(* variables *)
variables = {DNSPsi, DNSB[a], DNSalphaP, DNSSigma, DNSdSigma[a], DNSwB[a],
             DNSqgold, temp4}

constvariables = {Psi, B[a], alphaP, Sigma, dSigma[a], wB[a], x,y, qgold,
                  OmegaCrossR[a], xrdotor[a], omegMOmeg[a], xMxmax[a]} 

(* compute in this order *)
tocompute = {

  (* deriv of Sigma *)
  Cinstruction == "FirstDerivsOf_S(box,index_DNSdata_Sigma,
                                   Ind(\"DNSdata_Sigmax\"));",
  Cif == (MATTRtouch && (!wB0outside)),
    Cinstruction == "int biin = bi-6; /* works only for my CubSph setup */ \n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmax\"), biin,bi);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmay\"), biin,bi);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_Sigmaz\"), biin,bi);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBx\"), biin,bi);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBy\"), biin,bi);\n
    copy_Var_Box1ATlam1_to_Box2ATlam0(grid, Ind(\"DNSdata_wBz\"), biin,bi);",
  Cif == end,

  (* x,y,z *)
  Cinstruction == "if(box->x_of_X[1] != NULL) {",
  Cinstruction == "x = box->x_of_X[1]((void *) box, -1, X,Y,Z);",
  Cinstruction == "y = box->x_of_X[2]((void *) box, -1, X,Y,Z);",
  Cinstruction == "z = box->x_of_X[3]((void *) box, -1, X,Y,Z);",
  Cinstruction == "} else {",
  Cinstruction == "x = X;",
  Cinstruction == "y = Y;",
  Cinstruction == "z = Z;",
  Cinstruction == "}",

  (* Psi, B[a], alphaP, Sigma, dSigma[a], wB[a] by interpolation *)
  Cinstruction == "spec_Coeffs(box, DNSPsi, temp4);",
  Cinstruction == "Psi = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSB1, temp4);",
  Cinstruction == "B1 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSB2, temp4);",
  Cinstruction == "B2 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSB3, temp4);",
  Cinstruction == "B3 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSalphaP, temp4);",
  Cinstruction == "alphaP = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSSigma, temp4);",
  Cinstruction == "Sigma = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSdSigma1, temp4);",
  Cinstruction == "dSigma1 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSdSigma2, temp4);",
  Cinstruction == "dSigma2 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSdSigma3, temp4);",
  Cinstruction == "dSigma3 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSwB1, temp4);",
  Cinstruction == "wB1 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSwB2, temp4);",
  Cinstruction == "wB2 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, DNSwB3, temp4);",
  Cinstruction == "wB3 = spec_interpolate(box, temp4, X,Y,Z);",

   (* which star are we considering? *)
  Cinstruction == "VwApprox = corot = 0;",
  Cinstruction == "if(box->SIDE==STAR1) {",
    CC == C1,
    xmax == xmax1,
    Cinstruction == "if(corot1)    corot = 1;",
    Cinstruction == "if(VwApprox1) VwApprox = 1;",
    Cinstruction == "if(VwApprox) {
                       omegax1 = Getd(\"DNSdata_omegax1\");
                       omegay1 = Getd(\"DNSdata_omegay1\");
                       omegaz1 = Getd(\"DNSdata_omegaz1\"); }",
    omegMOmeg1 == omegax1,
    omegMOmeg2 == omegay1,
    omegMOmeg3 == omegaz1 - Omega,
  Cinstruction == "} else {",
    CC == C2,
    xmax == xmax2,
    Cinstruction == "if(corot2)    corot = 1;",
    Cinstruction == "if(VwApprox2) VwApprox = 1;",
    Cinstruction == "if(VwApprox) {
                       omegax2 = Getd(\"DNSdata_omegax2\");
                       omegay2 = Getd(\"DNSdata_omegay2\");
                       omegaz2 = Getd(\"DNSdata_omegaz2\"); }",
    omegMOmeg1 == omegax2,
    omegMOmeg2 == omegay2,
    omegMOmeg3 == omegaz2 - Omega,
  Cinstruction == "}",

  (**************)
  (* compute q: *)

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
      betaVRDown[a] == Psi4 delta[a,b] (beta[b] + vR[b]),
      oouzerosqr == alpha2 - betaVRDown[a] (beta[a] + vR[a]),
      Cif == (oouzerosqr==0),
        oouzerosqr == -1,
      Cif == end,
      Cif == (oouzerosqr<0),
        uzero == -1, (* -Sqrt[-1/oouzerosqr], *)
      Cif == else,
        uzero == Sqrt[1/oouzerosqr],
      Cif == end,

      (* more terms *)
      wBDown[a] == wB[a],
      wDown[a] == Psim2 wBDown[a],

      (* set q = h-1 *)
      h == (vR[a] wDown[a] - CC) uzero /(1.0 + vR[a] betaVRDown[a]/oouzerosqr),
      q == (h - 1.0),

    (* pure corotation with V^i=0 and w^i=0 *)
    Cif == else,

      (* vR[a] is zero for pure corotation *)
      vR[a] == 0,

      (* compute u^0 in rotating frame *)
      oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
      Cif == (oouzerosqr==0),
        oouzerosqr == -1,
      Cif == end,
      Cif == (oouzerosqr<0),
        uzero == -1, (* -Sqrt[-1/oouzerosqr], *)
      Cif == else,
        uzero == Sqrt[1/oouzerosqr],
      Cif == end,

      (* set q = h-1 *)
      q == -CC uzero - 1.0,

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
    Cif == qFromFields,
      twoalpha2wdSigmapw == 2 alpha2 w[c] (dSigma[c]+wDown[c]),
      betadSigmaMinusCC == beta[c] dSigma[c] - CC, 
      bb == betadSigmaMinusCC^2 + twoalpha2wdSigmapw,
      L2 == (bb + Sqrt[Abs[bb*bb + twoalpha2wdSigmapw^2]])/(2 alpha2),
      h == Sqrt[Abs[L2 - (dSigma[a]+wDown[a]) (DSigmaUp[a]+w[a])]],
    Cif == else,
      Cinstruction == "spec_Coeffs(box, DNSqgold, temp4);",
      Cinstruction == "qgold = spec_interpolate(box, temp4, X,Y,Z);",
      (* h == q + 1, *)
      hOLD == qgold + 1.0,
      hOLD2 == hOLD*hOLD,
      LOLD2 == hOLD2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
      uzerosqr == LOLD2/(alpha2 hOLD2),
      uzero == sqrt[uzerosqr],
      vR[a] == (w[a] + DSigmaUp[a])/(uzero*hOLD) - beta[a],
      h == -( CC + dSigma[a] vR[a] ) uzero,
    Cif == end,

    (* h == q + 1, *)
    q == (h - 1.0),
  Cif == end,

  Cinstruction == "/* end of computation */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "DNS_compute_new_q_atXYZ.c"

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

  pr["double DNS_compute_new_q_atXYZ(tGrid *grid, int bi, double X, double Y, double Z)\n"];
  pr["{\n"];

  pr["int VwApprox1 = Getv(\"DNSdata_rotationstate1\",\"VwApproximation\");\n"];
  pr["int VwApprox2 = Getv(\"DNSdata_rotationstate2\",\"VwApproximation\");\n"];
  pr["int corot1 = VwApprox1 || Getv(\"DNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = VwApprox2 || Getv(\"DNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int VwApprox, corot;\n"];
  pr["int qFromFields = Getv(\"DNSdata_new_q\",\"FromFields\");\n"];
  pr["int wB0outside  = Getv(\"DNSdata_wB_outside\",\"0\");\n"];
  pr["double C1 = Getd(\"DNSdata_C1\");\n"];
  pr["double C2 = Getd(\"DNSdata_C2\");\n"];
  pr["double Omega = Getd(\"DNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"DNSdata_x_CM\");\n"];
  pr["double ecc = Getd(\"DNSdata_ecc\");\n"];  
  pr["double rdot = Getd(\"DNSdata_rdot\");\n"];  
  pr["double xmax1 = Getd(\"DNSdata_actual_xmax1\");\n"];
  pr["double xmax2 = Getd(\"DNSdata_actual_xmax2\");\n"];
  pr["double omegax1,omegay1,omegaz1, omegax2,omegay2,omegaz2;\n"];
  pr["\n"];

  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["int MATTRtouch  = (box->MATTR== TOUCH);\n\n"];
  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  (* prdecvarname[{g[a,b]},    "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{DNSPsi},       "DNSdata_Psi"];
  prdecvarname[{DNSalphaP},    "DNSdata_alphaP"];
  prdecvarname[{DNSB[a]},      "DNSdata_Bx"];
  prdecvarname[{DNSq},         "DNSdata_q"];
  prdecvarname[{DNSwB[a]},     "DNSdata_wBx"];
  prdecvarname[{DNSSigma},     "DNSdata_Sigma"];
  prdecvarname[{DNSdSigma[a]}, "DNSdata_Sigmax"];
  prdecvarname[{DNSqgold},     "DNSdata_qgold"];
  prdecvarname[{temp4},        "DNSdata_temp4"];

  pr["double Psi, B1,B2,B3, alphaP;\n"];
  pr["double Sigma, dSigma1,dSigma2,dSigma3, wB1,wB2,wB3, x,y,z, qgold;\n"];

  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["return q;\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"
