(* DNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute residuals of DNS ham, mom, alphaP and Sigma eqns *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP ,FSigma,
              dPsi[a],   dB[a,b],   dalphaP[a],    dSigma[a],
             ddPsi[a,b],ddB[a,b,c],ddalphaP[a,b], ddSigma[a,b],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
	     g[a,b], alpha, beta[a], K[a,b], 
             q, wB[a], dq[a], dwB[a,b], VR[a], x, y, z, lam,dlam[a],
             dSigmadlam,dlSigmadlam, ddSigmadlam2,ddlSigmadlam2,
             dddSigmadlam3,dddlSigmadlam3, rhobar, CoordFac}

constvariables = {OmegaCrossR[a], xrdotor[a], omegMOmeg[a], xMxmax[a]}

(* compute in this order *)
tocompute = {

  (* get i,j,k and check if they are in current subbox *)
  Cinstruction == "if(blkinfo!=NULL) {
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int k = kOfInd_n1n2(ijk,n1,n2);
    int j = jOfInd_n1n2_k(ijk,n1,n2,k);
    int i = iOfInd_n1n2_jk(ijk,n1,n2,j,k);
    int sbi = blkinfo->sbi;
    int sbj = blkinfo->sbj;
    int sbk = blkinfo->sbk;
    int nsb1 = blkinfo->nsb1;
    int nsb2 = blkinfo->nsb2;
    int nsb3 = blkinfo->nsb3;
    int i1,i2, j1,j2, k1,k2;
    IndexRangesInSubbox(i1,i2, j1,j2, k1,k2, sbi,sbj,sbk, nsb1,nsb2,nsb3);
    if(i<i1 || i>=i2) continue;
    if(j<j1 || j>=j2) continue;
    if(k<k1 || k>=k2) continue;
  }",
 
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
  LBLB == LB[a,b] LBdo[a,b],

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi3   == Psi*Psi2,
  Psi4   == Psi2*Psi2,
  Psi5   == Psi*Psi4,
  Psim1  == 1/Psi,
  Psim2  == Psim1*Psim1,
  Psim3  == Psim2*Psim1,
  Psim4  == Psim2*Psim2,

  (* rest mass density rho0, pressure P, and total energy density rhoE *)
  Cinstruction == "EoS_T0->vars_from_hm1(q[ijk],&rho0, &P, &rhoE, &drho0dhm1);",

  (*****************************)
  (* BEGIN: corot/general case *)
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
    Psim8 == Psim4*Psim4,
    Psim6 == Psi2*Psim8,
    Psim5 == Psim6*Psi,
    Psim7 == Psim8*Psi,
    Psim9 == Psim8*Psim1,
    h == q + 1,
    h2 == h h,
    DSigmaUp[a] == Psim4 dSigma[a],
    dSigmaUp[a] == dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],
    L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
    uzerosqr == L2/(alpha2 h2),
    (* oouzerosqr == 1.0/uzerosqr, *)
    uzero == sqrt[uzerosqr],
    vR[a] == (w[a] + DSigmaUp[a])/(uzero*h) - beta[a],

    (* more terms which we need inside the stars *)
    drho0[a] == drho0dhm1 dq[a],
    dLnPsi[a] == dPsi[a]/Psi,
    dLnh[a] == dq[a] / h,
    dLnalphaP[a] == dalphaP[a]/alphaP,
    dalpha[a] == dalphaP[a]/Psi - alphaP dPsi[a]/Psi2,
    dLnalpha[a] == dalpha[a]/alpha,
    dL2[a] == 2*(Psim8 wBDown[c] dwB[c,a] +
                 Psim6 (dwB[c,a] dSigma[c] + wB[c] ddSigma[a,c]) +
                 Psim4 dSigmaUp[c] ddSigma[a,c])  - 
              (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
               4 Psim5 dSigma[c] dSigmaUp[c]) dPsi[a] + 2 h2 dLnh[a],
    duzerosqr[a] == (dL2[a] - 2 L2 (dalpha[a]/alpha + dLnh[a]))/(alpha2 h2),
    duzero[a] == duzerosqr[a]/(2 uzero),
    dbeta[a,b] == dB[a,b] + epsmatrix3d[b,a,3] Omega + delta[a,b] rdotor,
   
    (* dLnrhozalphaPsi2oh[a] == dLnrho0[a] + dLnalphaP[a] + dLnPsi[a] - dLnh[a],
       dLnrhozalphaPsi6uz[a] == dLnrho0[a] + dLnalphaP[a] + 5 dLnPsi[a] +
                                duzero[a]/uzero, *)
    dLnalphaPsi2oh[a] == dLnalphaP[a] + dLnPsi[a] - dLnh[a],
    dLnalphaoh[a]     == dLnalphaP[a] - dLnPsi[a] - dLnh[a],
    dLnalphaPsi6uz[a] == dLnalphaP[a] + 5 dLnPsi[a] + duzero[a]/uzero,
    drho0PLUSrho0dLnalphaPsi2oh[a] == drho0[a] + rho0 dLnalphaPsi2oh[a],
    drho0PLUSrho0dLnalphaoh[a]     == drho0[a] + rho0 dLnalphaoh[a],
    drho0PLUSrho0dLnalphaPsi6uz[a] == drho0[a] + rho0 dLnalphaPsi6uz[a],

    divwB == delta [b,c] dwB[b,c],
    divbeta == delta [b,c] dbeta[b,c],
  Cif == end,
  (* END: corot/general case *)
  (***************************)

  (* V^i =: VR[a] = vR[a] *)
  VR[a] == vR[a],  (* set VR on grid equal to local vR *)

  (* rest mass density rho0, pressure P, and total energy density rhoE 
     have been computed already *)

  (* fluid vars in 3+1 *)
  rho  == alpha2 (rhoE + P) uzerosqr - P,
  j[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),
  S    == 3P - rhoE + rho,

  (* dLnalphaPsim6[i] = \partial_i ln(alpha Psi^{-6}) 
			= \partial_i ln(alphaP Psi^{-7})
			= Psi^7 alphaP^{-1} \partial_i(alphaP Psi^{-7}) *)
  dLnalphaPsim6[a] == dalphaP[a]/alphaP - 7 dPsi[a]/Psi,

  (************************************************)
  (* decide if use non-linear or linear equations *)
  (************************************************)
  (* non-linear case: *)
  Cif == nonlin, (* non-linear case *)

    vecLapB[a] == delta[b,c] (ddB[a,b,c] + (1/3) ddB[b,c,a]),

    (* equations for Psi, B[a], alphaP, Sigma *)
    FPsi    == delta[b,c] ddPsi[b,c] + Psi5 LBLB/(32 alpha2) +
               2Pi ( Psi5  rho    (1 - CTSmod) +
                     Psim3 rhobar (CTSmod) ), 
    FB[a]   == vecLapB[a] - LB[a,b] dLnalphaPsim6[b] -
               16Pi alpha Psi4 j[a],
    FalphaP == delta[b,c] ddalphaP[b,c] - alphaP (
               (7/8) Psi4 LBLB/(4 alpha2) + 2Pi Psi4 (rho+2S) ),

    (*****************************)
    (* BEGIN: corot/general case *)
    (**************)
    (* corotation *)
    (**************)
    Cif == ( corot ),
      (* set Sigma to zero for corotation *)
      FSigma  == Sigma,
        (* Pedro's thing for inside stars (which I don't use anymore): 
           FSigma == delta[b,c] ddSigma[b,c] + 
                   (wB[a] + dSigma[a]) *
                   (dLnrho0[a] + dLnuzerosqr[a]/2 + dLnalphaP[a] + 5 dLnPsi[a]),
        *)
    (****************)
    (* general case *)
    (****************)
    Cif == else,
      Cinstruction == "if(MATTRinside) {", (* inside stars *)
        ddSigCoef == rho0 + rhMeps rho0max ((rho0max - rho0)/rho0max)^rhMpow,
        FSigma == ddSigCoef delta[b,c] ddSigma[b,c] + 
                  dSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (wB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divwB) -
                  h uzero Psi4 (rho0 divbeta +
                                beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]),
      Cinstruction == "} else if(MATTRtouch) {", (* touches star surface *)

        Cinstruction == "if(FakeMatterOutside) {", (* use fake matter *)
          Cinstruction == "if(LapSig) {",  (* Laplace Sigma *)
            FSigma == delta[b,c] ddSigma[b,c],
          Cinstruction == "} else {",
            (* use h=1, rho0 = -lam as fake matter *)
            hf == 1,       (* fake h *)
            Cif == FakeT0,
              rhof == -lam,  (* fake rho0=-lam *)
              drhof[a] == -dlam[a],
            Cif == else,
              rhof == -1,    (* fake rho0=-1 *)
              drhof[a] == 0,
            Cif == end,
            (* same terms as above but with dLnh = 0 *)
            dLnalphaPsi2[a] == dLnalphaP[a] + dLnPsi[a],
            dLnalpha[a]     == dLnalphaP[a] - dLnPsi[a],
            dLnalphaPsi6uz[a] == dLnalphaP[a] + 5 dLnPsi[a] + duzero[a]/uzero,
            drhofPLUSrhofdLnalphaPsi2[a] == drhof[a] + rhof dLnalphaPsi2[a],
            drhofPLUSrhofdLnalpha[a]     == drhof[a] + rhof dLnalpha[a],
            drhofPLUSrhofdLnalphaPsi6uz[a] == drhof[a] + rhof dLnalphaPsi6uz[a],
            (* FSigma with fake matter *)
            FSigma == rhof delta[b,c] ddSigma[b,c] + 
                      dSigmaUp[c] drhofPLUSrhofdLnalphaPsi2[c] +
                      Psim2 (wB[c] drhofPLUSrhofdLnalpha[c] + rhof divwB) -
                      hf uzero Psi4 (rhof divbeta +
                                  beta[c] drhofPLUSrhofdLnalphaPsi6uz[c]),

          Cinstruction == "}",
        Cinstruction == "} else {", (* continous Sigma *)
          FSigma == dddSigmadlam3 + 2 ddSigmadlam2 + dSigmadlam,
        Cinstruction == "} /* end !FakeMatterOutside */",

      Cinstruction == "} else {",   (* away from star surface *)
        FSigma  == Sigma,           (* set Sigma to zero *)
      Cinstruction == "}",

    Cif == end, (* END: corot/general case *)

    (* scale with CoordFac *)
    FPsi    ==  CoordFac * FPsi,
    FB[a]   ==  CoordFac * FB[a],
    FalphaP ==  CoordFac * FalphaP,
    FSigma  ==  CoordFac * FSigma,

  (************************************************)
  (* linear case *)
  Cif == else,

    alphaP2 == alphaP*alphaP,
    alphaP3 == alphaP2*alphaP,
    Psi6    == Psi4*Psi2,
    Psi7    == Psi4*Psi3,

    gdlB == delta[a,b] dlB[a,b],
    LlB[a,b] == dlB[a,b] + dlB[b,a] -(2/3) delta[a,b] gdlB,
    LlBdo[a,b] == delta[a,c] delta[b,d] LlB[c,d], 
    LlBLlB == LlB[a,b] LlBdo[a,b],
    vecLaplB[a] == delta[b,c] (ddlB[a,b,c] + (1/3) ddlB[b,c,a]),

    (* linearized alpha == alphaP/Psi  *)
    lalpha == lalphaP/Psi - alphaP lPsi/Psi2,

    (***************************)
    (* corotation/general case *)
    (***************************)
    (**************)
    (* corotation *)
    (**************)
    Cif == ( corot ),
      (* Since vR[a] == 0, *)
      lvR[a] == 0, 
      (* linearized oouzerosqr, recall:
         oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c])*)
      loouzerosqr == 2alpha lalpha - 
                     4 Psi3 lPsi delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]) -
                     Psi4 delta[b,c] 2 (lB[b] + lvR[b]) (beta[c] + vR[c]),
      luzerosqr == -uzerosqr^2 loouzerosqr,
    (****************)
    (* general case *)
    (****************)
    Cif == else,
      lLnalpha == lalpha/alpha,
      (* dLnalpha[a] == dalpha[a]/alpha, *)
      dlalpha[a] == dlalphaP[a]/Psi - dalphaP[a] lPsi/Psi2 - 
                    lalphaP dPsi[a]/Psi2 - alphaP dlPsi[a]/Psi2 +
                    2 alphaP dPsi[a] lPsi/Psi3,
      ldLnalpha[a] == dlalpha[a]/alpha - lalpha dalpha[a]/alpha2,
      (* h == q + 1, *)
      lh == 0,
      lLnh == 0,
      dlh[a] == 0,
      ldLnh[a] == 0,

      (* wB remains const under linearization *)
      lwB[a] == 0,
      dlwB[a,b] == 0,
      (* L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]), *)
      lL2 == 2*(Psim8 wBDown[c] lwB[c] +
                   Psim6 (lwB[c] dSigma[c] + wB[c] dlSigma[c]) +
                   Psim4 dSigmaUp[c] dlSigma[c])  - 
                (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
                 4 Psim5 dSigma[c] dSigmaUp[c]) lPsi + 2 h2 lLnh,
      luzerosqr == (lL2 - 2 L2 (lalpha/alpha + lLnh))/(alpha2 h2),
      luzero == luzerosqr/(2 uzero),
      lwBDown[a] == lwB[a],
      (* dSigmaUp[a] == dSigma[a], *)
      dlSigmaUp[a] == dlSigma[a],
      Psim10 == Psim9*Psim1,
      ldL2[a] == 2*(Psim8 (lwBDown[c] dwB[c,a] + wBDown[c] dlwB[c,a]) +
                    Psim6 (dlwB[c,a] dSigma[c] + dwB[c,a] dlSigma[c] +
                           lwB[c] ddSigma[a,c] + wB[c] ddlSigma[c,a]) +
                    Psim4 (dlSigmaUp[c] ddSigma[a,c] + 
                           dSigmaUp[c] ddlSigma[a,c] )) -
                 (16 Psim9 wBDown[c] dwB[c,a] +
                  12 Psim7 (dwB[c,a] dSigma[c] + wB[c] ddSigma[a,c]) +
                  8 Psim5 dSigmaUp[c] ddSigma[a,c] ) lPsi -
                 (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
                  4 Psim5 dSigma[c] dSigmaUp[c]) dlPsi[a] +
                 (72 Psim10 wBDown[c] wB[c] + 84 Psim8 wB[c] dSigma[c] +
                  20 Psim6 dSigma[c] dSigmaUp[c]) lPsi dPsi[a] -
                 (8 Psim9 2 wBDown[c] lwB[c] +
                  12 Psim7 (lwB[c] dSigma[c] + wB[c] dlSigma[c]) +
                  8 Psim5 dlSigma[c] dSigmaUp[c]) dPsi[a] +
                 2(h2 lLnh dLnh[a] + h dlh[a]),
      lduzerosqr[a] == -2(lLnalpha + lLnh) duzerosqr[a] + 
                       (ldL2[a] - 
                        2 lL2 (dLnalpha[a] + dLnh[a]))/(alpha2 h2) - 
                       2 L2 (ldLnalpha[a] + ldLnh[a])/(alpha2 h2),
      (* vR[a] == (wB[a] + dSigma[a])/(uzero*h) - beta[a], *)
      lvR[a] == (lwB[a] + dlSigma[a])/(uzero*h) - lB[a] +
                (wB[a] + dSigma[a]) (-luzero/(uzerosqr*h) - lh/(uzero*h2)),

      divlwB == delta [b,c] dlwB[b,c],
      divlbeta == delta [b,c] dlB[b,c],
      (* ldLnalpha[a] == dlalpha[a]/alpha - lalpha dalpha[a]/alpha2, *)
      (* rho0 == Power[q/kappa, n], *)
      (* drho0[a] == (n/n+1) rho0^2/P dq[a], *)
      (* ldrho0[a] == ??? *) 
      lrho0 == 0,
      ldrho0[a] == 0, (* since lq == 0 *)
      ldLnalphaP[a] == dlalphaP[a]/alphaP - lalphaP dalphaP[a]/alphaP2,
      ldLnPsi[a] == dlPsi[a] Psim1 - lPsi dPsi[a] Psim2,
    Cif == end,
    (****************************)
    (* END: corot./general case *)
    (****************************)

    (* rho  == alpha2 (rhoE + P) uzerosqr - P,
       j[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),
       S    == 3P - rhoE + rho, *)
    (* linearized fluid vars in 3+1 *)
    lrho  == 2 alpha lalpha (rhoE + P) uzerosqr + 
             alpha2 (rhoE + P) luzerosqr,
    lj[a] == lalpha (rhoE + P) uzerosqr (vR[a]+beta[a]) + 
             alpha (rhoE + P) luzerosqr (vR[a]+beta[a]) +
             alpha (rhoE + P) uzerosqr (lvR[a]+lB[a]),
    lS    == lrho,

    ldLnalphaPsim6[a] == dlalphaP[a]/alphaP - dalphaP[a] lalphaP/alphaP2 -
                         7 dlPsi[a]/Psi + 7 dPsi[a] lPsi/Psi2,

    (* linearized equations for Psi, B[a], alphaP, Sigma *)
    FlPsi    == delta[b,c] ddlPsi[b,c] + 7 Psi6 lPsi LBLB/(32 alphaP2) -
                (Psi7 LBLB/(16 alphaP3)) lalphaP + 
                (Psi5/(32 alpha2)) 2 LBdo[a,b] LlB[a,b] +
                2Pi ( ( 5 Psi4  lPsi rho + Psi5 lrho ) (1 - CTSmod) +
                      (-3 Psim4 lPsi rhobar) (CTSmod) ),
    FlB[a]   == vecLaplB[a] - LlB[a,b] dLnalphaPsim6[b] - 
                LB[a,b] ldLnalphaPsim6[b] - 16Pi lalphaP Psi3 j[a] -
                16Pi alphaP 3 Psi2 lPsi j[a] - 16Pi alpha Psi4 lj[a],
    FlalphaP == delta[b,c] ddlalphaP[b,c] - lalphaP (
                (-7/32) Psi6 LBLB/(alphaP2) + 2Pi Psi4 (rho+2S) ) - alphaP (
                (21/16)Psi5 lPsi LBLB/(alphaP2) +
                (7/16) (Psi4/alpha2)LBdo[a,b] LlB[a,b] + 2Pi (
                 4 Psi3 lPsi (rho+2S) + Psi4 (lrho+2lS) ) ), 

    (*****************************)
    (* BEGIN: corot/general case *)
    (**************)
    (* corotation *)
    (**************)
    Cif == ( corot ),
      (* set Sigma to zero for corotation *)
      FlSigma  == lSigma,
      (* Pedro's thing (which I don't use anymore): 
         Note: inside this comment q = P/rho0 and NOT q = h-1 
      /* ell. eqn. inside stars */
      Cif == (MATTRinside),
        dvR[a,b] == 0,
        dLnrho0[a] == (n/kappa) Power[q/kappa, n-1] dq[a],
        dalpha[a] == dalphaP[a]/Psi - alphaP dPsi[a]/Psi2,
        dbeta[a,b] == dB[a,b] + epsmatrix3d[b,a,3] Omega + delta[a,b] rdotor,
        doouzerosqr[a] == 2 alpha dalpha[a] -
                        4 Psi3 dPsi[a] delta[b,c] *
                        (beta[b] + vR[b]) (beta[c] + vR[c]) -
                        2 Psi4 delta[b,c] *
                        (beta[b] + vR[b]) (dbeta[c,a] + dvR[c,a]),
        duzerosqr[a] == -uzerosqr^2 doouzerosqr[a],
        dLnuzerosqr[a] == duzerosqr[a]/uzerosqr,
        dLnuzero[a] == dLnuzerosqr[a]/2,
        dLnalphaP[a] == dalphaP[a]/alphaP,
        dLnPsi[a] == dPsi[a]/Psi,

        dlalpha[a] == dlalphaP[a]/Psi - dalphaP[a] lPsi/Psi2 - 
                      lalphaP dPsi[a]/Psi2 - alphaP dlPsi[a]/Psi2 +
                      2 alphaP dPsi[a] lPsi/Psi3,
        ldoouzerosqr[a] == 2 lalpha dalpha[a] + 2 alpha dlalpha[a] -
                           (12 Psi2 lPsi dPsi[a] + 4 Psi3 dlPsi[a]) delta[b,c] *
                           (beta[b] + vR[b]) (beta[c] + vR[c]) -
                           4 Psi3 dPsi[a] delta[b,c] *
                           2 (beta[b] + vR[b]) (lB[c] + dlSigma[c]) -
                           8 Psi3 lPsi delta[b,c] *
                           (beta[b] + vR[b]) (dbeta[c,a] + dvR[c,a]) -
                           2 Psi4 delta[b,c] * (
                            (lB[b] + dlSigma[b]) (dbeta[c,a] + dvR[c,a]) + 
                            (beta[b] + vR[b]) (dlB[c,a] + ddlSigma[c,a]) ),
        lduzerosqr[a] == -2 uzerosqr luzerosqr doouzerosqr[a] - 
                          uzerosqr^2 ldoouzerosqr[a],

        FlSigma  == delta[b,c] ddlSigma[b,c] + 
                 dlSigma[a] * ( 
                  dLnrho0[a] + dLnuzerosqr[a]/2 + dLnalphaP[a] + 5 dLnPsi[a] )+
                 (wB[a] + dSigma[a]) * (
                  lduzerosqr[a]/(2 uzerosqr) - 
                  duzerosqr[a] luzerosqr /(2 uzerosqr*uzerosqr) +
                  dlalphaP[a]/alphaP - dalphaP[a] lalphaP/alphaP2 +
                  5 dlPsi[a]/Psi - 5 dPsi[a] lPsi/Psi2 ),
      Cif == else,
        FlSigma  == lSigma,
      Cif == end, 
      *)
    (*****************************)
    (* general case (not corot.) *)
    (*****************************)
    Cif == else,
      Cinstruction == "if(MATTRinside) {", (* inside stars *)
        lLnuzero == luzero/uzero,
        dLnuzero[a] == duzerosqr[a]/(2 uzerosqr),
        lduzero[a] == lduzerosqr[a]/(2 uzero) -
                      luzerosqr duzerosqr[a] / (4 uzero uzerosqr),
        ldLnuzero[a] == lduzero[a]/(uzero) - (lLnuzero) dLnuzero[a],
        lhuzeroPsi4 == lh uzero Psi4 + h luzero Psi4 + 4 h uzero Psi3 lPsi, 
        ldLnalphaPsi2oh[a] == ldLnalphaP[a] + ldLnPsi[a] - ldLnh[a],
        ldLnalphaoh[a]     == ldLnalphaP[a] - ldLnPsi[a] - ldLnh[a],
        ldLnalphaPsi6uz[a] == ldLnalphaP[a] + ldLnuzero[a] + 5 ldLnPsi[a],
        ldrho0PLUSrho0dLnalphaPsi2oh[a] == ldrho0[a] +
                                           lrho0 dLnalphaPsi2oh[a]+
                                           rho0 ldLnalphaPsi2oh[a],
        ldrho0PLUSrho0dLnalphaoh[a]     == ldrho0[a] +
                                           lrho0 dLnalphaoh[a] +
                                           rho0 ldLnalphaoh[a],
        ldrho0PLUSrho0dLnalphaPsi6uz[a] == ldrho0[a] +
                                           lrho0 dLnalphaPsi6uz[a] +
                                           rho0 ldLnalphaPsi6uz[a],
        ddSigCoef == rho0 + rhMeps rho0max ((rho0max - rho0)/rho0max)^rhMpow,
        lddSigCoef== lrho0 *
                     (1 - rhMeps rhMpow ((rho0max - rho0)/rho0max)^rhMpowm1),
        FlSigma == ddSigCoef delta[b,c] ddlSigma[b,c] + 
                  dlSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (lwB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divlwB) -
                  h uzero Psi4 (rho0 divlbeta +
                                lB[c] drho0PLUSrho0dLnalphaPsi6uz[c]) +
                  lddSigCoef delta[b,c] ddSigma[b,c] +
                  dSigmaUp[c] ldrho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (wB[c] ldrho0PLUSrho0dLnalphaoh[c] + lrho0 divwB) -
                  2 Psim3 lPsi (wB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divwB) -
                  lhuzeroPsi4 (rho0 divbeta +
                               beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]) -
                  h uzero Psi4 (lrho0 divbeta + 
                                beta[c] ldrho0PLUSrho0dLnalphaPsi6uz[c]),
(*
        FSigma == rho0 delta[b,c] ddSigma[b,c] + 
                  dSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (wB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divwB) -
                  h uzero Psi4 (rho0 divbeta +
                                beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]),

FlSigma == rho0 delta[b,c] ddlSigma[b,c] + 
           dlSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] -
           h luzero Psi4 (rho0 divbeta +
                          beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]) -
           h uzero Psi4 beta[c] (ldLnuzero[c]),
*)
      Cinstruction == "} else if(MATTRtouch) {", (* touches star surface *)
        Cinstruction == "if(FakeMatterOutside) {", (* use fake matter *)
          Cinstruction == "if(LapSig) {",  (* Laplace Sigma *)
            FlSigma == delta[b,c] ddlSigma[b,c],
          Cinstruction == "} else {",
            (* use h=1, rho0 = -lam as fake matter *)
            hf == 1,       (* fake h *)
            Cif == FakeT0,
              rhof == -lam,  (* fake rho0=-lam *)
              drhof[a] == -dlam[a],
            Cif == else,
              rhof == -1,    (* fake rho0=-1 *)
              drhof[a] == 0,
            Cif == end,
            (* same terms as above but with dLnh = 0 *)
            dLnalphaPsi2[a] == dLnalphaP[a] + dLnPsi[a],
            dLnalpha[a]     == dLnalphaP[a] - dLnPsi[a],
            dLnalphaPsi6uz[a] == dLnalphaP[a] + 5 dLnPsi[a] + duzero[a]/uzero,
            drhofPLUSrhofdLnalphaPsi2[a] == drhof[a] + rhof dLnalphaPsi2[a],
            drhofPLUSrhofdLnalpha[a]     == drhof[a] + rhof dLnalpha[a],
            drhofPLUSrhofdLnalphaPsi6uz[a] == drhof[a] + rhof dLnalphaPsi6uz[a],
            (* same terms as above but with lrho0=0, ldrho0=0, lh=0, ldLnh=0 *)
            lLnuzero == luzero/uzero,
            dLnuzero[a] == duzerosqr[a]/(2 uzerosqr),
            lduzero[a] == lduzerosqr[a]/(2 uzero) -
                          luzerosqr duzerosqr[a] / (4 uzero uzerosqr),
            ldLnuzero[a] == lduzero[a]/(uzero) - (lLnuzero) dLnuzero[a],
            lhfuzeroPsi4 == hf luzero Psi4 + 4 hf uzero Psi3 lPsi, 
            ldLnalphaPsi2[a] == ldLnalphaP[a] + ldLnPsi[a],
            ldLnalpha[a]     == ldLnalphaP[a] - ldLnPsi[a],
            ldLnalphaPsi6uz[a] == ldLnalphaP[a] + ldLnuzero[a] + 5 ldLnPsi[a],
            ldrhofPLUSrhofdLnalphaPsi2[a] == rhof ldLnalphaPsi2[a],
            ldrhofPLUSrhofdLnalpha[a]   == rhof ldLnalpha[a],
            ldrhofPLUSrhofdLnalphaPsi6uz[a] == rhof ldLnalphaPsi6uz[a],
            (* BC *)
            FlSigma == rhof delta[b,c] ddlSigma[b,c] + 
                       dlSigmaUp[c] drhofPLUSrhofdLnalphaPsi2[c] +
                       Psim2 (lwB[c] drhofPLUSrhofdLnalpha[c] + rhof divlwB) -
                       hf uzero Psi4 (rhof divlbeta +
                                      lB[c] drhofPLUSrhofdLnalphaPsi6uz[c]) +
                       dSigmaUp[c] ldrhofPLUSrhofdLnalphaPsi2[c] +
                       Psim2 (wB[c] ldrhofPLUSrhofdLnalpha[c]) -
                       2 Psim3 lPsi (wB[c] drhofPLUSrhofdLnalpha[c] + rhof divwB) -
                       lhfuzeroPsi4 (rhof divbeta +
                                     beta[c] drhofPLUSrhofdLnalphaPsi6uz[c]) -
                       hf uzero Psi4 (beta[c] ldrhofPLUSrhofdLnalphaPsi6uz[c]),

          Cinstruction == "}",
        Cinstruction == "} else {", (* continous Sigma *)
          FlSigma == dddlSigmadlam3 + 2 ddlSigmadlam2 + dlSigmadlam,
        Cinstruction == "} /* end !FakeMatterOutside */",

      Cinstruction == "} else {", (* away from star surface *)
        FlSigma  == lSigma,       (* set lSigma to zero *)
      Cinstruction == "}",

    Cif == end, (* END: corot/general case *)

    (* scale with CoordFac *)
    FlPsi    ==  CoordFac * FlPsi,
    FlB[a]   ==  CoordFac * FlB[a],
    FlalphaP ==  CoordFac * FlalphaP,
    FlSigma  ==  CoordFac * FlSigma,

  Cif == end, (* end of nonlin/linear case *)

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] :=  K[b,a] /; !OrderedQ[{a,b}]

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

CFunctionFile = "DNS_CTS.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"DNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Abs(x)     (fabs((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n"];
  pr["extern tEoS_T0 EoS_T0[1];"];
  pr["\n\n\n"];

  pr["void DNS_CTS(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];

  pr["int iStart = Set_pdb_iStart_AtPar(\"DNSdata_rotationstate1\");\n"];
  pr["int VwApprox1 = Getv(\"DNSdata_rotationstate1\",\"VwApproximation\");\n"];
  pr["int VwApprox2 = Getv(\"DNSdata_rotationstate2\",\"VwApproximation\");\n"];
  pr["int corot1 = VwApprox1 || Getv(\"DNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = VwApprox2 || Getv(\"DNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int VwApprox, corot;\n"];
  pr["int dqFromqg = Getv(\"DNSdata_q_derivs\",\"dqg\");\n"];
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
  pr["double CTSmod = Getv(\"DNSdata_CTSmod\",\"yes\");\n"];
  pr["double rhMeps = Getd(\"DNSdata_SigmaMod_eps\");\n"];
  pr["double rhMpow = Getd(\"DNSdata_SigmaMod_pow\");\n"];
  pr["double rhMpowm1 = rhMpow - 1.;\n"];
  pr["double qmax1 = Getd(\"DNSdata_qmax1\");\n"];
  pr["double qmax2 = Getd(\"DNSdata_qmax2\");\n"];
  pr["int FakeMatterOutside = Getv(\"DNSdata_Sigma_surface_BCs\",\"FakeMatterOutside\");\n"];
  pr["int FakeT0 = Getv(\"DNSdata_FakeMatterType\",\"rhoEQ-lam\");\n"];
  pr["int LapSig = Getv(\"DNSdata_FakeMatterType\",\"LaplaceSigmaOutside\");\n"];
  pr["\n"];

  pr["tGrid *grid = vlu->grid;\n"];
  pr["int bi;\n"];
  pr["tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vldu->vlPars;\n"];
  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["int MATTRinside = (box->MATTR== INSIDE);\n\n"];
  pr["int MATTRtouch  = (box->MATTR== TOUCH);\n\n"];
  pr["int isSTAR1     = (box->SIDE == STAR1);\n\n"];
  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{FPsi, FB[a], FalphaP ,FSigma}, "vlFu"];
  prdecvl[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];

  prdecvl[{FlPsi, FlB[a], FlalphaP ,FlSigma}, "vlJdu"];
  prdecvl[{ lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvl[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvlindices[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];
  prdecvlindices[{lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvlindices[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvarname[{dPsi[a]},       "DNSdata_Psix"];
  prdecvarname[{ddPsi[a,b]},    "DNSdata_Psixx"];
  prdecvarname[{dB[a,b]}, 	"DNSdata_Bxx"];
  prdecvarname[{ddB[a,b,c]},    "DNSdata_Bxxx"];
  prdecvarname[{dalphaP[a]},    "DNSdata_alphaPx"];
  prdecvarname[{ddalphaP[a,b]}, "DNSdata_alphaPxx"];
  prdecvarname[{dSigma[a]},     "DNSdata_Sigmax"];
  prdecvarname[{ddSigma[a,b]},  "DNSdata_Sigmaxx"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];
  prdecvarname[{z},      "z"];

  (* prdecvarname[{g[a,b]}, "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{alpha},   "alpha"];
  prdecvarname[{beta[a]}, "betax"];
  prdecvarname[{q},       "DNSdata_q"];
  prdecvarname[{wB[a]},   "DNSdata_wBx"];
  prdecvarname[{dq[a]},   "DNSdata_qx"];
  prdecvarname[{dwB[a,b]},"DNSdata_wBxx"];
  prdecvarname[{VR[a]},   "DNSdata_VRx"];

  prdecvarname[{lam},            "X"];
  prdecvarname[{dlam[a]},        "dXdx"];
  prdecvarname[{dSigmadlam},     "DNSdata_SigmaX"];
  prdecvarname[{ddSigmadlam2},   "DNSdata_SigmaXX"];
  prdecvarname[{dddSigmadlam3},  "DNSdata_SigmaXXX"];
  prdecvarname[{dlSigmadlam},    "DNSdata_lSigmaX"];
  prdecvarname[{ddlSigmadlam2},  "DNSdata_lSigmaXX"];
  prdecvarname[{dddlSigmadlam3}, "DNSdata_lSigmaXXX"];
  prdecvarname[{rhobar},         "DNSdata_rhobar"];

  prdecvarname[{CoordFac},     "DNSdata_CoordFac"];

  pr["\n"];
];    

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  (* some thread global C vars we already need here *)
  pr["double xmax, xC;\n"];
  pr["double omegMOmeg1, omegMOmeg2, omegMOmeg3;\n"];
  pr["double rho0max, Pmax, rhoEmax, drho0dhm1max;\n"];

  (* do nothing and continue if block is not in box bi *)
  pr["if(blkinfo!=NULL) if(bi!=blkinfo->bi) continue;\n"];

  pr["if(nonlin) {\n"];
    pr["D_and_DD_of_S(box, index_Psi, \ 
		      Ind(\"DNSdata_Psix\"), Ind(\"DNSdata_Psixx\"));\n"];
    pr["D_and_DD_of_Sa(box, index_B1, \
		       Ind(\"DNSdata_Bxx\"), Ind(\"DNSdata_Bxxx\"));\n"];
    pr["D_and_DD_of_S(box, index_alphaP, \
		      Ind(\"DNSdata_alphaPx\"), Ind(\"DNSdata_alphaPxx\"));\n"];
    pr["D_and_DD_of_S(box, index_Sigma, \
		      Ind(\"DNSdata_Sigmax\"), Ind(\"DNSdata_Sigmaxx\"));\n"];
    pr["spec_Deriv2(box, 1, Sigma, ddSigmadlam2);\n"];
    pr["spec_Deriv1(box, 1, ddSigmadlam2, dddSigmadlam3);\n"];
    pr["spec_Deriv1(box, 1, Sigma, dSigmadlam);\n"];
  pr["} else {\n"];
    pr["if(blkinfo!=NULL) {\n"];
    pr["  if(blkinfo->vari == index_lPsi) \
                         D_and_DD_of_S(box, index_lPsi, \ 
				       index_dlPsi1, index_ddlPsi11);\n"];
    pr["  if(blkinfo->vari == index_lB1) \
                         D_and_DD_of_S(box, index_lB1, \
				       index_dlB11, index_ddlB111);\n"];
    pr["  if(blkinfo->vari == index_lB2) \
                         D_and_DD_of_S(box, index_lB2, \
				       index_dlB21, index_ddlB211);\n"];
    pr["  if(blkinfo->vari == index_lB3) \
                         D_and_DD_of_S(box, index_lB3, \
				       index_dlB31, index_ddlB311);\n"];
    pr["  if(blkinfo->vari == index_lalphaP) \
                         D_and_DD_of_S(box, index_lalphaP, \
				       index_dlalphaP1, index_ddlalphaP11);\n"];
    pr["  if(blkinfo->vari == index_lSigma) {\n"];
    pr["    D_and_DD_of_S(box, index_lSigma, \
			  index_dlSigma1, index_ddlSigma11);\n"];
    pr["    spec_Deriv2(box, 1, lSigma, ddlSigmadlam2);\n"];
    pr["    spec_Deriv1(box, 1, ddlSigmadlam2, dddlSigmadlam3);\n"];
    pr["    spec_Deriv1(box, 1, lSigma, dlSigmadlam);\n"];
    pr["  } /* end blkinfo->vari == index_lSigma */\n"];
    pr["} /* end blkinfo!=NULL  */\n"];
    pr["else {\n"];
    pr["  D_and_DD_of_S(box, index_lPsi, \ 
			index_dlPsi1, index_ddlPsi11);\n"];
    pr["  D_and_DD_of_Sa(box, index_lB1, \
			 index_dlB11, index_ddlB111);\n"];
    pr["  D_and_DD_of_S(box, index_lalphaP, \
			index_dlalphaP1, index_ddlalphaP11);\n"];
    pr["  D_and_DD_of_S(box, index_lSigma, \
			index_dlSigma1, index_ddlSigma11);\n"];
    pr["  spec_Deriv2(box, 1, lSigma, ddlSigmadlam2);\n"];
    pr["  spec_Deriv1(box, 1, ddlSigmadlam2, dddlSigmadlam3);\n"];
    pr["  spec_Deriv1(box, 1, lSigma, dlSigmadlam);\n"];
    pr["} /* end else */\n"];
  pr["} /* end linear case */\n"];

  pr["FirstDerivsOf_Sa(box, Ind(\"DNSdata_wBx\"), \
					 Ind(\"DNSdata_wBxx\"));\n"];
  pr["if(dqFromqg) {\n"];
    pr["FirstDerivsOf_S(box,  Ind(\"DNSdata_qg\"), \
				     Ind(\"DNSdata_qx\"));\n"];
  pr["} else {\n"];
    pr["FirstDerivsOf_S(box,  Ind(\"DNSdata_q\"), \
				     Ind(\"DNSdata_qx\"));\n"];
  pr["} /* end if */\n"]; 

  (* which star are we considering? *)
  pr["VwApprox = corot = 0;\n"];
  pr["if(isSTAR1) {\n"];  
    pr["if(corot1)    corot = 1;\n"];
    pr["if(VwApprox1) VwApprox = 1;\n"];
    pr["xmax = xmax1;\n"];
    pr["omegMOmeg1 = omegax1;\n"];
    pr["omegMOmeg2 = omegay1;\n"];
    pr["omegMOmeg3 = omegaz1 - Omega;\n"];
    pr["EoS_T0->vars_from_hm1(qmax1, &rho0max, 
                                 &Pmax, &rhoEmax, &drho0dhm1max);\n"];
  pr["} else {\n"];
    pr["if(corot2)    corot = 1;\n"];
    pr["if(VwApprox2) VwApprox = 1;\n"];
    pr["xmax = xmax2;\n"];
    pr["omegMOmeg1 = omegax2;\n"];
    pr["omegMOmeg2 = omegay2;\n"];
    pr["omegMOmeg3 = omegaz2 - Omega;\n"];
    pr["EoS_T0->vars_from_hm1(qmax2, &rho0max, 
                                 &Pmax, &rhoEmax, &drho0dhm1max);\n"];
  pr["} /* end if */\n"];

  (* center of circle *)
  pr["xC = xCM + ecc * (xmax - xCM);\n"];

  (* loop of all points *)
  pr["SGRID_LEVEL3_Pragma(omp parallel for)\n"];
  pr["forallpoints(box, ijk) {\n"];
  pr["double rho0, P, rhoE, drho0dhm1; /* declare them */\n"];

  pr["/* Jetzt geht's los! */\n"];

];
(* auxillary variables are automatically inserted here *)

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
