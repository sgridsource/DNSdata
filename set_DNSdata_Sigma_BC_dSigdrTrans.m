(* set_DNSdata_Sigma_BC_dSigdrTrans.m
   Wolfgang Tichy  10/2018       *)

(* compute residuals of DNS ham, mom, alphaP and Sigma eqns *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP, FSigma,
             dSigma[a], ddSigma[a,b],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
             q, wB[a], dq[a], dlam[a], x,y,z, X,Y,Z, dSigmadlam,dlSigmadlam,
             ddSigmadlam2,ddlSigmadlam2, dddSigmadlam3,dddlSigmadlam3,
             dSigmaindlam, dlSigmaindlam, ddSigmaindlam2, ddlSigmaindlam2}

constvariables = {OmegaCrossR[a], xrdotor[a]}

(* compute in this order *)
tocompute = {

  (* do nothing and continue if block is not in box bi *)
  Cinstruction == "if(blkinfo != NULL) {
                     if(blkinfo->bi   != bi) continue;
                     if(blkinfo->vari != index_lSigma) continue;  }",

  (* Use Sigma=0 as BC if corot *)
  Cif == ( ((isSTAR1 && corot1) || ((!isSTAR1) && corot2)) && hasSSURF ),

    Cif == nonlin, (* non-linear case *)
      (* go over lam=0 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k);",
        FSigma  == Sigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",
      (* go over lam=1 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",
        FSigma  == Sigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",

    Cif == else,   (* linear case *)
      (* go over lam=0 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k);",
        FlSigma  == lSigma,
      Cinstruction == "} /* end forplane1 */",
      (* go over lam=1 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",
        FlSigma  == lSigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",

    Cif == end,

    Cinstruction == "continue; /* for corot we are done with this box */",
  Cif == end,

  (**********************)
  (* Start general case *)
  (**********************)

  (***********************************************)
  (* in box that touches matter for general case *)
  (***********************************************)
  Cif == ( MATTRtouch ),

    Cinstruction == "int    ind0, indin, biin, n1in,n2in,n3in;",
    Cinstruction == "double Sig, Sigin, lSig, lSigin;",
    Cinstruction == "double *Sigmain;",
    Cinstruction == "double *lSigmain;",
    Cinstruction == "double *dSigmaindlam;",
    Cinstruction == "double *dlSigmaindlam;",
    Cinstruction == "double *ddSigmaindlam2;",
    Cinstruction == "double *ddlSigmaindlam2;",
    Cinstruction == "biin = bi - 6; /* only if CubSp are arranged my way */",
    Cinstruction == "n1in = grid->box[biin]->n1;
                     n2in = grid->box[biin]->n2;
                     n3in = grid->box[biin]->n3;\n
             Sigmain      = grid->box[biin]->v[index_Sigma];
            lSigmain      = grid->box[biin]->v[index_lSigma];
            dSigmaindlam  = grid->box[biin]->v[index_DNSdata_SigmaX];
           dlSigmaindlam  = grid->box[biin]->v[index_DNSdata_lSigmaX];
           ddSigmaindlam2 = grid->box[biin]->v[index_DNSdata_SigmaXX];
          ddlSigmaindlam2 = grid->box[biin]->v[index_DNSdata_lSigmaXX];",
    Cinstruction == "if(n2in!=n2 || n3in!=n3) errorexit(\"we need n2in=n2 and n3in=n3\");",

    Cif == nonlin, (* non-linear case *)
 
     (* take derivs needed *)
      Cinstruction == "spec_Deriv2(box, 1, Sigma, ddSigmadlam2);",
      Cinstruction == "spec_Deriv1(box, 1, ddSigmadlam2, dddSigmadlam3);",
      Cinstruction == "spec_Deriv1(box, 1, Sigma, dSigmadlam);",

      (* set Sigma's equal at star surfaces, impose it at i=1 *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 1){ ijk=Index(i,j,k);",
      Cinstruction == "ind0   = Ind_n1n2(0,j,k,n1,n2);
                       indin  = Ind_n1n2(n1in-1,j,k,n1in,n2in);
                       Sig    = Sigma[ind0];
                       Sigin  = Sigmain[indin];",
        FSigma == Sig - Sigin,
      Cinstruction == "} /* endfor */",

      (* set dSigma/dr equal at star surfaces, 
         impose it at i=0. The n^i is the normal vec. *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      dSig == dSigmadlam,
      AA == Y, (* get A,B coords *)
      BB == Z,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(box, ijk,
                                                 0.,AA,BB, &rr, &drdlam);",
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      dSigin == dSigmaindlam,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(grid->box[biin], ijk,
                                                 1.,AA,BB, &rr, &drdlamin);",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      LoLin == drdlam/drdlamin,
      fTrans == LoLin OuterSigmaTransitionD1,
      FSigma == dSig - dSigin fTrans,
      Cinstruction == "} /* endfor */",

      (* set d^2Sigma/dr^2 equal at star surfaces, 
         impose it at i=2. The n^i is the normal vec. *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 2){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      AA == Y, (* get A,B coords *)
      BB == Z,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(box, ijk,
                                                 0.,AA,BB, &rr, &drdlam);",
      ddSig == ddSigmadlam2,
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      ddSigin == ddSigmaindlam2,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(grid->box[biin], ijk,
                                                 1.,AA,BB, &rr, &drdlamin);",
      Cinstruction == "ijk=Index(2,j,k); /* set index to i=2 */",
      LoLin == drdlam/drdlamin,
      fTrans == LoLin LoLin OuterSigmaTransitionD2,
      FSigma == ddSig - ddSigin fTrans,
      Cinstruction == "} /* endfor */",

    Cif == else,   (* linear case *)
      (* take derivs needed *)
      Cinstruction == "spec_Deriv2(box, 1, lSigma, ddlSigmadlam2);",
      Cinstruction == "spec_Deriv1(box, 1, ddlSigmadlam2, dddlSigmadlam3);",
      Cinstruction == "spec_Deriv1(box, 1, lSigma, dlSigmadlam);",

      (* set Sigma's equal at star surfaces, impose it at i=1 *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 1){ ijk=Index(i,j,k);",
      Cinstruction == "ind0   = Ind_n1n2(0,j,k,n1,n2);
                       indin  = Ind_n1n2(n1in-1,j,k,n1in,n2in);
                       lSig   = lSigma[ind0];
                       lSigin = lSigmain[indin];",
        FlSigma == lSig - lSigin,
      Cinstruction == "} /* endfor */",

      (* set dSigma/dr equal at star surfaces, 
         impose it at i=0. The n^i is the normal vec. *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      dlSig == dlSigmadlam,
      AA == Y, (* get A,B coords *)
      BB == Z,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(box, ijk,
                                                 0.,AA,BB, &rr, &drdlam);",
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      dlSigin == dlSigmaindlam,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(grid->box[biin], ijk,
                                                 1.,AA,BB, &rr, &drdlamin);",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      LoLin == drdlam/drdlamin,
      fTrans == LoLin OuterSigmaTransitionD1,
      FlSigma == dlSig - dlSigin fTrans,
      Cinstruction == "} /* endfor */",

      (* set d^2Sigma/dr^2 equal at star surfaces, 
         impose it at i=2. The n^i is the normal vec. *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 2){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      ddlSig == ddlSigmadlam2,
      AA == Y, (* get A,B coords *)
      BB == Z,
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(box, ijk,
                                                 0.,AA,BB, &rr, &drdlam);",
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      Cinstruction == "r_dr_dlam_of_lamAB_CubSph(grid->box[biin], ijk,
                                                 1.,AA,BB, &rr, &drdlamin);",
      ddlSigin == ddlSigmaindlam2,
      Cinstruction == "ijk=Index(2,j,k); /* set index to i=2 */",
      LoLin == drdlam/drdlamin,
      fTrans == LoLin LoLin OuterSigmaTransitionD2,
      FlSigma == ddlSig - ddlSigin fTrans,
      Cinstruction == "} /* endfor */",

    Cif == end, (* end linear case *)

  Cif == end, (* end case of boxes that touch matter *)

  (***********************************************************************)
  (* now treat the case where matter is inside the box with star surface *)
  (***********************************************************************)
  Cif == ( MATTRinside && hasSSURF ),
    (* if we get here bi is that of a cubedSph and there is no corot
       in this box *)
    Cif == dqFromqg,
      Cinstruction == "FirstDerivsOf_S(box, Ind(\"DNSdata_qg\"), \
					Ind(\"DNSdata_qx\"));",
    Cif == else,
      Cinstruction == "FirstDerivsOf_S(box, Ind(\"DNSdata_q\"), \
					Ind(\"DNSdata_qx\"));",
    Cif == end, 

    (* non-linear case: *)
    Cif == nonlin,
      Cinstruction == "FirstDerivsOf_S(box, index_Sigma, \
                                       Ind(\"DNSdata_Sigmax\"));",

      (* which star are we considering? *)
      Cif == (isSTAR1),
        xmax == xmax1,
      Cif == else,
        xmax == xmax2,
      Cif == end,

      (* go over lam=1 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",

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
        Psi2  == Psi*Psi,
        Psi4  == Psi2*Psi2,
        Psim2 == 1/Psi2,
        Psim4 == Psim2*Psim2,
        Psim6 == Psim4*Psim2,
        alpha == alphaP/Psi,
        alpha2 == alpha alpha,
        (* terms needed *)
        h == q + 1,
        h2 == h h,
        DSigmaUp[a] == Psim4 dSigma[a],
        dSigmaUp[a] == dSigma[a],
        w[a] == Psim6 wB[a],
        wBDown[a] == wB[a],
        wDown[a] == Psim2 wBDown[a],
        L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
        uzerosqr == L2/(alpha2 h2),
        uzero == sqrt[uzerosqr],

        (* do we use dlam instead of dq as surface normal? *)
        Cif == dQFromdlam,
          dQ[a] == dlam[a],
        Cif == else,
          dQ[a] == dq[a],
        Cif == end,

        (* actual Physical BC *)
        Cif == ImposeActualBC,
          FSigma == dSigmaUp[c] dQ[c] - h uzero Psi4 beta[c] dQ[c],
          (* add extra term with wB *)
          FSigma == FSigma + Psim2 wB[c] dQ[c],
        Cif == end,

      Cinstruction == "} /* end forplane1 */",


    (* linear case: *)
    Cif == else,
      Cinstruction == "FirstDerivsOf_S(box, index_lSigma, index_dlSigma1);",


      (* which star are we considering? *)
      Cif == (isSTAR1),
        xmax == xmax1,
      Cif == else,
        xmax == xmax2,
      Cif == end,

      (* go over lam=1 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",

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
        Psi2  == Psi*Psi,
        Psi3  == Psi2*Psi,
        Psi4  == Psi2*Psi2,
        Psim2 == 1/Psi2,
        Psim3 == 1/Psi3,
        Psim4 == Psim2*Psim2,
        Psim6 == Psim4*Psim2,
        Psim8 == Psim4*Psim4,
        Psim5 == Psim6*Psi,
        Psim7 == Psim5*Psim2,
        Psim9 == Psim7*Psim2,
        alpha == alphaP/Psi,
        alpha2 == alpha alpha,
        (* non-linear terms needed *)
        h == q + 1,
        h2 == h h,
        DSigmaUp[a] == Psim4 dSigma[a],
        dSigmaUp[a] == dSigma[a],
        w[a] == Psim6 wB[a],
        wBDown[a] == wB[a],
        wDown[a] == Psim2 wBDown[a],
        L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
        uzerosqr == L2/(alpha2 h2),
        uzero == sqrt[uzerosqr],

        (* linearized terms *)
        lq     == 0,
        dlQ[a] == 0,
        lh   == 0,
        lLnh == 0,
        (* wB remains const under linearization *)
        lwB[a] == 0,
        dlwB[a,b] == 0,
        lalpha == lalphaP/Psi - alphaP lPsi/Psi2,   
        (* dSigmaUp[a] == dSigma[a], *)
        dlSigmaUp[a] == dlSigma[a],
        lL2 == 2*(Psim8 wBDown[c] lwB[c] +
                     Psim6 (lwB[c] dSigma[c] + wB[c] dlSigma[c]) +
                     Psim4 dSigmaUp[c] dlSigma[c]  - 
                  (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
                   4 Psim5 dSigma[c] dSigmaUp[c]) lPsi + 2 h2 lLnh),
        luzerosqr == (lL2 - 2 L2 (lalpha/alpha + lLnh))/(alpha2 h2),
        luzero == luzerosqr/(2 uzero),
        lhuzeroPsi4beta[a] == h (luzero Psi4 beta[a] +
                                 4 uzero Psi3 lPsi beta[a] +
                                 uzero Psi4 lB[a]) + lh uzero Psi4 beta[a],

       (* do we use dlam instead of dq as surface normal? *)
        Cif == dQFromdlam,
          dQ[a] == dlam[a],
        Cif == else,
          dQ[a] == dq[a],
        Cif == end,

        (* actual Physical BC *)
        Cif == ImposeActualBC,
          FlSigma  == dlSigmaUp[c] dQ[c] - lhuzeroPsi4beta[c] dQ[c] +
                      dSigmaUp[c] dlQ[c] - h uzero Psi4 beta[c] dlQ[c],
          (* add extra term with wB *)
          FlSigma == FlSigma + Psim2 lwB[c] dQ[c] - 2 Psim3 lPsi wB[c] dQ[c] +
                               Psim2 wB[c] dlQ[c],
        Cif == end,

      Cinstruction == "} /* end forplane1 */",

    Cif == end, (* end of nonlin/linear case *)

  Cif == end,


  (* impose extra condition in one starbox *)
  Cif == ( MATTRinside && isVolAvBox ),

    Cif == nonlin, (* non-linear case *)

      (* we only impose InnerVolIntZero in one box *)
      Cif == InnerVolIntZero, (* (AddInnerVolIntToBC || InnerVolIntZero), *)
        Cinstruction == "VolAvSigma = BoxVolumeIntegral(box, index_Sigma);",
      Cif == end,
      Cif == InnerSumZero, (* (AddInnerSumToBC || InnerSumZero), *)
        Cinstruction == "VolAvSigma = 0.0;",
        Cinstruction == "forallpoints(box, ijk) {",
        Cinstruction == "VolAvSigma += Sigma[ijk];",
        Cinstruction == "} /* endfor */",
      Cif == end,

      (* find VolAvSigma0 that we would like to have in our condition *)
      Cif == (isSTAR1),
        Cinstruction == "VolAvSigma0 = VolAvSigma1;",
      Cif == else,
        Cinstruction == "VolAvSigma0 = VolAvSigma2;",
      Cif == end,
      Cinstruction == "//printf(\"VolAvSigma-VolAvSigma0=%g\\n\",VolAvSigma-VolAvSigma0);",

      (* impose conditions at this point: *)
      Cinstruction == "ijk = Index(n1/2, n2/2, n3/2);",
      Cinstruction == "//printf(\"(%d)\", ijk);",

      (* set Sigma to zero at ijk *)
      Cif == SigmaZeroAtPoint,
        FSigma == Sigma - VolAvSigma0,
      Cif == end,

      (* set VolAvSigma to zero at ijk *)
      Cif == (InnerVolIntZero || InnerSumZero),
        FSigma == VolAvSigma - VolAvSigma0,
      Cif == end,

    Cif == else,   (* linear case *)

      Cif == InnerVolIntZero, (* (AddInnerVolIntToBC || InnerVolIntZero), *)
        Cinstruction == "VolAvlSigma = BoxVolumeIntegral(box, index_lSigma);",
      Cif == end,
      Cif == InnerSumZero, (* (AddInnerSumToBC || InnerSumZero), *)
        Cinstruction == "VolAvlSigma = 0.0;",
        Cinstruction == "forallpoints(box, ijk) {",
        Cinstruction == "VolAvlSigma += lSigma[ijk];",
        Cinstruction == "} /* endfor */",
      Cif == end,
      Cinstruction == "//if(VolAvlSigma!=0.0) printf(\"box->b=%d VolAvlSigma=%g\\n\",box->b,VolAvlSigma);",

      (* impose conditions at this point: *)
      Cinstruction == "ijk = Index(n1/2, n2/2, n3/2);",
      Cinstruction == "//printf(\"|%d|\", ijk);",

      (* set Sigma to zero at ijk *)
      Cif == SigmaZeroAtPoint,
        FlSigma == lSigma,
      Cif == end,

      (* set VolAvSigma to zero at ijk *)
      Cif == (InnerVolIntZero || InnerSumZero),
        FlSigma == VolAvlSigma,
      Cif == end,

    Cif == end, (* end of nonlin/linear case *)
  Cif == end,

  (* set FSigma and FlSigma such that Sigma remains unchanged on inside *)
  Cif == ( MATTRinside && KeepInnerSigma ),
    Cif == nonlin, (* non-linear case *)
      Cinstruction == "forallpoints(box, ijk) {",
        FSigma  == 0,  (* Zero error ==> do not touch Sigma *)
      Cinstruction == "} /* endfor */",
    Cif == else,   (* linear case *)
      Cinstruction == "forallpoints(box, ijk) {",
        FlSigma  == lSigma,  (* make lin. matrix delta_ij ==> no correction *)
      Cinstruction == "} /* endfor */",
    Cif == end,
  Cif == end,

  (* for testing: set Sigma zero everywhere on outside *)
  Cif == ( (MATTRaway || MATTRtouch) && SigmaZeroInOuterBoxes),
    Cif == nonlin, (* non-linear case *)
      Cinstruction == "forallpoints(box, ijk) {",
        FSigma  == Sigma,  (* set Sigma=0 *)
      Cinstruction == "} /* endfor */",
    Cif == else,   (* linear case *)
      Cinstruction == "forallpoints(box, ijk) {",
        FlSigma  == lSigma,  (* set Sigma=0 *)
      Cinstruction == "} /* endfor */",
    Cif == end,
  Cif == end,

  Cinstruction == "/* end all */\n"
}


(* symmetries *)
ddPsi[a_,b_]     := ddPsi[b,a]    /; !OrderedQ[{a,b}]
ddB[a_,b_,c_]    := ddB[a,c,b]    /; !OrderedQ[{b,c}]
ddalphaP[a_,b_]  := ddalphaP[b,a] /; !OrderedQ[{a,b}]
ddSigma[a_,b_]   := ddSigma[b,a]  /; !OrderedQ[{a,b}]
ddSigmain[a_,b_] := ddSigmain[b,a]  /; !OrderedQ[{a,b}]
ddSig[a_,b_]     := ddSig[b,a]  /; !OrderedQ[{a,b}]
ddSigin[a_,b_]   := ddSigin[b,a]  /; !OrderedQ[{a,b}]

ddlPsi[a_,b_]     := ddlPsi[b,a]    /; !OrderedQ[{a,b}]
ddlB[a_,b_,c_]    := ddlB[a,c,b]    /; !OrderedQ[{b,c}]
ddlalphaP[a_,b_]  := ddlalphaP[b,a] /; !OrderedQ[{a,b}]
ddlSigma[a_,b_]   := ddlSigma[b,a]  /; !OrderedQ[{a,b}]
ddlSigmain[a_,b_] := ddlSigmain[b,a]  /; !OrderedQ[{a,b}]
ddlSig[a_,b_]     := ddlSig[b,a]  /; !OrderedQ[{a,b}]
ddlSigin[a_,b_]   := ddlSigin[b,a]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "set_DNSdata_Sigma_BC_dSigdrTrans.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"DNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void set_DNSdata_Sigma_BC_dSigdrTrans(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];

  pr["int VwApprox1 = Getv(\"DNSdata_rotationstate1\",\"VwApproximation\");\n"];
  pr["int VwApprox2 = Getv(\"DNSdata_rotationstate2\",\"VwApproximation\");\n"];
  pr["int corot1 = VwApprox1 || Getv(\"DNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = VwApprox2 || Getv(\"DNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int dqFromqg = Getv(\"DNSdata_q_derivs\",\"dqg\");\n"];
  pr["int dQFromdlam = Getv(\"DNSdata_drho0_inBC\",\"dlam\");\n"];
  pr["int SigmaZeroAtPoint = Getv(\"DNSdata_Sigma_surface_BCs\",\"ZeroAtPoint\");\n"];
  pr["//int AddInnerVolIntToBC = Getv(\"DNSdata_Sigma_surface_BCs\",\"AddInnerVolIntToBC\");\n"];
  pr["int InnerVolIntZero = Getv(\"DNSdata_Sigma_surface_BCs\",\"InnerVolIntZero\");\n"];
  pr["//int AddInnerSumToBC = Getv(\"DNSdata_Sigma_surface_BCs\",\"AddInnerSumToBC\");\n"];
  pr["int InnerSumZero = Getv(\"DNSdata_Sigma_surface_BCs\",\"InnerSumZero\");\n"];
  pr["int SigmaZeroInOuterBoxes = Getv(\"DNSdata_Sigma_surface_BCs\",\"ZeroInOuterBoxes\");\n"];
  pr["int noBCs = Getv(\"DNSdata_Sigma_surface_BCs\",\"none\");\n"];
  pr["int KeepInnerSigma = Getv(\"DNSdata_KeepInnerSigma\",\"yes\");\n"];
  pr["int ImposeActualBC = !Getv(\"DNSdata_Sigma_surface_BCs\",\"EllEqn\");\n"];
  pr["double Omega = Getd(\"DNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"DNSdata_x_CM\");\n"];
  pr["double ecc = Getd(\"DNSdata_ecc\");\n"];  
  pr["double rdot = Getd(\"DNSdata_rdot\");\n"];  
  pr["double xmax1 = Getd(\"DNSdata_actual_xmax1\");\n"];
  pr["double xmax2 = Getd(\"DNSdata_actual_xmax2\");\n"];
  pr["double cxmax1 = Getd(\"DNSdata_xmax1\");\n"];
  pr["double cxmax2 = Getd(\"DNSdata_xmax2\");\n"];
  pr["double VolAvSigma1 = Getd(\"DNSdata_desired_VolAvSigma1\");\n"];
  pr["double VolAvSigma2 = Getd(\"DNSdata_desired_VolAvSigma2\");\n"];
  pr["double VolAvSigma, VolAvSigma0, VolAvlSigma;\n"];
  pr["double OuterSigmaTransitionD1 = 1.0;\n"];
  pr["double OuterSigmaTransitionD2 = 1.0;\n"];
  pr["double rr, drdlam, drdlamin;\n"];
  pr["\n"];

  pr["tGrid *grid = vlu->grid;\n"];
  pr["int bi;\n"];
  pr["tVarBoxSubboxIndices *blkinfo = (tVarBoxSubboxIndices*) vldu->vlPars;\n"];
  pr["\n"];
  pr["\n"];

  pr["/* do nothing if noBCs, i.e. DNSdata_Sigma_surface_BCs = none */\n"];
  pr["if(noBCs) return;\n\n\n"];

  pr["/* parse some pars: */\n"];
  pr["/* check if DNSdata_InnerToOuterSigmaTransition is only C1 or C0 */\n"];
  pr["if(Getv(\"DNSdata_InnerToOuterSigmaTransition\",\"C1\"))\n"];
  pr["  OuterSigmaTransitionD2 = 0.0;\n"]; 
  pr["if(Getv(\"DNSdata_InnerToOuterSigmaTransition\",\"C0\"))\n"]; 
  pr["  OuterSigmaTransitionD2 = OuterSigmaTransitionD1 = 0.0;\n"]; 
  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["int n1 = box->n1;\n"];
  pr["int n2 = box->n2;\n"];
  pr["int n3 = box->n3;\n"];
  pr["int i,j,k, pln;\n\n"];
  pr["int isSTAR1     = (box->SIDE == STAR1);\n"];
  pr["int MATTRinside = (box->MATTR== INSIDE);\n"];
  pr["int MATTRtouch  = (box->MATTR== TOUCH);\n"];
  pr["int MATTRaway   = (box->MATTR== AWAY);\n"];
  pr["int hasSSURF    = (box->BOUND== SSURF);\n"];
  pr["//int isXinDom    = (box->CI->dom == box->SIDE - STAR1);\n"];
  pr["//int isVolAvBox  = (MATTRinside && hasSSURF && isXinDom);\n"];
  pr["int isCube = (box->CI->type == 0);\n"];
  pr["int isVolAvBox  = (MATTRinside && isCube);\n"];

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

  prdecvarname[{dSigma[a]},     "DNSdata_Sigmax"];
  prdecvarname[{ddSigma[a,b]},  "DNSdata_Sigmaxx"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];
  prdecvarname[{z},      "z"];

  prdecvarname[{X},      "X"];
  prdecvarname[{Y},      "Y"];
  prdecvarname[{Z},      "Z"];

  (* prdecvarname[{g[a,b]}, "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{q},       "DNSdata_q"];
  prdecvarname[{wB[a]},   "DNSdata_wBx"];
  prdecvarname[{dq[a]},   "DNSdata_qx"];
  prdecvarname[{dlam[a]}, "dXdx"];      

  prdecvarname[{dSigmadlam},     "DNSdata_SigmaX"];
  prdecvarname[{ddSigmadlam2},   "DNSdata_SigmaXX"];
  prdecvarname[{dddSigmadlam3},  "DNSdata_SigmaXXX"];
  prdecvarname[{dlSigmadlam},    "DNSdata_lSigmaX"];
  prdecvarname[{ddlSigmadlam2},  "DNSdata_lSigmaXX"];
  prdecvarname[{dddlSigmadlam3}, "DNSdata_lSigmaXXX"];

  pr["\n"];
];    

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

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
