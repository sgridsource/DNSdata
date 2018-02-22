(* DNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute residuals of DNS ham, mom, alphaP and Sigma eqns *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP, FSigma,
             dSigma[a], ddSigma[a,b],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
             q, wB[a], dq[a], dlam[a], x, y, z, dSigmadlam,dlSigmadlam,
             ddSigmadlam2,ddlSigmadlam2, dddSigmadlam3,dddlSigmadlam3,
             dSigmain[a], dlSigmain[a], ddSigmain[a,b], ddlSigmain[a,b]}

constvariables = {OmegaCrossR[a], xrdotor[a]}

(* compute in this order *)
tocompute = {

  (* do nothing if we are not in a box with a star surface *)
  Cif == (!hasSSURF),
    Cinstruction == "continue;",
  Cif == end,

  (* do nothing and continue if block is not in box bi *)
  Cinstruction == "if(blkinfo != NULL) {
                     if(blkinfo->bi   != bi) continue;
                     if(blkinfo->vari != index_lSigma) continue;  }",

  (* Use Sigma=0 as BC if corot *)
  Cif == ( (isSTAR1 && corot1) || ((!isSTAR1) && corot2) ),

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

  (* in box that touches matter for general case *)
  Cif == ( MATTRtouch ),

    Cinstruction == "int    ind0, indin, biin, n1in,n2in,n3in;",
    Cinstruction == "double Sig, Sigin, lSig, lSigin;",
    Cinstruction == "double *Sigmain;",
    Cinstruction == "double *lSigmain;",
    Cinstruction == "double *dSigmain1, *dSigmain2, *dSigmain3;",
    Cinstruction == "double *dlSigmain1, *dlSigmain2, *dlSigmain3;",
    Cinstruction == "double *ddSigmain11, *ddSigmain12, *ddSigmain13,
                            *ddSigmain22, *ddSigmain23, *ddSigmain33;",
    Cinstruction == "double *ddlSigmain11, *ddlSigmain12, *ddlSigmain13,
                            *ddlSigmain22, *ddlSigmain23, *ddlSigmain33;",
    Cinstruction == "biin = bi - 6; /* only if CubSp are arranged my way */",
    Cinstruction == "n1in = grid->box[biin]->n1;
                     n2in = grid->box[biin]->n2;
                     n3in = grid->box[biin]->n3;\n
               Sigmain    = grid->box[biin]->v[index_Sigma];
               lSigmain   = grid->box[biin]->v[index_lSigma];
               dSigmain1  = grid->box[biin]->v[index_DNSdata_Sigmax+0];
               dSigmain2  = grid->box[biin]->v[index_DNSdata_Sigmax+1];
               dSigmain3  = grid->box[biin]->v[index_DNSdata_Sigmax+2];
              dlSigmain1  = grid->box[biin]->v[index_dlSigma1];
              dlSigmain2  = grid->box[biin]->v[index_dlSigma2];
              dlSigmain3  = grid->box[biin]->v[index_dlSigma3];
              ddSigmain11 = grid->box[biin]->v[index_DNSdata_Sigmaxx+0];
              ddSigmain12 = grid->box[biin]->v[index_DNSdata_Sigmaxx+1];
              ddSigmain13 = grid->box[biin]->v[index_DNSdata_Sigmaxx+2];
              ddSigmain22 = grid->box[biin]->v[index_DNSdata_Sigmaxx+3];
              ddSigmain23 = grid->box[biin]->v[index_DNSdata_Sigmaxx+4];
              ddSigmain33 = grid->box[biin]->v[index_DNSdata_Sigmaxx+5];
             ddlSigmain11 = grid->box[biin]->v[index_ddlSigma11];
             ddlSigmain12 = grid->box[biin]->v[index_ddlSigma12];
             ddlSigmain13 = grid->box[biin]->v[index_ddlSigma13];
             ddlSigmain22 = grid->box[biin]->v[index_ddlSigma22];
             ddlSigmain23 = grid->box[biin]->v[index_ddlSigma23];
             ddlSigmain33 = grid->box[biin]->v[index_ddlSigma33];",
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

      (* set n^i dSigma/dx^i equal at star surfaces, 
         impose it at i=0. The n^i is the normal vec. *)
      Cif == (isSTAR1),
        xc == cxmax1,
      Cif == else,
        xc == cxmax2,
      Cif == end,
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      nv1 == x-xc, (* get normal vec n[a] *)
      nv2 == y,
      nv3 == z,
      nvm == Sqrt[nv1*nv1 + nv2*nv2 + nv3*nv3],
      nv[a] == nv[a]/nvm,
      dSig[a] == dSigma[a],
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      dSigin[a] == dSigmain[a],
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      FSigma == nv[a] (dSig[a] - dSigin[a] OuterSigmaTransitionD1),
      Cinstruction == "} /* endfor */",

      (* set n^i n^j d^2Sigma/(dx^i dx^j) equal at star surfaces, 
         impose it at i=2. The n^i is the normal vec. *)
      Cif == (isSTAR1),
        xc == cxmax1,
      Cif == else,
        xc == cxmax2,
      Cif == end,
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 2){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      nv1 == x-xc, (* get normal vec n[a] *)
      nv2 == y,
      nv3 == z,
      nvm == Sqrt[nv1*nv1 + nv2*nv2 + nv3*nv3],
      nv[a] == nv[a]/nvm,
      ddSig[a,b] == ddSigma[a,b],
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      ddSigin[a,b] == ddSigmain[a,b],
      Cinstruction == "ijk=Index(2,j,k); /* set index to i=2 */",        
      FSigma == nv[a] nv[b] (ddSig[a,b] - ddSigin[a,b] OuterSigmaTransitionD2),
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

      (* set n^i dSigma/dx^i equal at star surfaces, 
         impose it at i=0. The n^i is the normal vec. *)
      Cif == (isSTAR1),
        xc == cxmax1,
      Cif == else,
        xc == cxmax2,
      Cif == end,
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      nv1 == x-xc, (* get normal vec n[a] *)
      nv2 == y,
      nv3 == z,
      nvm == Sqrt[nv1*nv1 + nv2*nv2 + nv3*nv3],
      nv[a] == nv[a]/nvm,
      dlSig[a] == dlSigma[a],
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      dlSigin[a] == dlSigmain[a],
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",        
      FlSigma == nv[a] (dlSig[a] - dlSigin[a] OuterSigmaTransitionD1),
      Cinstruction == "} /* endfor */",

      (* set n^i n^j d^2Sigma/(dx^i dx^j) equal at star surfaces, 
         impose it at i=2. The n^i is the normal vec. *)
      Cif == (isSTAR1),
        xc == cxmax1,
      Cif == else,
        xc == cxmax2,
      Cif == end,
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 2){",
      Cinstruction == "ijk=Index(0,j,k); /* set index to i=0 */",
      nv1 == x-xc, (* get normal vec n[a] *)
      nv2 == y,
      nv3 == z,
      nvm == Sqrt[nv1*nv1 + nv2*nv2 + nv3*nv3],
      nv[a] == nv[a]/nvm,
      ddlSig[a,b] == ddlSigma[a,b],
      Cinstruction == "ijk=Ind_n1n2(n1in-1,j,k,n1in,n2in); /* ijk for other box */",
      ddlSigin[a,b] == ddlSigmain[a,b],
      Cinstruction == "ijk=Index(2,j,k); /* set index to i=2 */",        
      FlSigma == nv[a] nv[b] (ddlSig[a,b] - ddlSigin[a,b] OuterSigmaTransitionD2),
      Cinstruction == "} /* endfor */",

    Cif == end, (* end linear case *)

  Cif == end, (* end case of boxes that touch matter *)


  (* if we get here there is a star surface,
     now treat the case where matter is inside the box *)
  Cif == ( MATTRinside ),
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

      (* we only impose InnerVolIntZero, ... in every sixth box *)
      Cif == (AddInnerVolIntToBC || InnerVolIntZero),
        Cinstruction == "VolAvSigma = 0.0;",
        Cinstruction == "if(isVolAvBox) VolAvSigma =
                         BoxVolumeIntegral(grid->box[bi], index_Sigma);",
      Cif == end,

      Cif == (AddInnerSumToBC || InnerSumZero),
        Cinstruction == "VolAvSigma = 0.0;",
        Cinstruction == "if(isVolAvBox) forallpoints(box, ijk) {",
        Cinstruction == "VolAvSigma += Sigma[ijk];",
        Cinstruction == "} /* endfor */",
      Cif == end,

      (* modify VolAvSigma so that we later impose 
         VolAvSigma=VolAvSigma1/2 *)
      Cif == (isVolAvBox),
        Cinstruction == "VolAvSigma = VolAvSigma - VolAvSigma1;",
      Cif == else,
        Cinstruction == "VolAvSigma = VolAvSigma - VolAvSigma2;",
      Cif == end,

      Cinstruction == "//printf(\"VolAvSigma=%g\\n\",VolAvSigma);",

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
        (* add VolAvSigma=0 to BC *)
        Cif == ( (AddInnerVolIntToBC || AddInnerSumToBC) ),
          FSigma == FSigma + VolAvSigma,
        Cif == end,

      Cinstruction == "} /* end forplane1 */",

      (* impose extra condition in every 6th box *)
      Cinstruction == "if(isVolAvBox) {",
        (* set Sigma to zero at ijk=Index(n1-1,0,0) *)
        Cif == SigmaZeroAtPoint,
          Cinstruction == "ijk=Index(n1-1,0,0);",
          FSigma == Sigma,
        Cif == end,

        (* set VolAvSigma to zero,
           impose it at ijk=Index(n1-1,0,0) in *)
        Cif == (InnerVolIntZero || InnerSumZero),
          Cinstruction == "ijk=Index(n1-1,0,0);",
          FSigma == VolAvSigma,
        Cif == end,
      Cinstruction == "} /* end if(isVolAvBox) */",


    (* linear case: *)
    Cif == else,
      Cinstruction == "FirstDerivsOf_S(box, index_lSigma, index_dlSigma1);",

      Cif == (AddInnerVolIntToBC || InnerVolIntZero),
        Cinstruction == "VolAvlSigma = 0.0;",
        Cinstruction == "if(isVolAvBox) VolAvlSigma =
          BoxVolumeIntegral(grid->box[bi], index_lSigma);",
      Cif == end,

      Cif == (AddInnerSumToBC || InnerSumZero),
        Cinstruction == "VolAvlSigma = 0.0;",
        Cinstruction == "if(isVolAvBox) forallpoints(box, ijk) {",
        Cinstruction == "VolAvlSigma += lSigma[ijk];",
        Cinstruction == "} /* endfor */",
      Cif == end,
      Cinstruction == "//printf(\"VolAvlSigma=%g\\n\",VolAvlSigma);",

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
        alpha == alphaP/Psi,
        alpha2 == alpha alpha,

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

        (* add VolAvlSigma=0 to BC *)
        Cif == ( (AddInnerVolIntToBC || AddInnerSumToBC) ),
          FlSigma == FlSigma + VolAvlSigma,
        Cif == end,

      Cinstruction == "} /* end forplane1 */",

      (* impose extra condition in every 6th box *)
      Cinstruction == "if(isVolAvBox) {",
        (* set Sigma to zero at ijk=Index(n1-1,0,0) *)
        Cif == SigmaZeroAtPoint,
          Cinstruction == "ijk=Index(n1-1,0,0);",
          FlSigma == lSigma,
        Cif == end,

        (* set VolAvSigma to zero, impose it at ijk=Index(n1-1,0,0) *)
        Cif == (InnerVolIntZero || InnerSumZero),
          Cinstruction == "ijk=Index(n1-1,0,0);",
          FlSigma == VolAvlSigma,
        Cif == end,
      Cinstruction == "} /* end if(isVolAvBox) */",


    Cif == end, (* end of nonlin/linear case *)

  Cif == end,

  (* set FSigma and FlSigma such that Sigma remains unchanged on inside *)
  Cif == ( MATTRinside  && KeepInnerSigma ),
    Cinstruction == "int bb;   int star=grid->box[bi]->SIDE;",
    Cinstruction == "forallboxes(grid,bb) {",
      Cinstruction == "tBox *bo = grid->box[bb];",
      Cinstruction == "if(bo->SIDE != star) continue;",
      Cinstruction == "if(bo->MATTR != INSIDE) continue;",
      Cif == nonlin, (* non-linear case *)
        Cinstruction == "double *FSigma_bb = vlldataptr(vlFu, bo, 5);",
        Cinstruction == "/* Zero error ==> do not touch Sigma */";
        Cinstruction == "forallpoints(bo, ijk)
                           FSigma_bb[ijk] = 0.0;",
      Cif == else,   (* linear case *)
        Cinstruction == "double *FlSigma_bb = vlldataptr(vlJdu, bo, 5);",
        Cinstruction == "double *lSigma_bb  = vlldataptr( vldu, bo, 5);",
        Cinstruction == "/* make lin. matrix diagonal ==> no correction */";
        Cinstruction == "forallpoints(bo, ijk)
                           FlSigma_bb[ijk] = lSigma_bb[ijk];",
      Cif == end,
    Cinstruction == "} /* endfor bb */",
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

CFunctionFile = "set_DNSdata_Sigma_BCs.c"

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

  pr["void set_DNSdata_Sigma_BC(tVarList *vlFu, tVarList *vlu, \ 
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
  pr["int AddInnerVolIntToBC = Getv(\"DNSdata_Sigma_surface_BCs\",\"AddInnerVolIntToBC\");\n"];
  pr["int InnerVolIntZero = Getv(\"DNSdata_Sigma_surface_BCs\",\"InnerVolIntZero\");\n"];
  pr["int AddInnerSumToBC = Getv(\"DNSdata_Sigma_surface_BCs\",\"AddInnerSumToBC\");\n"];
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
  pr["double VolAvSigma, VolAvlSigma;\n"];
  pr["double OuterSigmaTransitionD1 = 1.0;\n"];
  pr["double OuterSigmaTransitionD2 = 1.0;\n"];
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
  pr["int isSTAR1     = (box->SIDE == STAR1);\n\n"];
  pr["int MATTRinside = (box->MATTR== INSIDE);\n\n"];
  pr["int MATTRtouch  = (box->MATTR== TOUCH);\n\n"];
  pr["int MATTRaway   = (box->MATTR== AWAY);\n\n"];
  pr["int hasSSURF    = (box->BOUND== SSURF);\n\n"];
  pr["int isXinDom    = (box->CI->dom == box->SIDE - STAR1);\n\n"];
  pr["int isVolAvBox  = (MATTRinside && hasSSURF && isXinDom);\n\n"];

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
