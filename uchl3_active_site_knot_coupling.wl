(* ::Package:: *)

(*
================================================================================
ACTIVE SITE - KNOT COUPLING ANALYSIS
UCH-L3: Does the knot mechanically couple to enzyme function?

Strategy: Even though our mutations didn't fully unknot the protein, we can
still test whether small topological perturbations affect the active site.

UCH-L3 active site (catalytic triad):
- Cys90 (nucleophile)
- His164 (general base)
- Asp179 (stabilize His)

We'll compare:
1. Wild-type (fully knotted)
2. W5P mutant (slightly perturbed knot, 93.2% core conservation)
3. W5P+L6P (most perturbed, 94.9% core conservation)

Question: Even minor knot perturbations - do they distort the active site?
================================================================================
*)

Print[""];
Print["================================================================================"];
Print["ACTIVE SITE - KNOT MECHANICAL COUPLING ANALYSIS"];
Print["UCH-L3 Deubiquitinating Enzyme"];
Print["================================================================================"];
Print[""];

(* ============================================================================
   HELPER FUNCTIONS
   ============================================================================ *)

ExtractCA[biomol_] := Module[{atomCoords},
  atomCoords = BioMoleculeValue[biomol, "AtomCoordinates"];
  Table[
    QuantityMagnitude[atomCoords["A"][[i]][[2]]],
    {i, Length[atomCoords["A"]]}
  ]
];

KabschAlign[mobile_, fixed_] := Module[{
    centroidMobile, centroidFixed,
    mobileCentered, fixedCentered,
    covarianceMatrix, svd, u, v, d, rotation, aligned
  },
  centroidMobile = Mean[mobile];
  centroidFixed = Mean[fixed];
  mobileCentered = Map[# - centroidMobile &, mobile];
  fixedCentered = Map[# - centroidFixed &, fixed];
  covarianceMatrix = Transpose[mobileCentered] . fixedCentered;
  svd = SingularValueDecomposition[covarianceMatrix];
  u = svd[[1]];
  v = svd[[3]];
  d = IdentityMatrix[3];
  If[Det[v . Transpose[u]] < 0, d[[3, 3]] = -1];
  rotation = v . d . Transpose[u];
  aligned = Map[(rotation . #) + centroidFixed &, mobileCentered];
  aligned
];

ComputeDistanceMatrix[coords_] := Module[{n},
  n = Length[coords];
  Table[EuclideanDistance[coords[[i]], coords[[j]]], {i, n}, {j, n}]
];

(* ============================================================================
   LOAD STRUCTURES
   ============================================================================ *)

Print["[STEP 1/5] Loading Structures"];
Print["----------------------------------------------"];

(* Load experimental *)
uchL3Exp = BioMolecule["1XD3"];
wtSeq = BioMoleculeValue[uchL3Exp, "BioSequences"][[1]];
fastaStr = ExportString[wtSeq, "FASTA"];
lines = StringSplit[fastaStr, "\n"];
wtSeqString = StringJoin[Rest[lines]];

Print["  Wild-type sequence loaded (", StringLength[wtSeqString], " residues)"];
Print[""];

(* Fold wild-type *)
Print["  Folding wild-type..."];
wtPred = BioMolecule[wtSeq, BioMoleculeFoldingMethod -> "Local"];
wtPredCA = ExtractCA[wtPred];
Print["  Done"];
Print[""];

(* Create W5P mutant *)
Print["  Creating and folding W5P mutant..."];
w5pSeqString = StringReplacePart[wtSeqString, "P", {5, 5}];
w5pSeq = BioSequence["Peptide", w5pSeqString];
w5pPred = BioMolecule[w5pSeq, BioMoleculeFoldingMethod -> "Local"];
w5pPredCA = ExtractCA[w5pPred];
Print["  Done"];
Print[""];

(* Create W5P+L6P double mutant *)
Print["  Creating and folding W5P+L6P double mutant..."];
doubleMutSeqString = StringReplacePart[
  StringReplacePart[wtSeqString, "P", {5, 5}],
  "P", {6, 6}
];
doubleMutSeq = BioSequence["Peptide", doubleMutSeqString];
doubleMutPred = BioMolecule[doubleMutSeq, BioMoleculeFoldingMethod -> "Local"];
doubleMutPredCA = ExtractCA[doubleMutPred];
Print["  Done"];
Print[""];

(* Align all to wild-type *)
w5pPredCAAligned = KabschAlign[w5pPredCA, wtPredCA];
doubleMutPredCAAligned = KabschAlign[doubleMutPredCA, wtPredCA];

Print["  All structures aligned to wild-type reference frame"];
Print[""];

(* ============================================================================
   DEFINE ACTIVE SITE
   ============================================================================ *)

Print["[STEP 2/5] Defining Active Site"];
Print["----------------------------------------------"];

(* UCH-L3 catalytic triad (approximate positions) *)
cys90 = 90;    (* Nucleophile *)
his164 = 164;  (* General base *)
asp179 = 179;  (* Stabilizer *)

activeSiteResidues = {cys90, his164, asp179};
activeSiteRegion = Range[88, 97];  (* Broader active site pocket *)

Print["  Catalytic triad:"];
Print["    Cys90:  ", StringTake[wtSeqString, {cys90, cys90}], " at position ", cys90];
Print["    His164: ", StringTake[wtSeqString, {his164, his164}], " at position ", his164];
Print["    Asp179: ", StringTake[wtSeqString, {asp179, asp179}], " at position ", asp179];
Print[""];
Print["  Active site region: residues 88-97"];
Print["  Sequence: ", StringTake[wtSeqString, {88, 97}]];
Print[""];

(* ============================================================================
   ACTIVE SITE GEOMETRY ANALYSIS
   ============================================================================ *)

Print["[STEP 3/5] Active Site Geometry Analysis"];
Print["----------------------------------------------"];

(* Extract catalytic triad coordinates *)
wtTriad = wtPredCA[[activeSiteResidues]];
w5pTriad = w5pPredCAAligned[[activeSiteResidues]];
doubleTriad = doubleMutPredCAAligned[[activeSiteResidues]];

(* Compute pairwise distances in catalytic triad *)
Print["  Catalytic triad distances:"];
Print[""];

dCysHisWT = EuclideanDistance[wtTriad[[1]], wtTriad[[2]]];
dCysAspWT = EuclideanDistance[wtTriad[[1]], wtTriad[[3]]];
dHisAspWT = EuclideanDistance[wtTriad[[2]], wtTriad[[3]]];

Print["  Wild-type:"];
Print["    Cys90-His164: ", Round[dCysHisWT, 0.1], " Angstroms"];
Print["    Cys90-Asp179: ", Round[dCysAspWT, 0.1], " Angstroms"];
Print["    His164-Asp179: ", Round[dHisAspWT, 0.1], " Angstroms"];
Print[""];

dCysHisW5P = EuclideanDistance[w5pTriad[[1]], w5pTriad[[2]]];
dCysAspW5P = EuclideanDistance[w5pTriad[[1]], w5pTriad[[3]]];
dHisAspW5P = EuclideanDistance[w5pTriad[[2]], w5pTriad[[3]]];

Print["  W5P mutant:"];
Print["    Cys90-His164: ", Round[dCysHisW5P, 0.1], " A (Delta: ",
      Round[dCysHisW5P - dCysHisWT, 0.2], " A)"];
Print["    Cys90-Asp179: ", Round[dCysAspW5P, 0.1], " A (Delta: ",
      Round[dCysAspW5P - dCysAspWT, 0.2], " A)"];
Print["    His164-Asp179: ", Round[dHisAspW5P, 0.1], " A (Delta: ",
      Round[dHisAspW5P - dHisAspWT, 0.2], " A)"];
Print[""];

dCysHisDouble = EuclideanDistance[doubleTriad[[1]], doubleTriad[[2]]];
dCysAspDouble = EuclideanDistance[doubleTriad[[1]], doubleTriad[[3]]];
dHisAspDouble = EuclideanDistance[doubleTriad[[2]], doubleTriad[[3]]];

Print["  W5P+L6P double mutant:"];
Print["    Cys90-His164: ", Round[dCysHisDouble, 0.1], " A (Delta: ",
      Round[dCysHisDouble - dCysHisWT, 0.2], " A)"];
Print["    Cys90-Asp179: ", Round[dCysAspDouble, 0.1], " A (Delta: ",
      Round[dCysAspDouble - dCysAspWT, 0.2], " A)"];
Print["    His164-Asp179: ", Round[dHisAspDouble, 0.1], " A (Delta: ",
      Round[dHisAspDouble - dHisAspWT, 0.2], " A)"];
Print[""];

(* ============================================================================
   ACTIVE SITE RMSD
   ============================================================================ *)

Print["[STEP 4/5] Active Site RMSD Analysis"];
Print["----------------------------------------------"];

(* Active site region RMSD *)
wtActiveSite = wtPredCA[[activeSiteRegion]];
w5pActiveSite = w5pPredCAAligned[[activeSiteRegion]];
doubleActiveSite = doubleMutPredCAAligned[[activeSiteRegion]];

activeSiteRMSDW5P = Sqrt[Mean[Map[
  SquaredEuclideanDistance[#[[1]], #[[2]]]&,
  Transpose[{wtActiveSite, w5pActiveSite}]
]]];

activeSiteRMSDDouble = Sqrt[Mean[Map[
  SquaredEuclideanDistance[#[[1]], #[[2]]]&,
  Transpose[{wtActiveSite, doubleActiveSite}]
]]];

Print["  Active site region (residues 88-97) RMSD:"];
Print["    W5P mutant: ", Round[activeSiteRMSDW5P, 0.02], " Angstroms"];
Print["    W5P+L6P mutant: ", Round[activeSiteRMSDDouble, 0.02], " Angstroms"];
Print[""];

(* ============================================================================
   KNOT-TO-ACTIVE-SITE DISTANCE ANALYSIS
   ============================================================================ *)

Print["[STEP 5/5] Knot-to-Active-Site Coupling Analysis"];
Print["----------------------------------------------"];

(* Measure distance from threading region to active site *)
threadingRegion = Range[5, 8];
knotLoop = Range[147, 161];  (* Loop that gets threaded *)

wtThreading = Mean[wtPredCA[[threadingRegion]]];
wtKnotLoop = Mean[wtPredCA[[knotLoop]]];
wtActiveSiteCentroid = Mean[wtActiveSite];

distThreadToActive = EuclideanDistance[wtThreading, wtActiveSiteCentroid];
distLoopToActive = EuclideanDistance[wtKnotLoop, wtActiveSiteCentroid];
distThreadToLoop = EuclideanDistance[wtThreading, wtKnotLoop];

Print["  Spatial relationships in wild-type:"];
Print["    Threading region -> Active site: ", Round[distThreadToActive, 0.1], " A"];
Print["    Knot loop -> Active site: ", Round[distLoopToActive, 0.1], " A"];
Print["    Threading -> Knot loop: ", Round[distThreadToLoop, 0.1], " A"];
Print[""];

(* Check if His164 is in the loop *)
Print["  Note: His164 (active site) is at position ", his164];
Print["        Knot loop spans residues 147-161"];
If[his164 >= 147 && his164 <= 161,
  Print["        >>> His164 IS PART OF THE KNOT LOOP!"];
  Print["        >>> Direct mechanical coupling expected!"],
  Print["        >>> His164 is near but not in the loop"];
  Print["        >>> Coupling may be indirect"]
];
Print[""];

(* ============================================================================
   FINAL ASSESSMENT
   ============================================================================ *)

Print["================================================================================"];
Print["ACTIVE SITE - KNOT COUPLING ASSESSMENT"];
Print["================================================================================"];
Print[""];

maxTriadDistW5P = Max[{
  Abs[dCysHisW5P - dCysHisWT],
  Abs[dCysAspW5P - dCysAspWT],
  Abs[dHisAspW5P - dHisAspWT]
}];

maxTriadDistDouble = Max[{
  Abs[dCysHisDouble - dCysHisWT],
  Abs[dCysAspDouble - dCysAspWT],
  Abs[dHisAspDouble - dHisAspWT]
}];

Print["  Maximum catalytic triad distortion:"];
Print["    W5P: ", Round[maxTriadDistW5P, 0.02], " Angstroms"];
Print["    W5P+L6P: ", Round[maxTriadDistDouble, 0.02], " Angstroms"];
Print[""];

Print["  Active site RMSD:"];
Print["    W5P: ", Round[activeSiteRMSDW5P, 0.02], " Angstroms"];
Print["    W5P+L6P: ", Round[activeSiteRMSDDouble, 0.02], " Angstroms"];
Print[""];

If[maxTriadDistDouble < 0.5 && activeSiteRMSDDouble < 0.5,
  Print["  VERDICT: NO SIGNIFICANT COUPLING"];
  Print["    Active site geometry preserved despite knot perturbations"];
  Print["    Knot and active site are structurally independent"];
  couplingVerdict = "NO_COUPLING",
  If[maxTriadDistDouble < 1.0 && activeSiteRMSDDouble < 1.0,
    Print["  VERDICT: WEAK COUPLING"];
    Print["    Minor active site changes with knot perturbations"];
    Print["    Some structural communication but not strong"];
    couplingVerdict = "WEAK_COUPLING",
    Print["  VERDICT: STRONG COUPLING"];
    Print["    Significant active site distortion from knot perturbations"];
    Print["    Knot is mechanically coupled to enzyme function"];
    couplingVerdict = "STRONG_COUPLING"
  ]
];
Print[""];

Print["BIOLOGICAL INTERPRETATION:"];
If[couplingVerdict == "NO_COUPLING",
  Print["  The knot is likely a STRUCTURAL SCAFFOLD"];
  Print["  - Provides global stability"];
  Print["  - Not directly involved in catalysis"];
  Print["  - Could potentially be engineered away without losing function"],
  If[couplingVerdict == "WEAK_COUPLING",
    Print["  The knot provides INDIRECT FUNCTIONAL SUPPORT"];
    Print["  - May affect dynamics or substrate binding"];
    Print["  - Removes some conformational degrees of freedom"];
    Print["  - Function might be impaired without it"],
    Print["  The knot is ESSENTIAL FOR FUNCTION"];
    Print["  - Directly controls active site geometry"];
    Print["  - Cannot be removed without destroying enzyme activity"];
    Print["  - This is why evolution preserved the knot"]
  ]
];
Print[""];

(* Export results *)
results = <|
  "Catalytic_Triad_Cys_His_WT" -> dCysHisWT,
  "Catalytic_Triad_Cys_Asp_WT" -> dCysAspWT,
  "Catalytic_Triad_His_Asp_WT" -> dHisAspWT,
  "Max_Triad_Distortion_W5P" -> maxTriadDistW5P,
  "Max_Triad_Distortion_W5P_L6P" -> maxTriadDistDouble,
  "ActiveSite_RMSD_W5P" -> activeSiteRMSDW5P,
  "ActiveSite_RMSD_W5P_L6P" -> activeSiteRMSDDouble,
  "Threading_To_ActiveSite_Distance" -> distThreadToActive,
  "KnotLoop_To_ActiveSite_Distance" -> distLoopToActive,
  "Coupling_Verdict" -> couplingVerdict
|>;

Export["D:\\uchl3_active_site_coupling_results.json", results, "JSON"];
Print["Results exported to D:\\uchl3_active_site_coupling_results.json"];
Print[""];

Print["================================================================================"];
Print["ACTIVE SITE - KNOT COUPLING ANALYSIS COMPLETE"];
Print["================================================================================"];
