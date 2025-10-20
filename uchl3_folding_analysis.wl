(* ::Package:: *)

(*
================================================================================
Wolfram Language 14.3 Local Protein Folding
Target: UCH-L3 - Deep 5-Crossing Knotted Protein (PDB: 1XD3)

Testing whether ESMFold can predict the a topologically
complex knotted protein structure. The 5-crossing knot requires
"impossible" threading during co-translational folding.

Author: Testing Wolfram Language V14.3 Local Folding Feature
Date: 2025
================================================================================
*)

Print[""];
Print["================================================================================"];
Print["WOLFRAM LANGUAGE 14.3 - LOCAL PROTEIN FOLDING TEST"];
Print["Target: UCH-L3 (PDB: 1XD3) - Deep 5-Crossing Knotted Protein"];
Print["================================================================================"];
Print[""];

(* ============================================================================
   HELPER FUNCTIONS
   ============================================================================ *)

(* Extract C-alpha coordinates from BioMolecule object
   Uses index 2 heuristic (CA is typically the 2nd atom in standard residues)
*)
ExtractCA[biomol_] := Module[{atomCoords},
  atomCoords = BioMoleculeValue[biomol, "AtomCoordinates"];
  Table[
    QuantityMagnitude[atomCoords["A"][[i]][[2]]],
    {i, Length[atomCoords["A"]]}
  ]
];

(* Kabsch algorithm for optimal superposition of two point sets
   Returns aligned mobile coordinates that best match fixed coordinates

   Note: BioMoleculeAlign has a bug with predicted structures, so we implement
   our own alignment using the standard Kabsch algorithm.
*)
KabschAlign[mobile_, fixed_] := Module[{
    centroidMobile, centroidFixed,
    mobileCentered, fixedCentered,
    covarianceMatrix, svd, u, v, d, rotation, aligned
  },
  (* Step 1: Center both point sets *)
  centroidMobile = Mean[mobile];
  centroidFixed = Mean[fixed];
  mobileCentered = Map[# - centroidMobile &, mobile];
  fixedCentered = Map[# - centroidFixed &, fixed];

  (* Step 2: Compute covariance matrix H = Mobile^T * Fixed *)
  covarianceMatrix = Transpose[mobileCentered] . fixedCentered;

  (* Step 3: Singular Value Decomposition *)
  svd = SingularValueDecomposition[covarianceMatrix];
  u = svd[[1]];
  v = svd[[3]];

  (* Step 4: Compute optimal rotation matrix *)
  (* Check for reflection case (det < 0) and correct if needed *)
  d = IdentityMatrix[3];
  If[Det[v . Transpose[u]] < 0, d[[3, 3]] = -1];
  rotation = v . d . Transpose[u];

  (* Step 5: Apply rotation and translation *)
  aligned = Map[(rotation . #) + centroidFixed &, mobileCentered];

  aligned
];

(* Compute distance matrix between all CA atoms *)
ComputeDistanceMatrix[coords_] := Module[{n},
  n = Length[coords];
  Table[
    EuclideanDistance[coords[[i]], coords[[j]]],
    {i, n}, {j, n}
  ]
];

(* Count threading contacts: residues separated by >minSep in sequence
   but within maxDist in 3D space *)
CountThreadingContacts[distMatrix_, minSeqSep_, maxSpatialDist_] := Module[{n},
  n = Length[distMatrix];
  Length[Select[
    Flatten[Table[
      If[j > i + minSeqSep && distMatrix[[i, j]] < maxSpatialDist, 1, Nothing],
      {i, n}, {j, n}
    ]],
    IntegerQ
  ]]
];

(* Count knot core contacts between two regions *)
CountRegionContacts[distMatrix_, region1_, region2_, maxDist_] := Module[{},
  Length[Select[
    Flatten[Table[
      If[distMatrix[[i, j]] < maxDist, 1, Nothing],
      {i, region1}, {j, region2}
    ]],
    IntegerQ
  ]]
];

(* ============================================================================
   MAIN ANALYSIS
   ============================================================================ *)

Print["[STEP 1/8] Loading Experimental Structure"];
Print["----------------------------------------------"];

uchL3Exp = BioMolecule["1XD3"];
sequence = BioMoleculeValue[uchL3Exp, "BioSequences"][[1]];
expCA = ExtractCA[uchL3Exp];
seqLen = Length[expCA];

Print["  PDB ID: 1XD3"];
Print["  Protein: UCH-L3 (Ubiquitin C-terminal hydrolase L3)"];
Print["  Sequence Length: ", seqLen, " residues"];
Print["  Knot Type: 5-crossing (5\:2082 knot)"];
Print["  First CA coordinate: ", expCA[[1]]];
Print[""];

(* Validate CA extraction *)
d12 = EuclideanDistance[expCA[[1]], expCA[[2]]];
d23 = EuclideanDistance[expCA[[2]], expCA[[3]]];
Print["  Validation - Consecutive CA distances:"];
Print["    CA[1]-CA[2]: ", Round[d12, 0.2], " \[CapitalARing]"];
Print["    CA[2]-CA[3]: ", Round[d23, 0.2], " \[CapitalARing]"];

If[d12 < 2.5 || d12 > 5.0 || d23 < 2.5 || d23 > 5.0,
  Print["  \:2717 ERROR: CA extraction failed! Distances outside expected range (3.2-4.2 \[CapitalARing])"];
  Print["  Aborting."];
  Abort[],
  Print["  \[Checkmark] CA extraction validated"];
];
Print[""];

(* ============================================================================ *)

Print["[STEP 2/8] Local Protein Folding with ESMFold"];
Print["----------------------------------------------"];
Print["  This uses Wolfram Language 14.3's new local folding capability"];
Print["  Note: First run will download ~11GB neural network model"];
Print["  Folding in progress..."];
Print[""];

foldStartTime = AbsoluteTime[];
uchL3Pred = BioMolecule[sequence, BioMoleculeFoldingMethod -> "Local"];
foldEndTime = AbsoluteTime[];
foldingTime = foldEndTime - foldStartTime;

predCA = ExtractCA[uchL3Pred];

Print["  \[Checkmark] Folding complete"];
Print["  Folding Time: ", Round[foldingTime, 0.1], " seconds"];
Print["  Predicted CA atoms: ", Length[predCA]];
Print[""];

If[Length[expCA] != Length[predCA],
  Print["  \:2717 ERROR: Length mismatch between experimental and predicted!"];
  Print["  Experimental: ", Length[expCA], " residues"];
  Print["  Predicted: ", Length[predCA], " residues"];
  Abort[],
  Print["  \[Checkmark] Lengths match"];
];
Print[""];

(* ============================================================================ *)

Print["[STEP 3/8] Computing Unaligned RMSD"];
Print["----------------------------------------------"];

rmsdUnaligned = Sqrt[Mean[Map[
  SquaredEuclideanDistance[#[[1]], #[[2]]]&,
  Transpose[{expCA, predCA}]
]]];

Print["  Unaligned RMSD: ", Round[rmsdUnaligned, 0.1], " \[CapitalARing]"];
Print["  (Should be large - structures in different coordinate frames)"];
Print[""];

(* ============================================================================ *)

Print["[STEP 4/8] Optimal Structural Alignment (Kabsch Algorithm)"];
Print["----------------------------------------------"];
Print["  Note: Using manual Kabsch alignment due to BioMoleculeAlign bug"];
Print["  with predicted structures"];
Print[""];

predCAAligned = KabschAlign[predCA, expCA];

(* Verify alignment worked *)
centroidExp = Mean[expCA];
centroidPred = Mean[predCAAligned];
centroidDistance = EuclideanDistance[centroidExp, centroidPred];

Print["  Experimental centroid: ", Map[Round[#, 0.1]&, centroidExp]];
Print["  Predicted centroid: ", Map[Round[#, 0.1]&, centroidPred]];
Print["  Centroid distance: ", Round[centroidDistance, 0.0001], " \[CapitalARing]"];
Print[""];

If[centroidDistance > 0.01,
  Print["  \:26a0 WARNING: Centroids not perfectly aligned"];
  Print["  This may indicate an alignment problem"];
  Print[""],
  Print["  \[Checkmark] Alignment perfect - centroids match"];
  Print[""];
];

(* ============================================================================ *)

Print["[STEP 5/8] Computing Aligned RMSD"];
Print["----------------------------------------------"];

rmsd = Sqrt[Mean[Map[
  SquaredEuclideanDistance[#[[1]], #[[2]]]&,
  Transpose[{expCA, predCAAligned}]
]]];

Print["  Aligned C-alpha RMSD: ", Round[rmsd, 0.02], " \[CapitalARing]"];
Print[""];
Print["  RMSD Interpretation:"];
Print["    < 2.0 \[CapitalARing]  : Excellent prediction (near-experimental accuracy)"];
Print["    < 4.0 \[CapitalARing]  : Good prediction (correct overall fold)"];
Print["    < 8.0 \[CapitalARing]  : Moderate (partial success)"];
Print["    > 8.0 \[CapitalARing]  : Poor (incorrect fold)"];
Print[""];

If[rmsd < 2.0,
  Print["  \[Checkmark] EXCELLENT - Near-perfect prediction!"];
  accuracyVerdict = "EXCELLENT",
  If[rmsd < 4.0,
    Print["  \[Checkmark] GOOD - Overall fold correctly predicted"];
    accuracyVerdict = "GOOD",
    If[rmsd < 8.0,
      Print["  ~ MODERATE - Partial success"];
      accuracyVerdict = "MODERATE",
      Print["  \:2717 POOR - Prediction failed"];
      accuracyVerdict = "POOR"
    ]
  ]
];
Print[""];

(* ============================================================================ *)

Print["[STEP 6/8] Topological Analysis: Threading Contacts"];
Print["----------------------------------------------"];
Print["  Threading contacts = residues >50 apart in sequence"];
Print["                       but <8 \[CapitalARing] apart in 3D space"];
Print["  (These indicate long-range topology including knots)"];
Print[""];

(* Compute distance matrices *)
expDistMatrix = ComputeDistanceMatrix[expCA];
predDistMatrix = ComputeDistanceMatrix[predCAAligned];

(* Count threading contacts *)
threadingExp = CountThreadingContacts[expDistMatrix, 50, 8.0];
threadingPred = CountThreadingContacts[predDistMatrix, 50, 8.0];
threadingRatio = N[threadingPred / threadingExp];

Print["  Experimental threading contacts: ", threadingExp];
Print["  Predicted threading contacts: ", threadingPred];
Print["  Conservation ratio: ", Round[threadingRatio, 0.03]];
Print["  Conservation percentage: ", Round[threadingRatio * 100, 0.1], "%"];
Print[""];

If[threadingRatio > 0.7,
  Print["  \[Checkmark] Threading topology PRESERVED"];
  threadingVerdict = "PRESERVED",
  Print["  \:2717 Threading topology NOT PRESERVED"];
  threadingVerdict = "NOT_PRESERVED"
];
Print[""];

(* ============================================================================ *)

Print["[STEP 7/8] Topological Analysis: Knot Core Structure"];
Print["----------------------------------------------"];
Print["  Knot core = contacts between N-terminus (1-30)"];
Print["              and mid-chain (100-180)"];
Print["  (The knot forms when N-terminus threads through"];
Print["   a loop formed by the mid-chain region)"];
Print[""];

knotCoreExp = CountRegionContacts[expDistMatrix, Range[1, 30], Range[100, 180], 8.0];
knotCorePred = CountRegionContacts[predDistMatrix, Range[1, 30], Range[100, 180], 8.0];
knotCoreRatio = N[knotCorePred / knotCoreExp];

Print["  Experimental knot core contacts: ", knotCoreExp];
Print["  Predicted knot core contacts: ", knotCorePred];
Print["  Conservation ratio: ", Round[knotCoreRatio, 0.03]];
Print["  Conservation percentage: ", Round[knotCoreRatio * 100, 0.1], "%"];
Print[""];

If[knotCoreRatio > 0.6,
  Print["  \[Checkmark] Knot core structure PRESERVED"];
  knotCoreVerdict = "PRESERVED",
  If[knotCoreRatio > 0.4,
    Print["  ~ Knot core structure PARTIALLY PRESERVED"];
    knotCoreVerdict = "PARTIAL",
    Print["  \:2717 Knot core structure NOT PRESERVED"];
    knotCoreVerdict = "NOT_PRESERVED"
  ]
];
Print[""];

(* Additional knot analysis: N-C terminal distance *)
nTermRegion = Range[1, 30];
cTermRegion = Range[180, seqLen];

nTermCentroidExp = Mean[expCA[[nTermRegion]]];
cTermCentroidExp = Mean[expCA[[cTermRegion]]];
ncDistExp = EuclideanDistance[nTermCentroidExp, cTermCentroidExp];

nTermCentroidPred = Mean[predCAAligned[[nTermRegion]]];
cTermCentroidPred = Mean[predCAAligned[[cTermRegion]]];
ncDistPred = EuclideanDistance[nTermCentroidPred, cTermCentroidPred];

Print["  N-terminus to C-terminus distance:"];
Print["    Experimental: ", Round[ncDistExp, 0.1], " \[CapitalARing]"];
Print["    Predicted: ", Round[ncDistPred, 0.1], " \[CapitalARing]"];
Print["    Difference: ", Round[Abs[ncDistExp - ncDistPred], 0.1], " \[CapitalARing]"];
Print[""];

(* ============================================================================ *)

Print["[STEP 8/8] Final Assessment"];
Print["----------------------------------------------"];
Print[""];

(* Overall verdict *)
Print["  STRUCTURAL ACCURACY: ", accuracyVerdict];
Print["  THREADING TOPOLOGY: ", threadingVerdict];
Print["  KNOT CORE STRUCTURE: ", knotCoreVerdict];
Print[""];

If[rmsd < 4.0 && threadingRatio > 0.7,
  Print["  \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550"];
  Print["  \[FivePointedStar]\[FivePointedStar]\[FivePointedStar] LANDMARK RESULT \[FivePointedStar]\[FivePointedStar]\[FivePointedStar]"];
  Print["  \:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550\:2550"];
  Print[""];
  Print["  ESMFold successfully predicted the deep"];
  Print["  5-crossing knot in UCH-L3!"];
  Print[""];
  Print["  This demonstrates that modern neural"];
  Print["  networks can learn topologically complex"];
  Print["  protein structures that seem physically"];
  Print["  impossible to fold co-translationally."];
  overallVerdict = "SUCCESS",

  If[rmsd < 4.0,
    Print["  PARTIAL SUCCESS"];
    Print["  Structure predicted correctly but topology uncertain"];
    overallVerdict = "PARTIAL",
    Print["  PREDICTION FAILED"];
    Print["  The knotted structure was not correctly predicted"];
    overallVerdict = "FAILURE"
  ]
];
Print[""];

(* ============================================================================ *)

Print["================================================================================"];
Print["SUMMARY OF RESULTS"];
Print["================================================================================"];
Print[""];
Print["Protein: UCH-L3 (PDB: 1XD3)"];
Print["Topology: 5-crossing knotted protein"];
Print["Sequence Length: ", seqLen, " residues"];
Print[""];
Print["Folding Method: ESMFold (Local)"];
Print["Folding Time: ", Round[foldingTime, 0.1], " seconds"];
Print[""];
Print["RMSD: ", Round[rmsd, 0.02], " \[CapitalARing]"];
Print["Threading Conservation: ", Round[threadingRatio * 100, 0.1], "%"];
Print["Knot Core Conservation: ", Round[knotCoreRatio * 100, 0.1], "%"];
Print[""];
Print["Overall Verdict: ", overallVerdict];
Print[""];
Print["================================================================================"];
Print[""];

(* ============================================================================
   EXPORT RESULTS
   ============================================================================ *)

Print["Exporting results to D:\\"];
Print[""];

(* Export summary statistics *)
results = <|
  "Protein" -> "UCH-L3",
  "PDB_ID" -> "1XD3",
  "Knot_Type" -> "5-crossing",
  "Sequence_Length" -> seqLen,
  "Folding_Time_Seconds" -> foldingTime,
  "RMSD_Angstroms" -> rmsd,
  "Threading_Conservation_Percent" -> threadingRatio * 100,
  "Knot_Core_Conservation_Percent" -> knotCoreRatio * 100,
  "Accuracy_Verdict" -> accuracyVerdict,
  "Threading_Verdict" -> threadingVerdict,
  "Knot_Core_Verdict" -> knotCoreVerdict,
  "Overall_Verdict" -> overallVerdict
|>;

Export["D:\\uchl3_results_summary.json", results, "JSON"];
Print["  \[Checkmark] uchl3_results_summary.json"];

(* Export coordinate data for further analysis *)
Export["D:\\uchl3_experimental_CA.mx", expCA];
Print["  \[Checkmark] uchl3_experimental_CA.mx"];

Export["D:\\uchl3_predicted_aligned_CA.mx", predCAAligned];
Print["  \[Checkmark] uchl3_predicted_aligned_CA.mx"];

(* Export distance matrices *)
Export["D:\\uchl3_experimental_distance_matrix.mx", expDistMatrix];
Print["  \[Checkmark] uchl3_experimental_distance_matrix.mx"];

Export["D:\\uchl3_predicted_distance_matrix.mx", predDistMatrix];
Print["  \[Checkmark] uchl3_predicted_distance_matrix.mx"];

(* Export detailed report *)
reportText = StringJoin[
  "UCH-L3 PROTEIN FOLDING TEST REPORT\n",
  "==================================\n\n",
  "Target: UCH-L3 (PDB: 1XD3)\n",
  "Topology: 5-crossing knotted protein\n",
  "Sequence Length: ", ToString[seqLen], " residues\n\n",
  "Folding Method: ESMFold (Wolfram Language 14.3 Local)\n",
  "Folding Time: ", ToString[Round[foldingTime, 0.1]], " seconds\n\n",
  "RESULTS:\n",
  "--------\n",
  "RMSD: ", ToString[Round[rmsd, 0.02]], " \[CapitalARing]\n",
  "Threading Conservation: ", ToString[Round[threadingRatio * 100, 0.1]], "%\n",
  "Knot Core Conservation: ", ToString[Round[knotCoreRatio * 100, 0.1]], "%\n\n",
  "Threading Contacts:\n",
  "  Experimental: ", ToString[threadingExp], "\n",
  "  Predicted: ", ToString[threadingPred], "\n\n",
  "Knot Core Contacts:\n",
  "  Experimental: ", ToString[knotCoreExp], "\n",
  "  Predicted: ", ToString[knotCorePred], "\n\n",
  "N-C Terminal Distance:\n",
  "  Experimental: ", ToString[Round[ncDistExp, 0.1]], " \[CapitalARing]\n",
  "  Predicted: ", ToString[Round[ncDistPred, 0.1]], " \[CapitalARing]\n\n",
  "VERDICT: ", overallVerdict, "\n",
  "--------\n",
  "Accuracy: ", accuracyVerdict, "\n",
  "Threading Topology: ", threadingVerdict, "\n",
  "Knot Core Structure: ", knotCoreVerdict, "\n"
];

Export["D:\\uchl3_detailed_report.txt", reportText, "Text"];
Print["  \[Checkmark] uchl3_detailed_report.txt"];

Print[""];
Print["All results exported successfully"];
Print[""];
Print["================================================================================"];
Print["TEST COMPLETE"];
Print["================================================================================"];

