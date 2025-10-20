Print[""];
Print["================================================================================"];
Print["UCH-L1 GEOMETRIC VALIDATION"];
Print["Verify predicted structure accuracy independent of writhe"];
Print["================================================================================"];
Print[""];

ExtractCA[biomol_] := Module[{atomCoords},
  atomCoords = BioMoleculeValue[biomol, "AtomCoordinates"];
  Table[QuantityMagnitude[atomCoords["A"][[i]][[2]]], {i, Length[atomCoords["A"]]}]
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

Print["[STEP 1/5] Loading Experimental UCH-L1"];
Print["----------------------------------------------"];
Print[""];

uchL1Exp = BioMolecule["2ETL"];
uchL1ExpCA = ExtractCA[uchL1Exp];
Print["  PDB: 2ETL"];
Print["  Residues: ", Length[uchL1ExpCA]];
Print[""];

Print["[STEP 2/5] Loading Predicted UCH-L1"];
Print["----------------------------------------------"];
Print[""];

uchL1PredCA = Import["D:\\uchl1_predicted_CA.mx"];
Print["  Predicted residues: ", Length[uchL1PredCA]];
Print[""];

If[Length[uchL1ExpCA] != Length[uchL1PredCA],
  Print["  ERROR: Length mismatch"];
  Abort[];
];

Print["[STEP 3/5] Kabsch Alignment"];
Print["----------------------------------------------"];
Print[""];

Print["  Aligning predicted to experimental structure..."];
uchL1PredAligned = KabschAlign[uchL1PredCA, uchL1ExpCA];
Print["  Done"];
Print[""];

Print["[STEP 4/5] Computing RMSD"];
Print["----------------------------------------------"];
Print[""];

rmsd = Sqrt[Mean[Map[
  SquaredEuclideanDistance[#[[1]], #[[2]]]&,
  Transpose[{uchL1ExpCA, uchL1PredAligned}]
]]];

Print["  Aligned RMSD: ", Round[rmsd, 0.02], " A"];
Print[""];

Print["  RMSD Interpretation:"];
Print["    < 2.0 A: Excellent (near-experimental accuracy)"];
Print["    < 4.0 A: Good (correct overall fold)"];
Print["    < 8.0 A: Moderate (partial success)"];
Print["    > 8.0 A: Poor (incorrect fold)"];
Print[""];

If[rmsd < 2.0,
  Print["  Assessment: EXCELLENT"];
  Print["    Predicted structure is geometrically accurate"];
  Print["    Writhe elevation (0.16 -> 0.54) is ANOMALOUS"];
  geometricVerdict = "excellent",
  If[rmsd < 4.0,
    Print["  Assessment: GOOD"];
    Print["    Overall fold correct"];
    Print["    Some local deviations may affect writhe"];
    geometricVerdict = "good",
    Print["  Assessment: POOR"];
    Print["    Prediction inaccurate"];
    Print["    Writhe elevation may reflect misprediction"];
    geometricVerdict = "poor"
  ]
];
Print[""];

Print["[STEP 5/5] Contact Pattern Comparison"];
Print["----------------------------------------------"];
Print[""];

Print["  Computing distance matrices..."];
expDistMatrix = ComputeDistanceMatrix[uchL1ExpCA];
predDistMatrix = ComputeDistanceMatrix[uchL1PredAligned];
Print["  Done"];
Print[""];

Print["  Threading contacts (>50 residues apart, <8A):"];
expThreading = CountThreadingContacts[expDistMatrix, 50, 8.0];
predThreading = CountThreadingContacts[predDistMatrix, 50, 8.0];

Print["    Experimental: ", expThreading];
Print["    Predicted: ", predThreading];
Print["    Conservation: ", Round[100.0 * predThreading / expThreading, 0.1], "%"];
Print[""];

Print["  Per-residue distance error distribution:"];
perResidueErrors = Table[
  Sqrt[SquaredEuclideanDistance[uchL1ExpCA[[i]], uchL1PredAligned[[i]]]],
  {i, Length[uchL1ExpCA]}
];

Print["    Mean: ", Round[Mean[perResidueErrors], 0.02], " A"];
Print["    Median: ", Round[Median[perResidueErrors], 0.02], " A"];
Print["    Max: ", Round[Max[perResidueErrors], 0.02], " A"];
Print["    Residues >5A error: ", Count[perResidueErrors, x_ /; x > 5.0]];
Print[""];

Print["================================================================================"];
Print["DIAGNOSTIC ANALYSIS"];
Print["================================================================================"];
Print[""];

Print["For comparison - UCH-L3 (knotted):"];
Print["  Predicted RMSD: 0.84 A"];
Print["  Threading conservation: 101.9%"];
Print["  Predicted writhe: 0.56"];
Print["  Experimental writhe: 0.53"];
Print[""];

Print["UCH-L1 (unknotted):"];
Print["  Predicted RMSD: ", Round[rmsd, 0.02], " A"];
Print["  Threading conservation: ", Round[100.0 * predThreading / expThreading, 0.1], "%"];
Print["  Predicted writhe: 0.54"];
Print["  Experimental writhe: 0.16"];
Print[""];

Print["================================================================================"];
Print["CONCLUSION"];
Print["================================================================================"];
Print[""];

If[geometricVerdict == "excellent" || geometricVerdict == "good",
  Print["FINDING: Geometric accuracy WITHOUT topological accuracy"];
  Print[""];
  Print["  UCH-L1 prediction is geometrically accurate (RMSD ", Round[rmsd, 0.02], " A)"];
  Print["  But writhe is severely elevated (0.16 -> 0.54)"];
  Print[""];
  Print["  This indicates ONE of the following:"];
  Print[""];
  Print["  Option 1: Writhe calculation artifact"];
  Print["    - Discrete curve sampling creates false crossings"];
  Print["    - Writhe integral is sensitive to numerical noise"];
  Print["    - C-alpha trace approximation inadequate for writhe"];
  Print[""];
  Print["  Option 2: Subtle backbone geometry differences"];
  Print["    - RMSD measures average deviation");
  Print["    - Writhe sensitive to global curve path"];
  Print["    - Small local errors accumulate in writhe integral"];
  Print[""];
  Print["  Option 3: ESMFold produces systematically 'tighter' structures"];
  Print["    - Backbone slightly more compact than reality");
  Print["    - Creates apparent self-crossings not in experiment"];
  Print["    - Bias affects all predictions similarly"];
  Print[""];
  Print["  CRITICAL: This affects UCH-L3 interpretation");
  Print["    - UCH-L3 predicted writhe (0.56) may have same artifact");
  Print["    - Cannot use writhe to validate topology learning");
  Print["    - Geometric accuracy != Topological accuracy");
  ,
  Print["FINDING: Poor geometric prediction"];
  Print["  UCH-L1 prediction is inaccurate (RMSD ", Round[rmsd, 0.02], " A)"];
  Print["  Writhe elevation may simply reflect misprediction");
  Print["  Cannot draw conclusions about topological learning");
];
Print[""];

Print["================================================================================"];
Print["EXPORT"];
Print["================================================================================"];
Print[""];

validation = <|
  "UCH_L1_RMSD" -> rmsd,
  "UCH_L1_Exp_Threading" -> expThreading,
  "UCH_L1_Pred_Threading" -> predThreading,
  "Threading_Conservation" -> predThreading / expThreading,
  "Geometric_Verdict" -> geometricVerdict,
  "Per_Residue_Errors" -> perResidueErrors,
  "Mean_Error" -> Mean[perResidueErrors],
  "Max_Error" -> Max[perResidueErrors]
|>;

Export["D:\\uchl1_geometric_validation.json", validation, "JSON"];
Print["  uchl1_geometric_validation.json"];

Export["D:\\uchl1_predicted_aligned_CA.mx", uchL1PredAligned];
Print["  uchl1_predicted_aligned_CA.mx"];

Print[""];
Print["================================================================================"];
Print["VALIDATION COMPLETE"];
Print["================================================================================"];
Print[""];
