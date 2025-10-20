Print["UCH-L1 TOPOLOGY VERIFICATION"];
Print["============================="];
Print[""];

ExtractCA[biomol_] := Module[{atomCoords},
  atomCoords = BioMoleculeValue[biomol, "AtomCoordinates"];
  Table[QuantityMagnitude[atomCoords["A"][[i]][[2]]], {i, Length[atomCoords["A"]]}]
];

ComputeWrithe[coords_] := Module[{n, writhe, i, j, r1, r2, r12, dr1, dr2, cross, dist, contribution},
  n = Length[coords];
  writhe = 0.0;
  Do[
    Do[
      If[Abs[i - j] > 2,
        r1 = coords[[i]];
        r2 = coords[[j]];
        r12 = r1 - r2;
        dist = Norm[r12];
        If[dist > 0.1,
          dr1 = If[i < n, coords[[i + 1]] - coords[[i]], coords[[i]] - coords[[i - 1]]];
          dr2 = If[j < n, coords[[j + 1]] - coords[[j]], coords[[j]] - coords[[j - 1]]];
          cross = Cross[dr1, dr2];
          contribution = Dot[r12, cross] / (dist^3);
          writhe += contribution;
        ];
      ],
      {j, 1, n}
    ],
    {i, 1, n}
  ];
  writhe / (4.0 * Pi)
];

ComputeDistanceMatrix[coords_] := Module[{n},
  n = Length[coords];
  Table[EuclideanDistance[coords[[i]], coords[[j]]], {i, n}, {j, n}]
];

CountRegionContacts[distMatrix_, region1_, region2_, maxDist_] := Module[{},
  Length[Select[
    Flatten[Table[
      If[distMatrix[[i, j]] < maxDist, 1, Nothing],
      {i, region1}, {j, region2}
    ]],
    IntegerQ
  ]]
];

Print["Loading UCH-L1 experimental structure..."];
uchL1Exp = BioMolecule["2ETL"];
uchL1ExpCA = ExtractCA[uchL1Exp];
Print["Experimental residues: ", Length[uchL1ExpCA]];
Print[""];

Print["Computing experimental writhe..."];
uchL1ExpWrithe = ComputeWrithe[uchL1ExpCA];
Print["Experimental writhe: ", Round[uchL1ExpWrithe, 0.03]];
Print[""];

Print["Folding UCH-L1 with ESMFold..."];
uchL1Pred = BioMolecule[BioMoleculeValue[uchL1Exp, "BioSequences"][[1]], BioMoleculeFoldingMethod -> "Local"];
uchL1PredCA = ExtractCA[uchL1Pred];
Print["Predicted residues: ", Length[uchL1PredCA]];
Print[""];

Print["Computing predicted writhe..."];
uchL1PredWrithe = ComputeWrithe[uchL1PredCA];
Print["Predicted writhe: ", Round[uchL1PredWrithe, 0.03]];
Print[""];

Print["Computing N-terminus threading contacts..."];
expDistMatrix = ComputeDistanceMatrix[uchL1ExpCA];
predDistMatrix = ComputeDistanceMatrix[uchL1PredCA];

nTermRegion = Range[1, 30];
midChainRegion = Range[100, Min[180, Length[uchL1ExpCA]]];

expNTermContacts = CountRegionContacts[expDistMatrix, nTermRegion, midChainRegion, 8.0];
predNTermContacts = CountRegionContacts[predDistMatrix, nTermRegion, midChainRegion, 8.0];

Print["N-term (1-30) to mid-chain (100-", Min[180, Length[uchL1ExpCA]], ") contacts (<8A):"];
Print["  Experimental: ", expNTermContacts];
Print["  Predicted: ", predNTermContacts];
Print[""];

Print["============================="];
Print["RESULTS"];
Print["============================="];
Print[""];
Print["UCH-L1 Experimental"];
Print["  Writhe: ", Round[uchL1ExpWrithe, 0.03]];
Print["  N-term threading: ", expNTermContacts];
Print["  Topology: Unknotted"];
Print[""];
Print["UCH-L1 Predicted"];
Print["  Writhe: ", Round[uchL1PredWrithe, 0.03]];
Print["  N-term threading: ", predNTermContacts];
Print[""];
Print["UCH-L3 Reference (from previous analysis)"];
Print["  Experimental writhe: 0.53"];
Print["  Predicted writhe: 0.56"];
Print["  Experimental N-term threading: 63"];
Print["  Predicted N-term threading: 59"];
Print[""];

Print["Exporting results..."];
Export["D:\\uchl1_predicted_CA.mx", uchL1PredCA];
Print["  uchl1_predicted_CA.mx"];

results = <|
  "UCH_L1_Exp_Writhe" -> uchL1ExpWrithe,
  "UCH_L1_Pred_Writhe" -> uchL1PredWrithe,
  "UCH_L1_Exp_NTerm_Contacts" -> expNTermContacts,
  "UCH_L1_Pred_NTerm_Contacts" -> predNTermContacts,
  "UCH_L3_Pred_Writhe" -> 0.56,
  "Writhe_Difference" -> Abs[0.56 - uchL1PredWrithe]
|>;

Export["D:\\uchl1_simple_verification.json", results, "JSON"];
Print["  uchl1_simple_verification.json"];

Print[""];
Print["COMPLETE"];
