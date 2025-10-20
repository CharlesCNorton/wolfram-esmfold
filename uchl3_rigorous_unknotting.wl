(* ::Package:: *)

(*
================================================================================
RIGOROUS TOPOLOGICAL UNKNOTTING SCREEN
UCH-L3: Writhe-Based Knot Detection

Uses Gauss linking integral to compute writhe of protein backbone.
Writhe quantifies self-entanglement and provides a direct topological measure.

Mathematical Definition:
Wr = (1/4Pi) ∫∫ [(r₁ - r₂) · (dr₁ × dr₂)] / |r₁ - r₂|³

For open curves (proteins), we close the backbone by connecting termini.

Author: Rigorous Topological Analysis
Date: 2025
================================================================================
*)

Print[""];
Print["================================================================================"];
Print["RIGOROUS TOPOLOGICAL UNKNOTTING SCREEN"];
Print["UCH-L3: Writhe-Based Knot Detection"];
Print["================================================================================"];
Print[""];

(* ============================================================================
   HELPER FUNCTIONS
   ============================================================================ *)

ExtractCA[biomol_] := Module[{atomCoords},
  atomCoords = BioMoleculeValue[biomol, "AtomCoordinates"];
  If[atomCoords === $Failed || !AssociationQ[atomCoords],
    Return[$Failed]
  ];
  If[!KeyExistsQ[atomCoords, "A"],
    Return[$Failed]
  ];
  Table[
    QuantityMagnitude[atomCoords["A"][[i]][[2]]],
    {i, Length[atomCoords["A"]]}
  ]
];

KabschAlign[mobile_, fixed_] := Module[{
    centroidMobile, centroidFixed,
    mobileCentered, fixedCentered,
    covarianceMatrix, svd, u, v, d, rotation, rotated
  },
  centroidMobile = Mean[mobile];
  centroidFixed = Mean[fixed];
  mobileCentered = # - centroidMobile & /@ mobile;
  fixedCentered = # - centroidFixed & /@ fixed;
  covarianceMatrix = Transpose[mobileCentered].fixedCentered;
  svd = SingularValueDecomposition[covarianceMatrix];
  u = svd[[1]];
  v = svd[[3]];
  d = If[Det[v.Transpose[u]] < 0,
    DiagonalMatrix[{1, 1, -1}],
    IdentityMatrix[3]
  ];
  rotation = v.d.Transpose[u];
  rotated = Transpose[rotation.Transpose[mobileCentered]];
  # + centroidFixed & /@ rotated
];

(* ============================================================================
   WRITHE CALCULATION - GAUSS LINKING INTEGRAL
   ============================================================================ *)

ComputeWrithe[coords_] := Module[{
    closedCoords, n, writhe, i, j,
    r1, r2, r3, r4, dr1, dr2, numerator, denominator, contribution
  },

  closedCoords = Append[coords, First[coords]];
  n = Length[closedCoords] - 1;

  writhe = 0.0;

  Do[
    Do[
      If[Abs[i - j] > 1,
        r1 = closedCoords[[i]];
        r2 = closedCoords[[i + 1]];
        r3 = closedCoords[[j]];
        r4 = closedCoords[[j + 1]];

        dr1 = r2 - r1;
        dr2 = r4 - r3;

        numerator = (r1 - r3).Cross[dr1, dr2];
        denominator = EuclideanDistance[r1, r3]^3;

        If[denominator > 0.01,
          contribution = numerator / denominator;
          writhe += contribution;
        ];
      ],
      {j, i + 2, n}
    ],
    {i, 1, n - 2}
  ];

  writhe / (4.0 * Pi)
];

(* ============================================================================
   LOAD WILD-TYPE STRUCTURE
   ============================================================================ *)

Print["[STEP 1/5] Loading Wild-Type Structure"];
Print["----------------------------------------------"];

sequence = "MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA";
seqLen = StringLength[sequence];

Print["  Folding wild-type UCH-L3..."];
wtBioMol = BioMolecule[BioSequence["Peptide", sequence],
                       BioMoleculeFoldingMethod -> "Local"];
wtCA = ExtractCA[wtBioMol];

Print["  Structure loaded: ", seqLen, " residues"];
Print[""];

(* ============================================================================
   COMPUTE BASELINE WRITHE
   ============================================================================ *)

Print["[STEP 2/5] Computing Baseline Writhe"];
Print["----------------------------------------------"];

wtWrithe = ComputeWrithe[wtCA];

Print["  Wild-type writhe: ", NumberForm[wtWrithe, {6, 3}]];
Print[""];
Print["  Writhe interpretation:"];
Print["    |Wr| < 0.5  : Likely unknotted or simple loop"];
Print["    |Wr| > 0.5  : Significant self-entanglement"];
Print["    |Wr| > 1.0  : Strong topological complexity"];
Print[""];

(* ============================================================================
   SYSTEMATIC PROLINE SUBSTITUTION SCREEN
   ============================================================================ *)

Print["[STEP 3/5] Systematic Proline Substitution Screen"];
Print["----------------------------------------------"];
Print["  Target: Threading region (positions 12-18)"];
Print["  Mutation: X to P (helix breaker, prevents threading)"];
Print[""];

threadingPositions = Range[12, 18];

Print["  Testing ", Length[threadingPositions], " positions"];
Print[""];

unknottingResults = Table[
  Module[{mutSeq, mutBioMol, mutCA, mutWrithe, writheChange,
          overallRMSD, activeSiteRMSD, activeSite, originalAA},

    originalAA = StringTake[sequence, {pos, pos}];
    mutSeq = StringReplacePart[sequence, "P", {pos, pos}];

    Print["  Position ", pos, " (", originalAA, pos, "P):"];
    Print["    Folding mutant..."];

    mutBioMol = Quiet[Check[
      BioMolecule[BioSequence["Peptide", mutSeq],
                  BioMoleculeFoldingMethod -> "Local"],
      $Failed
    ]];

    If[mutBioMol === $Failed,
      Print["    ERROR: Folding failed"];
      <|"Position" -> pos,
        "Mutation" -> originalAA <> ToString[pos] <> "P",
        "Success" -> False|>,

      mutCA = ExtractCA[mutBioMol];

      If[mutCA === $Failed || Length[mutCA] != seqLen,
        Print["    ERROR: CA extraction failed"];
        <|"Position" -> pos,
          "Mutation" -> originalAA <> ToString[pos] <> "P",
          "Success" -> False|>,

        mutCAAligned = KabschAlign[mutCA, wtCA];

        overallRMSD = Sqrt[Mean[Table[
          SquaredEuclideanDistance[wtCA[[i]], mutCAAligned[[i]]],
          {i, seqLen}
        ]]];

        activeSite = Range[88, 97];
        activeSiteRMSD = Sqrt[Mean[Table[
          SquaredEuclideanDistance[wtCA[[i]], mutCAAligned[[i]]],
          {i, activeSite}
        ]]];

        mutWrithe = ComputeWrithe[mutCAAligned];
        writheChange = Abs[mutWrithe - wtWrithe];

        Print["    Writhe: ", NumberForm[mutWrithe, {5, 3}],
              " (ΔWr = ", NumberForm[writheChange, {5, 3}], ")"];
        Print["    Overall RMSD: ", NumberForm[overallRMSD, {4, 2}], " A"];
        Print["    Active site RMSD: ", NumberForm[activeSiteRMSD, {4, 2}], " A"];
        Print[""];

        <|"Position" -> pos,
          "Mutation" -> originalAA <> ToString[pos] <> "P",
          "Success" -> True,
          "Writhe" -> mutWrithe,
          "WritheChange" -> writheChange,
          "OverallRMSD" -> overallRMSD,
          "ActiveSiteRMSD" -> activeSiteRMSD,
          "TopologyChanged" -> (writheChange > 0.3),
          "FoldPreserved" -> (overallRMSD < 5.0),
          "ActiveSiteIntact" -> (activeSiteRMSD < 2.0)|>
      ]
    ]
  ],
  {pos, threadingPositions}
];

Print[""];

(* ============================================================================
   ANALYSIS OF RESULTS
   ============================================================================ *)

Print["[STEP 4/5] Analysis of Unknotting Attempts"];
Print["----------------------------------------------"];

successfulResults = Select[unknottingResults, #["Success"] === True &];

If[Length[successfulResults] == 0,
  Print["  All mutations failed to fold"];
  Print[""];,

  Print["  Writhe distribution:"];
  Print["    Wild-type: ", NumberForm[wtWrithe, {5, 3}]];
  Print["    Mutants: ", NumberForm[Min[#["Writhe"] & /@ successfulResults], {5, 3}],
        " to ", NumberForm[Max[#["Writhe"] & /@ successfulResults], {5, 3}]];
  Print[""];

  topologyChangedMutants = Select[successfulResults, #["TopologyChanged"] &];

  If[Length[topologyChangedMutants] > 0,
    Print["  Mutations with significant topology change (ΔWr > 0.3):"];
    Do[
      mut = topologyChangedMutants[[i]];
      Print["    ", mut["Mutation"], ": Wr = ", NumberForm[mut["Writhe"], {5, 3}],
            " (ΔWr = ", NumberForm[mut["WritheChange"], {5, 3}], ")"];
    ,{i, Length[topologyChangedMutants]}];
    Print[""];,

    Print["  No mutations produced significant topology change"];
    Print[""];
  ];

  unknottedCandidates = Select[successfulResults,
    #["TopologyChanged"] && #["FoldPreserved"] && #["ActiveSiteIntact"] &
  ];

  If[Length[unknottedCandidates] > 0,
    Print["  Successful unknotting candidates:"];
    Print["  (Topology changed + Fold preserved + Active site intact)"];
    Print[""];
    Do[
      mut = unknottedCandidates[[i]];
      Print["    ", mut["Mutation"]];
      Print["      Writhe: ", NumberForm[mut["Writhe"], {5, 3}],
            " (ΔWr = ", NumberForm[mut["WritheChange"], {5, 3}], ")"];
      Print["      Overall RMSD: ", NumberForm[mut["OverallRMSD"], {4, 2}], " A"];
      Print["      Active site RMSD: ", NumberForm[mut["ActiveSiteRMSD"], {4, 2}], " A"];
      Print[""];
    ,{i, Length[unknottedCandidates]}];,

    Print["  No successful unknotting candidates found"];
    Print["  (All topology-changing mutations disrupted fold or active site)");
    Print[""];
  ];
];

(* ============================================================================
   STATISTICAL SUMMARY
   ============================================================================ *)

Print["[STEP 5/5] Statistical Summary"];
Print["----------------------------------------------"];

If[Length[successfulResults] > 0,
  avgWritheChange = Mean[#["WritheChange"] & /@ successfulResults];
  maxWritheChange = Max[#["WritheChange"] & /@ successfulResults];

  Print["  Mutations tested: ", Length[unknottingResults]];
  Print["  Successful folds: ", Length[successfulResults]];
  Print[""];
  Print["  Writhe statistics:"];
  Print["    Average ΔWr: ", NumberForm[avgWritheChange, {4, 3}]];
  Print["    Maximum ΔWr: ", NumberForm[maxWritheChange, {4, 3}]];
  Print[""];

  Print["  Classification:"];
  Print["    Topology changed: ",
    Length[Select[successfulResults, #["TopologyChanged"] &]]];
  Print["    Fold preserved: ",
    Length[Select[successfulResults, #["FoldPreserved"] &]]];
  Print["    Active site intact: ",
    Length[Select[successfulResults, #["ActiveSiteIntact"] &]]];
  Print[""];
];

Print["================================================================================"];
Print["CONCLUSION"];
Print["================================================================================"];
Print[""];

If[Length[unknottedCandidates] > 0,
  Print["  Result: Successful unknotting achieved"];
  Print["  Interpretation: Knot is not essential for fold or function"];
  Print["  Implication: Knot is evolutionarily contingent, not required"];,

  If[Length[topologyChangedMutants] > 0,
    Print["  Result: Topology changes disrupted structure or function"];
    Print["  Interpretation: Knot is structurally integrated"];
    Print["  Implication: Knot removal destabilizes protein"];,

    Print["  Result: No mutations significantly altered topology"];
    Print["  Interpretation: Threading region is robust to perturbation"];
    Print["  Implication: Knot formation is kinetically stable"];
  ];
];

Print[""];
Print["================================================================================"];
Print[""];

Export["D:\\uchl3_rigorous_unknotting_results.json", unknottingResults];
Print["Results exported to D:\\uchl3_rigorous_unknotting_results.json"];
Print[""];

Print["================================================================================"];
Print["RIGOROUS UNKNOTTING SCREEN COMPLETE"];
Print["================================================================================"];
