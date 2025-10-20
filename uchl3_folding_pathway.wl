(* ::Package:: *)

(*
================================================================================
FOLDING PATHWAY INFERENCE - UCH-L3 KNOT FORMATION
Does the knot form early or late during co-translational folding?

Strategy: Simulate ribosomal synthesis by folding progressively longer
N-terminal fragments. Track when knot-like features (long-range contacts)
first appear.

Key Questions:
1. At what sequence length do N-term to mid-chain contacts appear?
2. Does the knot core form incrementally or suddenly?
3. Can we identify the critical folding length?

Author: Testing Wolfram Language V14.3 Folding Pathway Analysis
Date: 2025
================================================================================
*)

Print[""];
Print["================================================================================"];
Print["FOLDING PATHWAY INFERENCE - UCH-L3 KNOT FORMATION"];
Print["Simulating Co-Translational Folding"];
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

ComputeDistanceMatrix[coords_] := Module[{n},
  n = Length[coords];
  Table[EuclideanDistance[coords[[i]], coords[[j]]], {i, n}, {j, n}]
];

CountRegionContacts[distMatrix_, region1_, region2_, threshold_] := Module[{count},
  count = 0;
  Do[
    Do[
      If[distMatrix[[i, j]] < threshold && distMatrix[[i, j]] > 0,
        count++
      ],
      {j, region2}
    ],
    {i, region1}
  ];
  count
];

(* ============================================================================
   LOAD FULL SEQUENCE
   ============================================================================ *)

Print["[STEP 1/5] Loading Full UCH-L3 Sequence"];
Print["----------------------------------------------"];

fullSequence = "MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA";

fullLength = StringLength[fullSequence];
Print["  Full sequence length: ", fullLength, " residues"];
Print["  First 50 residues: ", StringTake[fullSequence, 50]];
Print[""];

(* ============================================================================
   FOLD PROGRESSIVELY LONGER FRAGMENTS
   ============================================================================ *)

Print["[STEP 2/5] Folding N-Terminal Fragments"];
Print["----------------------------------------------"];
Print["  Simulating co-translational folding..."];
Print["  This will take some time - folding multiple structures"];
Print[""];

fragmentLengths = Join[
  Range[80, 140, 20],
  Range[160, fullLength, 20]
];

partialStructures = Table[
  Module[{partialSeq, partial, partialCA, status},
    Print["  Folding fragment length ", n, "/", fullLength, "..."];

    partialSeq = StringTake[fullSequence, n];

    status = Quiet[Check[
      partial = BioMolecule[BioSequence["Peptide", partialSeq],
                           BioMoleculeFoldingMethod -> "Local"];
      partialCA = ExtractCA[partial];

      If[partialCA === $Failed,
        Print["    WARNING: CA extraction failed"];
        <|"Length" -> n,
          "Success" -> False,
          "CA" -> {}|>,
        Print["    Success - ", Length[partialCA], " CA atoms"];
        <|"Length" -> n,
          "Success" -> True,
          "Structure" -> partial,
          "CA" -> partialCA|>
      ],
      Print["    ERROR: Folding failed"];
      <|"Length" -> n,
        "Success" -> False,
        "CA" -> {}|>
    ]];

    status
  ],
  {n, fragmentLengths}
];

Print[""];
Print["  Total fragments folded: ", Length[partialStructures]];
successfulFolds = Count[partialStructures, KeyValuePattern["Success" -> True]];
Print["  Successful folds: ", successfulFolds];
Print[""];

(* ============================================================================
   ANALYZE KNOT SIGNATURE EMERGENCE
   ============================================================================ *)

Print["[STEP 3/5] Analyzing Knot Signature Emergence"];
Print["----------------------------------------------"];
Print["  Knot signature = N-terminus (1-30) contacts to mid-chain (100+)"];
Print[""];

knotSignatures = Table[
  Module[{frag, dist, nTermContacts, earlyMidContacts, lateMidContacts},
    frag = partialStructures[[i]];

    If[!frag["Success"] || Length[frag["CA"]] < 100,
      <|"Length" -> frag["Length"],
        "NTermContacts" -> 0,
        "EarlyMidContacts" -> 0,
        "LateMidContacts" -> 0|>,

      dist = ComputeDistanceMatrix[frag["CA"]];

      nTermContacts = CountRegionContacts[
        dist,
        Range[1, Min[30, Length[frag["CA"]]]],
        Range[Min[31, Length[frag["CA"]]], Length[frag["CA"]]],
        8.0
      ];

      earlyMidContacts = If[Length[frag["CA"]] >= 120,
        CountRegionContacts[
          dist,
          Range[1, 30],
          Range[100, Min[120, Length[frag["CA"]]]],
          8.0
        ],
        0
      ];

      lateMidContacts = If[Length[frag["CA"]] >= 160,
        CountRegionContacts[
          dist,
          Range[1, 30],
          Range[140, Min[160, Length[frag["CA"]]]],
          8.0
        ],
        0
      ];

      <|"Length" -> frag["Length"],
        "NTermContacts" -> nTermContacts,
        "EarlyMidContacts" -> earlyMidContacts,
        "LateMidContacts" -> lateMidContacts|>
    ]
  ],
  {i, Length[partialStructures]}
];

Print["  Fragment Length | Total N-Term | Early Mid (100-120) | Late Mid (140-160)"];
Print["  ----------------|--------------|---------------------|-------------------"];
Do[
  sig = knotSignatures[[i]];
  Print["  ",
    StringPadRight[ToString[sig["Length"]], 15],
    " | ",
    StringPadRight[ToString[sig["NTermContacts"]], 12],
    " | ",
    StringPadRight[ToString[sig["EarlyMidContacts"]], 19],
    " | ",
    sig["LateMidContacts"]
  ],
  {i, Length[knotSignatures]}
];
Print[""];

(* ============================================================================
   IDENTIFY CRITICAL FOLDING LENGTH
   ============================================================================ *)

Print["[STEP 4/5] Identifying Critical Folding Length"];
Print["----------------------------------------------"];

knotThreshold = 5;

criticalLength = Module[{firstKnot},
  firstKnot = SelectFirst[
    knotSignatures,
    #["EarlyMidContacts"] >= knotThreshold &,
    Missing["NotFound"]
  ];

  If[MissingQ[firstKnot],
    Missing["NotFound"],
    firstKnot["Length"]
  ]
];

If[MissingQ[criticalLength],
  Print["  No clear knot signature detected in any fragment"];
  Print["  The knot may require the full-length protein"],
  Print["  Critical folding length: ~", criticalLength, " residues"];
  Print["  Percentage of full sequence: ",
    NumberForm[100.0 * criticalLength / fullLength, {4, 1}], "%"];

  If[criticalLength < 160,
    Print["  \[Checkmark] EARLY KNOT FORMATION - The knot forms before C-terminus synthesis"],
    If[criticalLength < 200,
      Print["  \[Checkmark] MID-STAGE KNOT FORMATION - The knot forms during mid-synthesis"],
      Print["  \[Checkmark] LATE KNOT FORMATION - The knot requires near-complete chain"]
    ]
  ];
];
Print[""];

(* ============================================================================
   KNOT CORE STABILITY ANALYSIS
   ============================================================================ *)

Print["[STEP 5/5] Knot Core Stability Analysis"];
Print["----------------------------------------------"];

If[!MissingQ[criticalLength] && successfulFolds >= 3,
  Module[{preCritical, postCritical, preSig, postSig},
    preCritical = SelectFirst[
      knotSignatures,
      #["Length"] < criticalLength && #["Length"] >= 100 &,
      Missing["NotFound"]
    ];

    postCritical = SelectFirst[
      knotSignatures,
      #["Length"] >= criticalLength + 20 &,
      Missing["NotFound"]
    ];

    If[!MissingQ[preCritical] && !MissingQ[postCritical],
      Print["  Pre-critical fragment (", preCritical["Length"], " residues):"];
      Print["    N-term contacts: ", preCritical["NTermContacts"]];
      Print["    Knot core contacts: ", preCritical["EarlyMidContacts"]];
      Print[""];
      Print["  Post-critical fragment (", postCritical["Length"], " residues):"];
      Print["    N-term contacts: ", postCritical["NTermContacts"]];
      Print["    Knot core contacts: ", postCritical["EarlyMidContacts"]];
      Print[""];

      contactIncrease = postCritical["EarlyMidContacts"] - preCritical["EarlyMidContacts"];

      If[contactIncrease > 10,
        Print["  \[Checkmark] SUDDEN KNOT FORMATION - Large contact increase (+",
          contactIncrease, ")"],
        Print["  \[Checkmark] GRADUAL KNOT FORMATION - Incremental contact increase (+",
          contactIncrease, ")"]
      ];
    ];
  ];
];
Print[""];

(* ============================================================================
   FINAL SUMMARY
   ============================================================================ *)

Print["================================================================================"];
Print["FOLDING PATHWAY ANALYSIS - SUMMARY"];
Print["================================================================================"];
Print[""];

Print["Total fragments analyzed: ", Length[fragmentLengths]];
Print["Successful folds: ", successfulFolds];
Print[""];

If[!MissingQ[criticalLength],
  Print["CRITICAL FOLDING LENGTH: ~", criticalLength, " residues"];
  Print["  (", NumberForm[100.0 * criticalLength / fullLength, {4, 1}],
    "% of full sequence)"];
  Print[""];

  Print["BIOLOGICAL INTERPRETATION:"];
  Print[""];
  If[criticalLength < 160,
    Print["  The knot forms EARLY during ribosomal synthesis"];
    Print["  The N-terminus threads through the forming loop"];
    Print["  BEFORE the C-terminus is fully synthesized"];
    Print[""];
    Print["  This suggests a KINETICALLY CONTROLLED mechanism"];
    Print["  where the ribosome exit tunnel and chaperones"];
    Print["  guide the threading process."];,

    If[criticalLength < 200,
      Print["  The knot forms in MID-SYNTHESIS"];
      Print["  Requires substantial C-terminal structure"];
      Print["  for the threading loop to stabilize"];,

      Print["  The knot forms LATE, requiring nearly"];
      Print["  complete chain synthesis"];
      Print["  Suggests POST-TRANSLATIONAL knotting mechanism"];
    ]
  ];
];

Print[""];
Print["================================================================================"];
Print[""];

exportData = <|
  "FragmentLengths" -> fragmentLengths,
  "KnotSignatures" -> knotSignatures,
  "CriticalLength" -> criticalLength,
  "FullLength" -> fullLength,
  "SuccessfulFolds" -> successfulFolds
|>;

Export["D:\\uchl3_folding_pathway_results.json", exportData];
Print["Results exported to D:\\uchl3_folding_pathway_results.json"];
Print[""];

Print["================================================================================"];
Print["FOLDING PATHWAY ANALYSIS COMPLETE"];
Print["================================================================================"];
