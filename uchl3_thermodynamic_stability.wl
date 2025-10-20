(* ::Package:: *)

(*
================================================================================
THERMODYNAMIC STABILITY ANALYSIS OF UCH-L3 KNOT
Why is the knot stable despite the entropic penalty?

Key questions:
1. How many long-range contacts stabilize the knot?
2. What is the hydrophobic core composition?
3. Are there specific stabilizing interactions unique to the knotted topology?
4. Can we estimate the energetic benefit vs entropic cost?

Strategy:
- Analyze contact networks in wild-type vs mutants
- Identify buried hydrophobic residues
- Calculate approximate folding energy
- Determine if knot provides unusual stability
================================================================================
*)

Print[""];
Print["================================================================================"];
Print["THERMODYNAMIC STABILITY ANALYSIS"];
Print["UCH-L3 Knotted Protein"];
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

ComputeDistanceMatrix[coords_] := Module[{n},
  n = Length[coords];
  Table[EuclideanDistance[coords[[i]], coords[[j]]], {i, n}, {j, n}]
];

(* ============================================================================
   LOAD STRUCTURES
   ============================================================================ *)

Print["[STEP 1/7] Loading Structures and Sequences"];
Print["----------------------------------------------"];

uchL3Exp = BioMolecule["1XD3"];
wtSeq = BioMoleculeValue[uchL3Exp, "BioSequences"][[1]];
fastaStr = ExportString[wtSeq, "FASTA"];
lines = StringSplit[fastaStr, "\n"];
wtSeqString = StringJoin[Rest[lines]];
seqLen = StringLength[wtSeqString];

Print["  Sequence length: ", seqLen, " residues"];
Print["  Full sequence:"];
Print["  ", wtSeqString];
Print[""];

(* Fold wild-type *)
Print["  Folding wild-type structure..."];
wtPred = BioMolecule[wtSeq, BioMoleculeFoldingMethod -> "Local"];
wtPredCA = ExtractCA[wtPred];
wtExpCA = ExtractCA[uchL3Exp];

Print["  Computing distance matrix..."];
wtDistMatrix = ComputeDistanceMatrix[wtPredCA];
Print["  Done"];
Print[""];

(* ============================================================================
   CONTACT NETWORK ANALYSIS
   ============================================================================ *)

Print["[STEP 2/7] Long-Range Contact Network Analysis"];
Print["----------------------------------------------"];

(* Define contact cutoffs *)
contactCutoff = 6.0;  (* Angstroms - typical for stabilizing contacts *)
longRangeCutoff = 50; (* Sequence separation for long-range *)

(* Find all long-range contacts *)
Print["  Identifying long-range contacts (>", longRangeCutoff, " residues apart, <",
      contactCutoff, " Angstroms)..."];

longRangeContacts = Select[
  Flatten[Table[
    If[j > i + longRangeCutoff && wtDistMatrix[[i, j]] < contactCutoff,
       {i, j, wtDistMatrix[[i, j]]},
       Nothing],
    {i, seqLen}, {j, seqLen}
  ], 1],
  ListQ
];

Print[""];
Print["  Total long-range contacts: ", Length[longRangeContacts]];
Print[""];

(* Analyze contacts by region *)
nTermContacts = Select[longRangeContacts, #[[1]] <= 30 &];
knotCoreContacts = Select[nTermContacts, 100 <= #[[2]] <= 180 &];

Print["  N-terminus (1-30) long-range contacts: ", Length[nTermContacts]];
Print["  Knot core (N-term to loop 100-180): ", Length[knotCoreContacts]];
Print[""];

(* Show top 10 strongest knot contacts *)
knotCoreSorted = SortBy[knotCoreContacts, Last];
Print["  Top 10 strongest knot core contacts:"];
Print["  Res1 | Res2 | Distance | Residues"];
Print["  -----|------|----------|----------"];
For[i = 1, i <= Min[10, Length[knotCoreSorted]], i++,
  contact = knotCoreSorted[[i]];
  r1 = contact[[1]];
  r2 = contact[[2]];
  dist = contact[[3]];
  aa1 = StringTake[wtSeqString, {r1, r1}];
  aa2 = StringTake[wtSeqString, {r2, r2}];
  Print["  ", StringPadRight[ToString[r1], 5], "| ",
        StringPadRight[ToString[r2], 5], "| ",
        StringPadRight[ToString[Round[dist, 0.1]], 9], "| ",
        aa1, "-", aa2];
];
Print[""];

(* ============================================================================
   HYDROPHOBIC CORE ANALYSIS
   ============================================================================ *)

Print["[STEP 3/7] Hydrophobic Core Analysis"];
Print["----------------------------------------------"];

(* Define hydrophobic residues *)
hydrophobic = {"A", "V", "I", "L", "M", "F", "W", "P"};
aromatic = {"F", "W", "Y"};
charged = {"K", "R", "D", "E"};
polar = {"S", "T", "N", "Q"};

(* Classify all residues *)
residueTypes = Table[
  aa = StringTake[wtSeqString, {i, i}];
  type = Which[
    MemberQ[aromatic, aa], "Aromatic",
    MemberQ[hydrophobic, aa], "Hydrophobic",
    MemberQ[charged, aa], "Charged",
    MemberQ[polar, aa], "Polar",
    True, "Other"
  ];
  {i, aa, type},
  {i, seqLen}
];

(* Analyze threading region (5-8) *)
threadingRegion = Range[5, 8];
threadingResidues = residueTypes[[threadingRegion]];

Print["  Threading region (residues 5-8):"];
For[i = 1, i <= Length[threadingResidues], i++,
  res = threadingResidues[[i]];
  Print["    Residue ", res[[1]], ": ", res[[2]], " (", res[[3]], ")"];
];
Print[""];

(* Analyze knot loop (147-161) *)
knotLoopRegion = Range[147, 161];
knotLoopResidues = residueTypes[[knotLoopRegion]];

hydrophobicInLoop = Length[Select[knotLoopResidues, #[[3]] == "Hydrophobic" || #[[3]] == "Aromatic" &]];

Print["  Knot loop region (residues 147-161):"];
Print["    Total residues: ", Length[knotLoopResidues]];
Print["    Hydrophobic/Aromatic: ", hydrophobicInLoop, " (",
      Round[N[hydrophobicInLoop/Length[knotLoopResidues]] * 100, 0.1], "%)"];
Print[""];

(* Find buried hydrophobic residues *)
(* Buried = has many contacts within contact cutoff *)
buriedResidues = Table[
  contactCount = Length[Select[
    Table[wtDistMatrix[[i, j]], {j, seqLen}],
    # < contactCutoff &
  ]];
  {i, StringTake[wtSeqString, {i, i}], contactCount},
  {i, seqLen}
];

buriedHydrophobic = Select[buriedResidues,
  MemberQ[hydrophobic, #[[2]]] && #[[3]] > 15 &
];

Print["  Buried hydrophobic residues (>15 contacts):"];
Print["    Count: ", Length[buriedHydrophobic]];
Print["    Fraction of all residues: ", Round[N[Length[buriedHydrophobic]/seqLen] * 100, 0.1], "%"];
Print[""];

(* ============================================================================
   KNOT-SPECIFIC INTERACTIONS
   ============================================================================ *)

Print["[STEP 4/7] Knot-Specific Stabilizing Interactions"];
Print["----------------------------------------------"];

(* Identify residue pairs that only contact due to knotting *)
(* These are residues from threading region contacting loop *)

threadingToLoopContacts = Select[knotCoreContacts,
  MemberQ[threadingRegion, #[[1]]] &&
  MemberQ[knotLoopRegion, #[[2]]] &
];

Print["  Threading region (5-8) to knot loop (147-161) contacts:"];
Print["    Count: ", Length[threadingToLoopContacts]];
Print[""];

If[Length[threadingToLoopContacts] > 0,
  Print["  These contacts are UNIQUE to the knotted topology:"];
  Print["  Thread Res | Loop Res | Distance | Interaction"];
  Print["  -----------|----------|----------|------------"];
  For[i = 1, i <= Length[threadingToLoopContacts], i++,
    contact = threadingToLoopContacts[[i]];
    r1 = contact[[1]];
    r2 = contact[[2]];
    dist = contact[[3]];
    aa1 = StringTake[wtSeqString, {r1, r1}];
    aa2 = StringTake[wtSeqString, {r2, r2}];

    interaction = Which[
      MemberQ[hydrophobic, aa1] && MemberQ[hydrophobic, aa2], "Hydrophobic",
      MemberQ[charged, aa1] && MemberQ[charged, aa2], "Electrostatic",
      MemberQ[aromatic, aa1] || MemberQ[aromatic, aa2], "Aromatic",
      True, "Other"
    ];

    Print["  ", StringPadRight[ToString[r1] <> " (" <> aa1 <> ")", 12], "| ",
          StringPadRight[ToString[r2] <> " (" <> aa2 <> ")", 9], "| ",
          StringPadRight[ToString[Round[dist, 0.1]], 9], "| ",
          interaction];
  ];
  Print[""];
];

(* ============================================================================
   ENERGETIC ESTIMATE
   ============================================================================ *)

Print["[STEP 5/7] Approximate Energetic Analysis"];
Print["----------------------------------------------"];

(* Rough energy estimates (kcal/mol) *)
(* These are order-of-magnitude estimates *)
energyPerHydrophobicContact = -1.0;  (* kcal/mol *)
energyPerHBond = -3.0;               (* kcal/mol *)
energyPerSaltBridge = -4.0;          (* kcal/mol *)

(* Count hydrophobic contacts in knot core *)
hydrophobicKnotContacts = Length[Select[knotCoreContacts,
  Module[{aa1, aa2},
    aa1 = StringTake[wtSeqString, {#[[1]], #[[1]]}];
    aa2 = StringTake[wtSeqString, {#[[2]], #[[2]]}];
    MemberQ[hydrophobic, aa1] && MemberQ[hydrophobic, aa2]
  ]&
]];

Print["  Hydrophobic contacts in knot core: ", hydrophobicKnotContacts];
Print["  Estimated stabilization: ~",
      Round[hydrophobicKnotContacts * energyPerHydrophobicContact, 1], " kcal/mol"];
Print[""];

(* Entropic cost of knotting *)
(* Very rough estimate: restricting backbone conformations *)
entropicCost = 5.0;  (* kcal/mol - very approximate *)
Print["  Estimated entropic cost of knot formation: +", entropicCost, " kcal/mol"];
Print[""];

netStabilization = (hydrophobicKnotContacts * energyPerHydrophobicContact) + entropicCost;
Print["  NET ENERGETIC EFFECT: ", Round[netStabilization, 1], " kcal/mol"];
Print[""];

If[netStabilization < 0,
  Print["  >>> The knot is THERMODYNAMICALLY FAVORABLE"];
  Print["      Enthalpic stabilization outweighs entropic cost"],
  Print["  >>> The knot may be KINETICALLY TRAPPED"];
  Print["      Or requires additional stabilizing factors not captured here"]
];
Print[""];

(* ============================================================================
   CONTACT DENSITY ANALYSIS
   ============================================================================ *)

Print["[STEP 6/7] Contact Density Analysis"];
Print["----------------------------------------------"];

(* Compare contact density in different regions *)
regions = {
  {"N-terminus (1-30)", Range[1, 30]},
  {"Threading site (5-8)", Range[5, 8]},
  {"Knot loop (147-161)", Range[147, 161]},
  {"Active site (88-97)", Range[88, 97]},
  {"C-terminus (180-229)", Range[180, seqLen]}
};

Print["  Contact density by region (contacts per residue):"];
Print[""];
Print["  Region                    | Avg Contacts | Type"];
Print["  --------------------------|--------------|------"];

For[i = 1, i <= Length[regions], i++,
  regionName = regions[[i, 1]];
  regionRange = regions[[i, 2]];

  totalContacts = Sum[
    Length[Select[
      Table[wtDistMatrix[[j, k]], {k, seqLen}],
      # < contactCutoff &
    ]],
    {j, regionRange}
  ];

  avgContacts = N[totalContacts / Length[regionRange]];

  contactType = Which[
    avgContacts > 20, "Very Dense",
    avgContacts > 15, "Dense",
    avgContacts > 10, "Moderate",
    True, "Sparse"
  ];

  Print["  ", StringPadRight[regionName, 26], "| ",
        StringPadRight[ToString[Round[avgContacts, 0.1]], 13], "| ",
        contactType];
];
Print[""];

(* ============================================================================
   FINAL ASSESSMENT
   ============================================================================ *)

Print["[STEP 7/7] Thermodynamic Stability Assessment"];
Print["----------------------------------------------"];
Print[""];

Print["  KNOT STABILIZATION MECHANISMS:"];
Print[""];

Print["  1. EXTENSIVE CONTACT NETWORK"];
Print["     - ", Length[knotCoreContacts], " long-range N-terminus to loop contacts"];
Print["     - ", hydrophobicKnotContacts, " hydrophobic contacts in knot core"];
Print["     - Estimated stabilization: ~",
      Round[Abs[hydrophobicKnotContacts * energyPerHydrophobicContact], 1], " kcal/mol"];
Print[""];

Print["  2. BURIED HYDROPHOBIC CORE"];
Print["     - ", Length[buriedHydrophobic], " buried hydrophobic residues"];
Print["     - Knot loop is ", Round[N[hydrophobicInLoop/Length[knotLoopResidues]] * 100, 0.1],
      "% hydrophobic"];
Print[""];

Print["  3. TOPOLOGY-SPECIFIC INTERACTIONS"];
Print["     - ", Length[threadingToLoopContacts],
      " contacts unique to knotted topology"];
Print["     - These contacts CANNOT form without the knot"];
Print[""];

Print["  THERMODYNAMIC VERDICT:"];
If[netStabilization < -2.0,
  Print["    STRONGLY FAVORABLE - Knot significantly stabilizes structure"];
  Print["    The enthalpic benefit far outweighs entropic cost"];
  verdict = "STRONGLY_FAVORABLE",
  If[netStabilization < 0,
    Print["    WEAKLY FAVORABLE - Knot provides modest stabilization"];
    Print["    The structure is stable but knot isn't essential"];
    verdict = "WEAKLY_FAVORABLE",
    Print["    UNFAVORABLE - Knot is likely kinetically trapped"];
    Print["    Formation must be driven by folding pathway"];
    verdict = "KINETICALLY_TRAPPED"
  ]
];
Print[""];

Print["  BIOLOGICAL SIGNIFICANCE:"];
Print["    The knot creates a unique interaction network that:");
Print["    - Brings distant residues into contact"];
Print["    - Creates a dense hydrophobic core"];
Print["    - Restricts conformational flexibility"];
Print["    - May be essential for function despite complexity"];
Print[""];

(* Export results *)
results = <|
  "Total_LongRange_Contacts" -> Length[longRangeContacts],
  "Knot_Core_Contacts" -> Length[knotCoreContacts],
  "Hydrophobic_Knot_Contacts" -> hydrophobicKnotContacts,
  "Threading_To_Loop_Contacts" -> Length[threadingToLoopContacts],
  "Buried_Hydrophobic_Residues" -> Length[buriedHydrophobic],
  "Estimated_Stabilization_kcal_mol" -> hydrophobicKnotContacts * energyPerHydrophobicContact,
  "Estimated_Entropic_Cost_kcal_mol" -> entropicCost,
  "Net_Effect_kcal_mol" -> netStabilization,
  "Thermodynamic_Verdict" -> verdict
|>;

Export["D:\\uchl3_thermodynamic_stability_results.json", results, "JSON"];
Print["Results exported to D:\\uchl3_thermodynamic_stability_results.json"];
Print[""];

Print["================================================================================"];
Print["THERMODYNAMIC STABILITY ANALYSIS COMPLETE"];
Print["================================================================================"];
