Print[""];
Print["================================================================================"];
Print["TOPOLOGICAL WRITHE COMPUTATION"];
Print["UCH-L3 (PDB: 1XD3) - 5-crossing knotted protein"];
Print["================================================================================"];
Print[""];

ComputeWrithe[coords_] := Module[{
    n, writhe, i, j, r1, r2, r12, dr1, dr2, cross, dist, contribution
  },

  n = Length[coords];
  writhe = 0.0;

  Print["  Computing double integral over ", n, " backbone points"];
  Print["  Estimated time: ", Round[n^2/1000000], " seconds"];
  Print[""];

  Do[
    If[Mod[i, 10] == 0,
      Print["  Progress: ", Round[100.0 * i/n, 0.1], "%"];
    ];

    Do[
      If[Abs[i - j] > 2,
        r1 = coords[[i]];
        r2 = coords[[j]];
        r12 = r1 - r2;
        dist = Norm[r12];

        If[dist > 0.1,
          dr1 = If[i < n,
                   coords[[i + 1]] - coords[[i]],
                   coords[[i]] - coords[[i - 1]]];
          dr2 = If[j < n,
                   coords[[j + 1]] - coords[[j]],
                   coords[[j]] - coords[[j - 1]]];

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

Print["[STEP 1/3] Loading Structure Data"];
Print["----------------------------------------------"];
Print[""];

Print["  Loading experimental CA coordinates"];
expCA = Import["D:\\uchl3_experimental_CA.mx"];
Print["  Residues: ", Length[expCA]];
Print[""];

Print["  Loading predicted aligned CA coordinates"];
predCAAligned = Import["D:\\uchl3_predicted_aligned_CA.mx"];
Print["  Residues: ", Length[predCAAligned]];
Print[""];

If[Length[expCA] != Length[predCAAligned],
  Print["  ERROR: Length mismatch between structures"];
  Abort[];
];

Print["  Structures loaded"];
Print[""];

Print["[STEP 2/3] Computing Writhe - Experimental Structure"];
Print["----------------------------------------------"];
Print[""];

expWritheStart = AbsoluteTime[];
expWrithe = ComputeWrithe[expCA];
expWritheTime = AbsoluteTime[] - expWritheStart;

Print[""];
Print["  Experimental writhe: ", Round[expWrithe, 0.01]];
Print["  Computation time: ", Round[expWritheTime, 0.1], " seconds"];
Print[""];

Print["  Analysis:"];
Print["    Magnitude: ", Round[Abs[expWrithe], 0.01]];
Print["    Sign: ", If[expWrithe > 0, "positive (right-handed)", "negative (left-handed)"]];
Print["    Expected for 5_2 knot: |Wr| approx 0.5-2"];
Print[""];

Print["[STEP 3/3] Computing Writhe - Predicted Structure"];
Print["----------------------------------------------"];
Print[""];

predWritheStart = AbsoluteTime[];
predWrithe = ComputeWrithe[predCAAligned];
predWritheTime = AbsoluteTime[] - predWritheStart;

Print[""];
Print["  Predicted writhe: ", Round[predWrithe, 0.01]];
Print["  Computation time: ", Round[predWritheTime, 0.1], " seconds"];
Print[""];

Print["================================================================================"];
Print["COMPARISON"];
Print["================================================================================"];
Print[""];

Print["Experimental writhe:  ", Round[expWrithe, 0.01]];
Print["Predicted writhe:     ", Round[predWrithe, 0.01]];
Print[""];
Print["Absolute difference:  ", Round[Abs[expWrithe - predWrithe], 0.01]];
Print["Relative error:       ", Round[100 * Abs[expWrithe - predWrithe] / Abs[expWrithe], 0.1], "%"];
Print[""];

If[Sign[expWrithe] == Sign[predWrithe],
  Print["Chirality: preserved"];
  Print["  Both structures: ", If[expWrithe > 0, "right-handed", "left-handed"]];
  chiralityMatch = True,
  Print["Chirality: not preserved"];
  Print["  Experimental: ", If[expWrithe > 0, "right-handed", "left-handed"]];
  Print["  Predicted: ", If[predWrithe > 0, "right-handed", "left-handed"]];
  chiralityMatch = False
];
Print[""];

writheError = Abs[expWrithe - predWrithe];
relativeError = Abs[(expWrithe - predWrithe) / expWrithe];

Print["Assessment:"];
If[writheError < 0.1 && chiralityMatch,
  Print["  Topological invariant accurately reproduced"];
  Print["  Error within computational precision"];
  verdict = "accurate",
  If[writheError < 0.5 && chiralityMatch,
    Print["  Topological invariant approximately reproduced"];
    Print["  Measurable deviation present"];
    verdict = "approximate",
    Print["  Significant deviation in topological invariant"];
    If[!chiralityMatch, Print["  Chirality mismatch indicates topology error"]];
    verdict = "mismatch"
  ]
];
Print[""];

Print["================================================================================"];
Print["INTERPRETATION"];
Print["================================================================================"];
Print[""];

Print["Writhe definition:"];
Print["  Topological measure quantifying 3D curve self-linking"];
Print["  For knotted curves: Wr reflects crossing number and conformation"];
Print["  Sign indicates chirality (handedness)"];
Print[""];

Print["Theoretical expectation for 5_2 knot:"];
Print["  |Wr| typically 0.5-2 for protein conformations"];
Print["  (Lower than ideal mathematical knots due to loose structure)"];
Print["  Sign encodes chirality"];
Print[""];

Print["Observed results:"];
Print["  Experimental |Wr| = ", Round[Abs[expWrithe], 0.01]];
Print["  Predicted |Wr| = ", Round[Abs[predWrithe], 0.01]];
Print["  Difference: ", Round[writheError, 0.01],
      " (", Round[relativeError * 100, 0.1], "% relative error)"];
Print[""];

If[verdict == "accurate",
  Print["Conclusion:"];
  Print["  Predicted structure reproduces experimental topology"];
  Print["  Neural network prediction preserves topological invariant"];
  Print["  Writhe agreement confirms knot type preservation"];
];

If[verdict == "approximate",
  Print["Conclusion:"];
  Print["  Predicted structure approximately reproduces topology"];
  Print["  Some topological features preserved"];
  Print["  Quantitative deviation warrants further analysis"];
];

If[verdict == "mismatch",
  Print["Conclusion:"];
  Print["  Significant topological deviation detected"];
  Print["  Further investigation required"];
  If[!chiralityMatch,
    Print["  Chirality inversion suggests fundamental topology error"]];
];
Print[""];

Print["================================================================================"];
Print["EXPORT"];
Print["================================================================================"];
Print[""];

writheResults = <|
  "ExperimentalWrithe" -> expWrithe,
  "PredictedWrithe" -> predWrithe,
  "AbsoluteDifference" -> Abs[expWrithe - predWrithe],
  "RelativeError" -> relativeError,
  "ChiralityMatch" -> chiralityMatch,
  "Verdict" -> verdict
|>;

Export["D:\\uchl3_writhe_analysis.json", writheResults, "JSON"];
Print["  uchl3_writhe_analysis.json"];

writheReport = StringJoin[
  "Topological Writhe Analysis\n",
  "============================\n\n",
  "Protein: UCH-L3 (PDB: 1XD3)\n",
  "Knot type: 5-crossing (5_2)\n",
  "Residues: ", ToString[Length[expCA]], "\n\n",
  "Results\n",
  "-------\n",
  "Experimental writhe:  ", ToString[Round[expWrithe, 0.01]], "\n",
  "Predicted writhe:     ", ToString[Round[predWrithe, 0.01]], "\n",
  "Absolute difference:  ", ToString[Round[Abs[expWrithe - predWrithe], 0.01]], "\n",
  "Relative error:       ", ToString[Round[relativeError * 100, 0.1]], "%\n\n",
  "Chirality match: ", ToString[chiralityMatch], "\n",
  "Assessment: ", verdict, "\n"
];

Export["D:\\uchl3_writhe_report.txt", writheReport, "Text"];
Print["  uchl3_writhe_report.txt"];

Print[""];
Print["================================================================================"];
Print["ANALYSIS COMPLETE"];
Print["================================================================================"];
