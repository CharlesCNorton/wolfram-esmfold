Print[""];
Print["================================================================================"];
Print["TOPOLOGICAL WRITHE LOCALIZATION ANALYSIS"];
Print["UCH-L3 - Sliding Window Knot Core Identification"];
Print["================================================================================"];
Print[""];

Print["Objective:"];
Print["  Map writhe distribution along backbone using sliding windows"];
Print["  Identify minimal knotted region"];
Print["  Verify ESMFold places knot in correct sequence location"];
Print[""];

ComputeSegmentWrithe[coords_, startIdx_, endIdx_] := Module[{
    segCoords, n, writhe, i, j, r1, r2, r12, dr1, dr2, cross, dist, contribution
  },

  segCoords = coords[[startIdx ;; endIdx]];
  n = Length[segCoords];
  writhe = 0.0;

  Do[
    Do[
      If[Abs[i - j] > 2,
        r1 = segCoords[[i]];
        r2 = segCoords[[j]];
        r12 = r1 - r2;
        dist = Norm[r12];

        If[dist > 0.1,
          dr1 = If[i < n,
                   segCoords[[i + 1]] - segCoords[[i]],
                   segCoords[[i]] - segCoords[[i - 1]]];
          dr2 = If[j < n,
                   segCoords[[j + 1]] - segCoords[[j]],
                   segCoords[[j]] - segCoords[[j - 1]]];

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

Print["[STEP 1/5] Loading Structures"];
Print["----------------------------------------------"];
Print[""];

Print["  Loading experimental CA coordinates"];
expCA = Import["D:\\uchl3_experimental_CA.mx"];
nRes = Length[expCA];
Print["  Residues: ", nRes];
Print[""];

Print["  Loading predicted aligned CA coordinates"];
predCAAligned = Import["D:\\uchl3_predicted_aligned_CA.mx"];
Print["  Residues: ", Length[predCAAligned]];
Print[""];

Print["[STEP 2/5] Configuring Sliding Window Parameters"];
Print["----------------------------------------------"];
Print[""];

windowSize = 100;
stepSize = 10;
nWindows = Floor[(nRes - windowSize) / stepSize] + 1;

Print["  Window size: ", windowSize, " residues"];
Print["  Step size: ", stepSize, " residues"];
Print["  Number of windows: ", nWindows];
Print["  Coverage: residues 1-", nRes];
Print[""];

Print["  Rationale:"];
Print["    100-residue windows capture knot core structure"];
Print["    10-residue steps provide good spatial resolution"];
Print["    Expected knot core around residues 5-160"];
Print[""];

Print["[STEP 3/5] Computing Writhe Profile - Experimental Structure"];
Print["----------------------------------------------"];
Print[""];

expWritheProfile = {};
expWindowCenters = {};

Print["  Computing writhe for ", nWindows, " windows..."];
Print[""];

Do[
  startIdx = 1 + (i - 1) * stepSize;
  endIdx = Min[startIdx + windowSize - 1, nRes];
  centerRes = Round[(startIdx + endIdx) / 2];

  If[Mod[i, 5] == 0,
    Print["  Progress: ", Round[100.0 * i / nWindows, 0.1], "% (window ", i, "/", nWindows, ")"];
  ];

  wr = ComputeSegmentWrithe[expCA, startIdx, endIdx];
  AppendTo[expWritheProfile, wr];
  AppendTo[expWindowCenters, centerRes];
  ,
  {i, 1, nWindows}
];

Print[""];
Print["  Experimental profile computed"];
Print["  Mean writhe: ", Round[Mean[expWritheProfile], 0.01]];
Print["  Max writhe: ", Round[Max[expWritheProfile], 0.01],
      " at residue ", expWindowCenters[[Position[expWritheProfile, Max[expWritheProfile]][[1, 1]]]]];
Print["  Min writhe: ", Round[Min[expWritheProfile], 0.01],
      " at residue ", expWindowCenters[[Position[expWritheProfile, Min[expWritheProfile]][[1, 1]]]]];
Print[""];

Print["[STEP 4/5] Computing Writhe Profile - Predicted Structure"];
Print["----------------------------------------------"];
Print[""];

predWritheProfile = {};
predWindowCenters = {};

Print["  Computing writhe for ", nWindows, " windows..."];
Print[""];

Do[
  startIdx = 1 + (i - 1) * stepSize;
  endIdx = Min[startIdx + windowSize - 1, nRes];
  centerRes = Round[(startIdx + endIdx) / 2];

  If[Mod[i, 5] == 0,
    Print["  Progress: ", Round[100.0 * i / nWindows, 0.1], "% (window ", i, "/", nWindows, ")"];
  ];

  wr = ComputeSegmentWrithe[predCAAligned, startIdx, endIdx];
  AppendTo[predWritheProfile, wr];
  AppendTo[predWindowCenters, centerRes];
  ,
  {i, 1, nWindows}
];

Print[""];
Print["  Predicted profile computed"];
Print["  Mean writhe: ", Round[Mean[predWritheProfile], 0.01]];
Print["  Max writhe: ", Round[Max[predWritheProfile], 0.01],
      " at residue ", predWindowCenters[[Position[predWritheProfile, Max[predWritheProfile]][[1, 1]]]]];
Print["  Min writhe: ", Round[Min[predWritheProfile], 0.01],
      " at residue ", predWindowCenters[[Position[predWritheProfile, Min[predWritheProfile]][[1, 1]]]]];
Print[""];

Print["[STEP 5/5] Identifying Knot Core Regions"];
Print["----------------------------------------------"];
Print[""];

expMaxIdx = Position[expWritheProfile, Max[expWritheProfile]][[1, 1]];
predMaxIdx = Position[predWritheProfile, Max[predWritheProfile]][[1, 1]];

expKnotCenter = expWindowCenters[[expMaxIdx]];
predKnotCenter = predWindowCenters[[predMaxIdx]];

Print["  Experimental knot core center: residue ", expKnotCenter];
Print["  Predicted knot core center: residue ", predKnotCenter];
Print["  Localization error: ", Abs[expKnotCenter - predKnotCenter], " residues"];
Print[""];

threshold = 0.5 * Max[expWritheProfile];
expCoreWindows = Select[Range[Length[expWritheProfile]], expWritheProfile[[#]] > threshold &];
predCoreWindows = Select[Range[Length[predWritheProfile]], predWritheProfile[[#]] > threshold &];

expCoreStart = expWindowCenters[[Min[expCoreWindows]]];
expCoreEnd = expWindowCenters[[Max[expCoreWindows]]];
predCoreStart = predWindowCenters[[Min[predCoreWindows]]];
predCoreEnd = predWindowCenters[[Max[predCoreWindows]]];

Print["  Knot core region (>50% max writhe):"];
Print["    Experimental: residues ", expCoreStart, "-", expCoreEnd, " (span: ", expCoreEnd - expCoreStart, ")"];
Print["    Predicted: residues ", predCoreStart, "-", predCoreEnd, " (span: ", predCoreEnd - predCoreStart, ")"];
Print[""];

overlapStart = Max[expCoreStart, predCoreStart];
overlapEnd = Min[expCoreEnd, predCoreEnd];
If[overlapEnd >= overlapStart,
  overlap = overlapEnd - overlapStart;
  expSpan = expCoreEnd - expCoreStart;
  predSpan = predCoreEnd - predCoreStart;
  overlapFraction = overlap / Mean[{expSpan, predSpan}];
  Print["  Core region overlap: ", Round[100 * overlapFraction, 0.1], "%"];
  Print["  Overlap span: residues ", overlapStart, "-", overlapEnd, " (", overlap, " residues)"];
  ,
  Print["  Core region overlap: 0%"];
  Print["  No overlap detected"];
];
Print[""];

Print["================================================================================"];
Print["CORRELATION ANALYSIS"];
Print["================================================================================"];
Print[""];

profileCorrelation = Correlation[expWritheProfile, predWritheProfile];
Print["  Pearson correlation: ", Round[profileCorrelation, 0.001]];
Print[""];

If[profileCorrelation > 0.9,
  Print["  Interpretation: Excellent agreement"];
  Print["    ESMFold accurately reproduces writhe distribution"],
  If[profileCorrelation > 0.7,
    Print["  Interpretation: Good agreement"];
    Print["    ESMFold captures major topological features"],
    Print["  Interpretation: Moderate agreement"];
    Print["    Some topological localization differences present"]
  ]
];
Print[""];

profileRMSD = Sqrt[Mean[(expWritheProfile - predWritheProfile)^2]];
Print["  Profile RMSD: ", Round[profileRMSD, 0.001]];
Print["  Mean absolute difference: ", Round[Mean[Abs[expWritheProfile - predWritheProfile]], 0.001]];
Print[""];

Print["================================================================================"];
Print["BIOLOGICAL INTERPRETATION"];
Print["================================================================================"];
Print[""];

Print["Known UCH-L3 topology:"];
Print["  Threading region: residues ~5-30"];
Print["  Knot loop: residues ~100-160"];
Print["  Expected knot core center: ~80-100"];
Print[""];

Print["Observed knot core locations:"];
Print["  Experimental: residue ", expKnotCenter];
Print["  Predicted: residue ", predKnotCenter];
Print[""];

If[Abs[expKnotCenter - predKnotCenter] < 20,
  Print["  Assessment: Accurate knot localization"];
  Print["    ESMFold places knot in correct sequence position"];
  Print["    Localization error <20 residues confirms topology preserved"];
  ,
  Print["  Assessment: Knot localization deviation detected"];
  Print["    Localization error >20 residues warrants investigation"];
];
Print[""];

Print["Minimal knotted region:"];
Print["  Experimental core span: ", expCoreEnd - expCoreStart, " residues"];
Print["  Predicted core span: ", predCoreEnd - predCoreStart, " residues"];
Print[""];

If[Abs[(expCoreEnd - expCoreStart) - (predCoreEnd - predCoreStart)] < 30,
  Print["  Assessment: Core size well-preserved"];
  Print["    ESMFold accurately reproduces knot compactness"];
  ,
  Print["  Assessment: Core size difference detected"];
  Print["    Knot may be more/less compact in predicted structure"];
];
Print[""];

Print["================================================================================"];
Print["EXPORT"];
Print["================================================================================"];
Print[""];

localizationResults = <|
  "WindowSize" -> windowSize,
  "StepSize" -> stepSize,
  "NumberOfWindows" -> nWindows,
  "ExperimentalWritheProfile" -> expWritheProfile,
  "PredictedWritheProfile" -> predWritheProfile,
  "WindowCenters" -> expWindowCenters,
  "ExperimentalKnotCenter" -> expKnotCenter,
  "PredictedKnotCenter" -> predKnotCenter,
  "LocalizationError" -> Abs[expKnotCenter - predKnotCenter],
  "ExperimentalCoreRegion" -> {expCoreStart, expCoreEnd},
  "PredictedCoreRegion" -> {predCoreStart, predCoreEnd},
  "ProfileCorrelation" -> profileCorrelation,
  "ProfileRMSD" -> profileRMSD
|>;

Export["D:\\uchl3_writhe_localization.json", localizationResults, "JSON"];
Print["  uchl3_writhe_localization.json"];

localizationReport = StringJoin[
  "Topological Writhe Localization Analysis\n",
  "=========================================\n\n",
  "Protein: UCH-L3 (PDB: 1XD3)\n",
  "Analysis: Sliding window writhe profile\n",
  "Window size: ", ToString[windowSize], " residues\n",
  "Step size: ", ToString[stepSize], " residues\n\n",
  "Knot Core Localization\n",
  "----------------------\n",
  "Experimental center: residue ", ToString[expKnotCenter], "\n",
  "Predicted center: residue ", ToString[predKnotCenter], "\n",
  "Localization error: ", ToString[Abs[expKnotCenter - predKnotCenter]], " residues\n\n",
  "Core Region Boundaries\n",
  "----------------------\n",
  "Experimental: ", ToString[expCoreStart], "-", ToString[expCoreEnd],
    " (", ToString[expCoreEnd - expCoreStart], " residues)\n",
  "Predicted: ", ToString[predCoreStart], "-", ToString[predCoreEnd],
    " (", ToString[predCoreEnd - predCoreStart], " residues)\n\n",
  "Profile Statistics\n",
  "------------------\n",
  "Correlation: ", ToString[Round[profileCorrelation, 0.001]], "\n",
  "RMSD: ", ToString[Round[profileRMSD, 0.001]], "\n"
];

Export["D:\\uchl3_writhe_localization_report.txt", localizationReport, "Text"];
Print["  uchl3_writhe_localization_report.txt"];

profileData = Transpose[{expWindowCenters, expWritheProfile, predWritheProfile}];
Export["D:\\uchl3_writhe_profiles.csv", profileData, "CSV"];
Print["  uchl3_writhe_profiles.csv"];

Print[""];
Print["================================================================================"];
Print["ANALYSIS COMPLETE"];
Print["================================================================================"];
Print[""];

Print["Key findings:"];
Print["  1. Knot localization error: ", Abs[expKnotCenter - predKnotCenter], " residues"];
Print["  2. Profile correlation: ", Round[profileCorrelation, 0.001]];
Print["  3. Core region defined with high-writhe windows"];
Print["  4. Topological features spatially mapped along sequence"];
Print[""];
