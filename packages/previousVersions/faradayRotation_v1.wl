(* ::Package:: *)

dropBoxOn=FileNameJoin[{"Research","DataRuns"}];
If[SystemInformation["Kernel","OperatingSystem"]=="Unix",
researchFolder=FileNameJoin[{"/","home","karl","Research"}]; ,(* Home coputer is Linux, so needs "/" at beginning *)
researchFolder=FileNameJoin[{"C:","Users","kahrendsen2","Box Sync","Research"}];(* Lab computer is Windows, so needs "C:" at beginning"*)
];
SetDirectory[researchFolder];

SetDirectory[FileNameJoin[{"mathematicaRb","packages"}]];
Import["fileManipulation.wl"];
Import["physicalConstants.wl"];
ResetDirectory[];

Q=(re fge k \[Mu] )/h;
nDensC=(2  (\[Nu]0*1*^9)^2)/(c  Q);
nDensCLitaker=(2  \[Nu]0*1*^9)/(c Q);

CalculateNumberDensity[fitParameter_,integratedMagneticField_]:=nDensC*fitParameter/integratedMagneticField;

CalculateNumberDensityLitaker[fitParameter_,integratedMagneticField_]:=nDensCLitaker*fitParameter/integratedMagneticField;
(* The detuning should be input as GHznDe
The angle should be input as radians
*)
CalculateFittingParameters[faradayRotationOrderedPairs_]:=Module[{startPoints,nlmFit,fitParam},
Clear[a2];
startPoints={{a2,1},{a4,1}};
$Assumptions=a2>0;
$Assumptions=a4>0;
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{(a2 ((\[Delta]+\[Nu]0)*1*^9)^2)/(\[Delta]*1*^9)^2+(a4 ((\[Delta]+\[Nu]0)*1*^9)^2)/(\[Delta]*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]]
];

CalculateFittingParametersLitaker[faradayRotationOrderedPairs_]:=Module[{startPoints,nlmFit,fitParam},
Clear[a2];
startPoints={{a2,600},{a4,0}};
$Assumptions=a2>0;
$Assumptions=a4>0;
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{(a2(\[Delta]+\[Nu]0)*1*^9)/(\[Delta]*1*^9)^2+(a4(\[Delta]+\[Nu]0)*1*^9)/(\[Delta]*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]];
nlmFit
];

(* A file is split into two datasets, the header and the data. The header is information that shouldn't change during the run. *)
ConvertSectionedFaradayDataFileToDatasets[fileName_,wavemeterOffset_:1.0587]:=Module[{wavelengthColumnName,wavelengthLowestReasonable,import,headerInfo,datasetRaw,dataset},
wavelengthColumnName="WAV";


import=ImportFile[fileName];
headerInfo=import[[1]];
datasetRaw=import[[2]];
dataset=datasetRaw;

(* This if-statement evaluates only when there is a column named wavelength. *)
If[AnyTrue[Normal[Keys[dataset[[1]]]],StringContainsQ[wavelengthColumnName]],
	wavelengthLowestReasonable=793;
	dataset=Select[datasetRaw,#[wavelengthColumnName]>wavelengthLowestReasonable&];
	(* Creates a new column where the detuning is calculated *)
	dataset=dataset[All,<|#,"DETUNING"->c/100/(#[wavelengthColumnName])-\[Nu]0-wavemeterOffset|>&];
];

<|"header"->headerInfo,"dataset"->dataset|>
];

FourierAnalyzeSectionedFaradayDatasets[faradayDataset_,columnName_,sectionsColumnName_]:=Module[{},
];
