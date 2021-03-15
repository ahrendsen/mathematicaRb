(* ::Package:: *)

(* ::Title:: *)
(*Faraday Rotation Functions*)


(* ::Chapter:: *)
(*Initialization*)


dropBoxOn=FileNameJoin[{"Research","DataRuns"}];
If[SystemInformation["Kernel","OperatingSystem"]=="Unix",
	researchFolder=FileNameJoin[{$HomeDirectory,"Research"}],(* Home coputer is Linux, so needs "/" at beginning *)
	researchFolder=FileNameJoin[{$HomeDirectory,"Box Sync","Research"}](* Lab computer is Windows, so needs "C:" at beginning"*)
];
SetDirectory[researchFolder];

SetDirectory[FileNameJoin[{"mathematicaRb","packages"}]];
Import["fileManipulation.wl"];
Import["dataManipulation.wl"];
Import["apparatus.wl"];
Import["physicalConstants.wl"];
ResetDirectory[];
ResetDirectory[];



(* ::Chapter::Closed:: *)
(*Number Density*)


(* ::Subchapter:: *)
(*Define Constants*)


Q=(re fge k \[Mu] )/h;
nDensC=(2  (\[Nu]0*1*^9)^2)/(c  Q);
nDensCLitaker=(2  \[Nu]0*1*^9)/(c Q);

	detCT="DETUNING";
	freqCT="FREQ";
	angCT="ANGLE";
	m1CT="magnet1Voltage(V)"; (* magnet 1 Column Title *)
	m2CT="magnet2Voltage(V)"; (* magnet 2 Column Title *)
	pdCT="PUMP"; (* Photodiode Column Title *)
	wavelengthCT="WAV";
	\[Delta]lowCut=-8;(* The lowest detuning that we want to keep from the dataset (default: 6, number density dependent) *)
	\[Delta]highCut=8;(* The lowest detuning that we want to keep from the dataset (default: 6, number density dependent) *)
	wavelengthLowestReasonable=793;
	wavemeterOffset=1.06; (* This value is subtracted from all wavelengths read in from data files.*) 


(* ::Subchapter:: *)
(*Equations to work with Data*)


CalculateNumberDensitySinglePoint[detuning_,integratedMagneticField_,angleShift_]:=nDensC/((\[Nu]0*1*^9)^2)* (detuning*1*^9)^2 * angleShift/integratedMagneticField;

CalculateNumberDensityMagneticFieldSlope[detuning_,slope_]:=nDensC/((\[Nu]0*1*^9)^2)* (detuning*1*^9)^2 * slope;

CalculateNumberDensity[fitParameter_,integratedMagneticField_]:=nDensC*fitParameter/integratedMagneticField;

CalculateNumberDensityLitaker[fitParameter_,integratedMagneticField_]:=nDensCLitaker*fitParameter/integratedMagneticField;
(* The detuning should be input as GHznDe
The angle should be input as radians
*)
CalculateFittingParameters[faradayRotationOrderedPairs_]:=Module[{startPoints,nlmFit,fitParam},
Clear[a2];
startPoints={{a2,1},{a4,1}};
(*
$Assumptions=a2>0;
$Assumptions=a4>0;
*)
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{(a2 ((\[Delta]+\[Nu]0)*1*^9)^2)/(\[Delta]*1*^9)^2+(a4 ((\[Delta]+\[Nu]0)*1*^9)^2)/(\[Delta]*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]]
];

CalculateFittingParametersLitaker[faradayRotationOrderedPairs_]:=Module[{startPoints,nlmFit,fitParam},
Clear[a2];
startPoints={{a2,600},{a4,0}};
(*
$Assumptions=a2>0;
$Assumptions=a4>0;
*)
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{(a2(\[Delta]+\[Nu]0)*1*^9)/(\[Delta]*1*^9)^2+(a4(\[Delta]+\[Nu]0)*1*^9)/(\[Delta]*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]];
nlmFit
];


(* ::Subchapter:: *)
(*Equations to Organize Data*)


(* A file is split into two datasets, the header and the data. The header is information that shouldn't change during the run. *)
ConvertSectionedFaradayDataFileToDatasets[fileName_,wavemeterOffset_:1.0587]:=Module[{wavelengthLowestReasonable,import,headerInfo,datasetRaw,dataset},
import=ImportFile[fileName];
headerInfo=import[[1]];
datasetRaw=import[[2]];
dataset=datasetRaw;

(* This if-statement evaluates only when there is a column named wavelength. *)
If[AnyTrue[Normal[Keys[dataset[[1]]]],StringContainsQ[wavelengthCT]],
	dataset=Select[dataset,#[wavelengthCT]>wavelengthLowestReasonable&];
	(* Creates a new column where the detuning is calculated *)
	dataset=dataset[All,<|#,detCT->c/100/(#[wavelengthCT])-\[Nu]0-wavemeterOffset|>&];
];

<|"header"->headerInfo,"dataset"->dataset|>
];

Clear[FourierAnalyzeSectionedFaradayDatasets];
FourierAnalyzeSectionedFaradayDatasets[faradayDataset_,
										columnName_,
										wavemeterOffset_:1.0586857110296357]:=
		Module[{lineData,allData,workingDataset,probeVoltages,
		listOrdPair,angles,i,discreteFourierData,transformedData,
		sn,cs,angle},
		
	workingDataset=faradayDataset;
	
	probeVoltages=Normal[Keys[workingDataset]];
	allData=<||>;
	lineData=<||>;

	For[i=1,i<=Length[probeVoltages],i++,
		discreteFourierData=Normal[workingDataset[i][All,columnName]];
		transformedData=DFT[discreteFourierData];
		AppendTo[lineData,"Fourier Data"->transformedData];
		AppendTo[lineData,detCT->probeVoltages[[i]]];
		AppendTo[allData,probeVoltages[[i]]->lineData];
	];
	allData
];

(* Returns the angle (in radians) of a given frequency component of 
a set of fourier components. The angle is determined by the arc tangent of the
sine and cosine components. *)
Clear[GetAnglesFromFourierDataset];
GetAnglesFromFourierDataset[fourierDataset_,frequencyComponent_]:=Module[{intensity,lineData,i,sectionHeaders,angle,sn,cs,allData},
	allData=<||>;
	sectionHeaders=Normal[Keys[fourierDataset]];
	For[i=1,i<=Length[sectionHeaders],i++,
		lineData=<||>;
		sn=fourierDataset[[i]]["Fourier Data","Sin Coefficients"];
		cs=fourierDataset[[i]]["Fourier Data","Cos Coefficients"];
		angle=-.5*ArcTan[sn[Key[frequencyComponent]]/cs[Key[frequencyComponent]]];
		AppendTo[lineData,"DETUNING"->sectionHeaders[[i]]];
		AppendTo[lineData,"ANGLE"->angle];
		AppendTo[allData,sectionHeaders[[i]]->lineData];
	];
	allData
];

(* Appends the angle (in radians) of a given frequency component of 
a set of fourier components. The angle is determined by the arc tangent of the
sine and cosine components. This function should replace GetAngles as I move forward with the code. *)
Clear[AppendAnglesToFourierDataset];
AppendAnglesToFourierDataset[fourierDataset_,frequencyComponent_]:=Module[{intensity,lineData,i,sectionHeaders,angle,sn,cs,allData},
	allData=<||>;
	sectionHeaders=Normal[Keys[fourierDataset]];
	For[i=1,i<=Length[sectionHeaders],i++,
		lineData=<||>;
		sn=fourierDataset[[i]]["Fourier Data","Sin Coefficients"];
		cs=fourierDataset[[i]]["Fourier Data","Cos Coefficients"];
		intensity=cs[Key[0]];
		angle=-.5*ArcTan[sn[Key[frequencyComponent]]/cs[Key[frequencyComponent]]];
		AppendTo[lineData,"DETUNING"->sectionHeaders[[i]]];
		AppendTo[lineData,"ANGLE"->angle];
		AppendTo[lineData,"INTENSITY"->intensity];
		AppendTo[allData,sectionHeaders[[i]]->lineData];
	];
	allData
];

ConvertDatasetWithRowHeadersToOrderedPairs[dataset_]:=Transpose[{Map[ToExpression,Keys[dataset]],Values[dataset]}];

SetMinAngleOfDatasetToZero[dataset_,detuningColumnName_:"DETUNING",angleColumnName_:"ANGLE"]:=Module[{min,orderedPairs},
min=MinMax[Transpose[orderedPairs][[2]]][[1]];
Transpose[{Transpose[orderedPairs][[1]],Transpose[orderedPairs][[2]]-min}]
];

ProcessFaradayRotationNumberDensityFile[fileName_]:=Module[{importedFile,header,m1,m2,
							magneticField,workingDataset,d,a,datasetNoAbs,minAngle,dataset,ordPair,nlmFit,a2err,n,nerr,return,lp,plot,viz},
	importedFile=ConvertSectionedFaradayDataFileToDatasets[fileName];
	header=importedFile[[1]];
	m1=ToExpression[header[m1CT]];
	m2=ToExpression[header[m2CT]];
	magneticField=CalculateIntegratedMagneticField[m1,m2];
	
	workingDataset=importedFile[[2]];
	d=Dataset[FourierAnalyzeSectionedFaradayDatasets[workingDataset,pdCT]];
	a=AppendAnglesToFourierDataset[d,4];
	datasetNoAbs=Select[Dataset[a],Abs[#[detCT]]>\[Delta]lowCut&];
	
	ProcessFaradayRotationNumberDensityDataset[datasetNoAbs,magneticField]
];

ProcessFaradayRotationBPDNumberDensityFile[fileName_]:=Module[{importedFile,header,m1,m2,
							magneticField,workingDataset,d,a,datasetNoAbs,minAngle,dataset,ordPair,nlmFit,a2err,n,nerr,return,lp,plot,viz},
	importedFile=ImportFile[fileName];
	header=importedFile[[1]];
	m1=ToExpression[header[m1CT]];
	m2=ToExpression[header[m2CT]];
	magneticField=CalculateIntegratedMagneticField[m1,m2];

	dataset=importedFile[[2]];
	(* This if-statement evaluates only when there is a column named wavelength. *)
	If[AnyTrue[Normal[Keys[dataset[[1]]]],StringContainsQ[wavelengthCT]],
		dataset=Select[dataset,#[wavelengthCT]>wavelengthLowestReasonable&];
		(* Creates a new column where the detuning is calculated *)
		dataset=dataset[All,<|#,detCT->c/100/(#[wavelengthCT])-\[Nu]0-wavemeterOffset|>&];
	];
	If[AnyTrue[Normal[Keys[dataset[[1]]]],StringContainsQ[freqCT]],
		dataset=Select[dataset,#[freqCT]>wavelengthLowestReasonable&];
		(* Creates a new column where the detuning is calculated *)
		dataset=dataset[All,<|#,detCT->#[freqCT]-\[Nu]0-wavemeterOffset,angCT->-#[angCT]|>&];
	];
	

	return=ProcessFaradayRotationNumberDensityDataset[dataset,magneticField];
	AppendTo[return,"Comments"->header["Comments"]]
];

ProcessFaradayRotationNumberDensityDataset[dataset_,magneticField_]:=Module[{importedFile,header,m1,m2,
							workingDataset,d,a,datasetNoAbs,minAngle,ordPair,nlmFit,a2err,n,nerr,return,lp,plot,viz},
							
	
	return=<||>;


	If[ToString[Head[Normal[dataset]]]=="List",
	ordPair=Normal[Values[dataset[All,{detCT,angCT}]]],
	ordPair=Normal[Values[dataset[Values,{detCT,angCT}]]]
	];
	AppendTo[return,"RawDetuning"->Transpose[ordPair][[1]]];
	AppendTo[return,"RawAngles"->Transpose[ordPair][[2]]];
	

	workingDataset=Select[dataset,#[detCT]<\[Delta]lowCut||#[detCT]>\[Delta]highCut&];
	minAngle=Normal[workingDataset[All,angCT]][[-1]];
	workingDataset=workingDataset[All,<|#,angCT->#[angCT]-minAngle|>&];
	AppendTo[return,"Working Dataset"->workingDataset];
	If[ToString[Head[Normal[workingDataset]]]=="List",
	ordPair=Normal[Values[workingDataset[All,{detCT,angCT}]]],
	ordPair=Normal[Values[workingDataset[Values,{detCT,angCT}]]]
	];
	AppendTo[return,"Detuning"->Transpose[ordPair][[1]]];
	AppendTo[return,"Angles"->Transpose[ordPair][[2]]];
	AppendTo[return,"Ordered Pairs"->ordPair];
	AppendTo[return,"Zero Angle"->minAngle];                                           


	nlmFit=CalculateFittingParameters[ordPair];
	AppendTo[return,"BestFitParameters"->nlmFit["BestFitParameters"]];
	a2=a2/.nlmFit["BestFitParameters"];
	a2err=nlmFit["ParameterErrors"][[1]];
	n=CalculateNumberDensity[a2,magneticField];
	nerr=CalculateNumberDensity[a2err,magneticField];
	AppendTo[return,<|"n_Rbx10^12"->n/1*^12,"n_Rberrx10^12"->nerr/1*^12|>];


	lp=ListPlot[Legended[ordPair,"Data"],PlotRange->Full,FrameLabel->{"Detuning (GHz)","\[CapitalDelta]\[Theta](rad)"}];
	plot=Plot[Legended[nlmFit[\[Delta]]/.nlmFit["BestFitParameters"],"Fit"],{\[Delta],-60,60}];
	viz=Show[{lp,plot}];

	AppendTo[return,"graph"->viz];
	AppendTo[return,"best fit function"->nlmFit["Function"][x]];
	return
];

(*
ListPlot[Transpose[{Map[ToExpression,Keys[a]],Values[a]}],PlotRange\[Rule]Full]
*)


(* ::Chapter::Closed:: *)
(*Polarization *)


(* ::Subchapter:: *)
(*Define Constants*)


(* I desperately want to treat the polarization formula in the same way that I have treated the number density formula. That is, to 
follow along with its derivation and account for missing factors. This will take a good deal of time though, and for now I will simply
use what I have been doing that has been giving me reasonable results. *)
(* Subscript[P, Rb]=polC (\[Delta](Subscript[\[Theta], p]-Subscript[\[Theta], B]))/(Subscript[n, Rb]*L)*)

(*\[CapitalDelta]\[Theta]= (L * Subscript[n, Rb] * Subscript[P, Rb])/polC 1/\[Delta]*)

(*slope=((L * Subscript[n, Rb]*Subscript[P, Rb])/polC) --> Subscript[P, Rb]=polCslope/(Subscript[n, Rb]*L)*)

(* THESE ARE ALMOST THE SAME VALUE*)
polC=(16 Pi gHzToHz)/(\[CapitalGamma] (\[Lambda]0*1*^-7)^2);(*MUNIR*)

polC=(4 gHzToHz)/(re*c*fge);(*KARL*)


(* ::Subchapter:: *)
(*Equations to work with Data*)


SinglePointPolarization[angleShift_,nrb_,detuning_]:=polC*detuning*angleShift/(nrb * cellLength);


TransformPolDataForPlotting[orderedPairs_]:=Module[{detunings,angleShifts},
	{detunings,angleShifts}=Transpose[orderedPairs];
	Transpose[{1/detunings,angleShifts}]
];

(* The ordered pairs should have detuning as the first item (GHz) and angle shift as the second item (rad) *)
Clear[CalculateRbPolarizationFit];
CalculateRbPolarizationFit[orderedPairs_]:=Module[{fittableData},
	fittableData=TransformPolDataForPlotting[orderedPairs];
	LinearModelFit[fittableData,x,x]
];

(* The function below should be used when the difference between S+ and S- is being used for the angle. *)
Clear[CalculateRubidiumPolarization];
CalculateRubidiumPolarizationTwoPumps[linearFitSlope_,numberDensity_]:=Module[{pRb},
	(* First we convert the angle shifts to inverse of that quantity so that we can do a linear fit *)
	pRb=polC (linearFitSlope/2)/(numberDensity * cellLength)
];

(* The Function below should be used when the difference between S+ (S-) and no pump is being used for the angle *)
Clear[CalculateRubidiumPolarization];
CalculateRubidiumPolarizationOnePump[linearFitSlope_,numberDensity_]:=Module[{pRb},
	(* First we convert the angle shifts to inverse of that quantity so that we can do a linear fit *)
	pRb=polC linearFitSlope/(numberDensity * cellLength)
];




(* ::Subchapter:: *)
(*Equations to Organize Data*)


CalculateAngleDifferenceOrdPairsFromFileNames[fileName1_,fileName2_,dataColumnName_]:=
	Module[
		{indexToRemove,differenceFirst,differenceLast,d,d1,d2,
		detuningLimit,det,a1,a2,frequencyComponent},
		
		detuningLimit=15;
		detuningLimitUpper=detuningLimit;
		detuningLimitLower=detuningLimit;
		frequencyComponent=4;
		det="DETUNING";
		d=ConvertSectionedFaradayDataFileToDatasets[fileName1][[2]];
		d=Dataset[FourierAnalyzeSectionedFaradayDatasets[d,dataColumnName]];
		d=Dataset[GetAnglesFromFourierDataset[d,frequencyComponent]];
		d=Normal[Values[Values[Select[d,Abs[#[det]]>detuningLimit&]]]];
		{d1,a1}=Transpose[d];

		d=ConvertSectionedFaradayDataFileToDatasets[fileName2][[2]];
		d=Dataset[FourierAnalyzeSectionedFaradayDatasets[d,dataColumnName]];
		d=Dataset[GetAnglesFromFourierDataset[d,frequencyComponent]];
		d=Normal[Values[Values[Select[d,Abs[#[det]]>detuningLimit&]]]];
		{d2,a2}=Transpose[d];

		If[Length[d1]!=  Length[d2],
			differenceFirst=d1[[1]]-d2[[1]];
			differenceLast=d1[[-1]]-d2[[-1]];
			If[Abs[differenceFirst]>Abs[differenceLast],
				indexToRemove=1;,
				indexToRemove=-1;
			];
			If[Length[d1]>Length[d2],
				a1=Drop[a1,indexToRemove];
				d1=Drop[d1,indexToRemove],
				a2=Drop[a2,indexToRemove];
				d2=Drop[d2,indexToRemove];
			];
		];
		Transpose[{d1,(a1-a2)}]
];


(* ::Subchapter:: *)
(*One File, One Function*)


(* The number density should be input in pcc*10^12 *)
Clear[ProcessFaradayRotationPolarizationFilesTwoPumps];
ProcessFaradayRotationPolarizationFilesTwoPumps[fileName1_,fileName2_,dataColumnName_,numberDensity_]:=Module[{ordPair,return,lm,pol,viz},
	return=<||>;
	ordPair=CalculateAngleDifferenceOrdPairsFromFileNames[fileName1,fileName2,dataColumnName];
	AppendTo[return,"Ordered Pairs"->ordPair];
	lm=CalculateRbPolarizationFit[ordPair];
	AppendTo[return,"Model Fit"->lm];
	pol=CalculateRubidiumPolarizationTwoPumps[lm["BestFitParameters"][[2]],numberDensity*1*^12];
	AppendTo[return,"Polarization"->pol];
	viz=Show[ListPlot[TransformPolDataForPlotting[ordPair],FrameLabel->{"1/\[Delta](\!\(\*SuperscriptBox[\(GHz\), \(-1\)]\))","\[CapitalDelta]\[Theta] (rad)"}],Plot[ lm["Function"][x],{x,-.1,.1}],FrameLabel->{"1/\[Delta](\!\(\*SuperscriptBox[\(GHz\), \(-1\)]\))","\[CapitalDelta]\[Theta] (rad)"},FrameLabel->{"1/\[Delta](\!\(\*SuperscriptBox[\(GHz\), \(-1\)]\))","\[CapitalDelta]\[Theta] (rad)"}];
	AppendTo[return,"Graph"->viz];
	return
];

(* The number density should be input in pcc*10^12 *)
Clear[ProcessFaradayRotationPolarizationFilesOnePump];
ProcessFaradayRotationPolarizationFilesOnePump[fileName1_,fileName2_,dataColumnName_,numberDensity_]:=
	Module[{ordPair,return,lm,pol,viz},
	
		return=<||>;
		ordPair=CalculateAngleDifferenceOrdPairsFromFileNames[fileName1,fileName2,dataColumnName];
		AppendTo[return,"Ordered Pairs"->ordPair];
		lm=CalculateRbPolarizationFit[ordPair];
		AppendTo[return,"Model Fit"->lm];
		pol=CalculateRubidiumPolarizationOnePump[lm["BestFitParameters"][[2]],numberDensity*1*^12];
		AppendTo[return,"Polarization"->pol];
		viz=Show[ListPlot[TransformPolDataForPlotting[ordPair]],Plot[lm["Function"][x],{x,-.1,.1}],FrameLabel->{"1/\[Delta](\!\(\*SuperscriptBox[\(GHz\), \(-1\)]\))","\[CapitalDelta]\[Theta](rad)"}];
		AppendTo[return,"Graph"->viz];
		AppendTo[return,"Dataset Notes"->"The ordered pairs are the detuning (GHz) and difference in angle (rad) between the two pump types."];
		return
];
