(* ::Package:: *)

(* ::Title:: *)
(*Faraday Rotation Functions*)


(* ::Chapter:: *)
(*Initialization*)


SetDirectory[$PACKAGES];
Import["fileManipulation.wl"];
Import["dataManipulation.wl"];
Import["apparatus.wl"];
Import["physicalConstants.wl"];
Import["generalPlottingFunctions.wl"];
ResetDirectory[];


(* ::Section:: *)
(*Define Constants*)


Q=(re fge k \[Mu] )/h;
nDensC=(2  (\[Nu]0*1*^9)^2)/(c  Q);
nDensCLitaker=(2  \[Nu]0*1*^9)/(c Q);
fileNameConst="File";
faradayRotationDensityModel[\[Delta]_,a2_]:=(a2 ((\[Delta]+\[Nu]0)*1*^9)^2)/(\[Delta]*1*^9)^2;
detCT="DETUNING";
freqCT="FREQ";
angCT="ANGLE";
m1CT="magnet1Voltage(V)"; (* magnet 1 Column Title *)
m2CT="magnet2Voltage(V)"; (* magnet 2 Column Title *)
pdCT="PUMP"; (* Photodiode Column Title *)
wavelengthCT="WAV";
\[Delta]Cut=4;
\[Delta]lowCut=\[Delta]Cut;(* The lowest detuning that we want to keep from the dataset (default: 6, number density dependent) *)
\[Delta]highCut=-\[Delta]Cut;(* The lowest detuning that we want to keep from the dataset (default: 6, number density dependent) *)
wavelengthLowestReasonable=793;
wavemeterOffset=1.06; (* This value is subtracted from all wavelengths read in from data files.*) 


(* ::Chapter:: *)
(*Single Faraday Rotation*)


(* ::Section:: *)
(*Read in Files*)


GetAngleFromFaradayRotationFile[fileName_]:=Module[
{f,d,h,fc,s,c,angle},
f=ImportFile[fileName];
h=f[[1]];
d=f[[2]];
fc=DFT[d[All,"HORIZ"]//Normal];
s=fc["Sin Coefficients"];
c=fc["Cos Coefficients"];
angle=-.5*ArcTan[s[4],c[4]]
];
SetAttributes[GetAngleFromFaradayRotationFile,Listable];


(* ::Chapter:: *)
(*Faraday Scans*)


(* ::Section:: *)
(*Density*)


(* ::Subsection::Closed:: *)
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
startPoints={{\[Theta]0,0},{a2,1},{a4,1}};
(*
$Assumptions=a2>0;
$Assumptions=a4>0;
*)
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{\[Theta]0+(*.012/62*\[Delta]-.006+*)(a2 ((\[Delta]+\[Nu]0)*1*^9)^2)/(\[Delta]*1*^9)^2+(a4 ((\[Delta]+\[Nu]0)*1*^9)^2)/((\[Delta])*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]]
];

CalculateFittingParametersFreeLineCenter[faradayRotationOrderedPairs_]:=Module[{startPoints,nlmFit,fitParam},
Clear[a2];
startPoints={{\[Theta]0,0},{d,1},{a2,1},{a4,1}};
(*
$Assumptions=a2>0;
$Assumptions=a4>0;
*)
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{\[Theta]0+(*.012/62*\[Delta]-.006+*)(a2 ((\[Delta]+d+\[Nu]0)*1*^9)^2)/((\[Delta]+d)*1*^9)^2+(a4 ((\[Delta]+d+\[Nu]0)*1*^9)^2)/((\[Delta]+d)*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]]
];

CalculateFittingParametersLitaker[faradayRotationOrderedPairs_]:=Module[{startPoints,nlmFit,fitParam},
Clear[a2];
startPoints={{\[Theta]0,0},{a2,600},{a4,0}};
(*
$Assumptions=a2>0;
$Assumptions=a4>0;
*)
nlmFit =NonlinearModelFit[Normal[faradayRotationOrderedPairs],{\[Theta]0+(a2(\[Delta]+\[Nu]0)*1*^9)/(\[Delta]*1*^9)^2+(a4(\[Delta]+\[Nu]0)*1*^9)/(\[Delta]*1*^9)^4,a2>0,a4>0},startPoints,\[Delta]];
nlmFit
];


(* ::Subsection:: *)
(*Equations to Organize Data*)


(* A file is split into two datasets, the header and the data. The header is information that shouldn't change during the run. *)
ConvertSectionedFaradayDataFileToDatasets[fileName_,wavemeterOffset_:1.0587]:=
Module[{wavelengthLowestReasonable,import,headerInfo,datasetRaw,dataset},
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

<|"header"->headerInfo,"dataset"->dataset,"Time"->ExtractTimeInfoFromFileNameString[headerInfo[fileNameConst]]|>
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
		
		(*Correct angles that are likely wrapping.*)
		(*angle=If[angle<-.1,angle=\[Pi]/4-angle,angle];*)
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
							magneticField,workingDataset,d,a,datasetNoAbs,minAngle,
							dataset,ordPair,nlmFit,a2err,n,nerr,return,lp,plot,viz,
							fullRotationDataset,allRawDataPlot,
							results},
	importedFile=ConvertSectionedFaradayDataFileToDatasets[fileName,wavemeterOffset];
	header=importedFile[[1]];
	m1=ToExpression[header[m1CT]];
	m2=ToExpression[header[m2CT]];
	(* Sometimes the magnets aren't recorded properly, in which case, hopefully the 
	approximate voltage is recorded some place and can be input below. *)
	(*m1=65.4;
	m2=64.9;*)
	magneticField=CalculateIntegratedMagneticField[m1,m2];
	
	workingDataset=importedFile[[2]];
	fullRotationDataset=workingDataset[Values,All,{"PUMP"}][Values]//Flatten;
	allRawDataPlot=ListPlot[fullRotationDataset,PlotRange->{{0,500},Automatic}];
	d=Dataset[FourierAnalyzeSectionedFaradayDatasets[workingDataset,pdCT]];
	a=AppendAnglesToFourierDataset[d,4];
	datasetNoAbs=Select[Dataset[a],Abs[#[detCT]]>\[Delta]lowCut&];
	
	results=ProcessFaradayRotationNumberDensityDataset[datasetNoAbs,magneticField];
	AppendTo[results,header];
	AppendTo[results,"Time"->ExtractTimeInfoFromFileNameString[header[fileNameConst]]];
	AppendTo[results,<|"allRawDataPlot"->allRawDataPlot,"allRawData"->Normal[fullRotationDataset]|>]
];
SetAttributes[ProcessFaradayRotationNumberDensityFile,Listable];

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
	AppendTo[return,<|"Comments"->header["Comments"],"Time"->ExtractTimeInfoFromFileNameString[header[fileNameConst]]|>]
];

ProcessFaradayRotationNumberDensityDataset[dataset_,magneticField_]:=Module[
{importedFile,header,m1,m2,
workingDataset,d,a,a2Val,a2PosInPar,
datasetNoAbs,minAngle,ordPair,nlmFit,nlmFit2,
a2err,n,nerr,return,lp,
plot,viz,graphNoLabel},
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
	AppendTo[return,"Working Dataset"->Normal[workingDataset]];
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
	a2Val=a2/.nlmFit["BestFitParameters"];
	a2PosInPar=Position[Keys[nlmFit["BestFitParameters"]],a2][[1,1]]; (* The position of a2 in the list of parameters *)
	a2err=nlmFit["ParameterErrors"][[a2PosInPar]];

	n=Around[CalculateNumberDensity[a2Val,magneticField],CalculateNumberDensity[a2err,magneticField]];
	AppendTo[return,<|"n_Rbx10^12"->n/(1*^12),"n_RbPercentErr"->n[[2]]/n[[1]]|>];


	lp=ListPlot[{Legended[ordPair,"Data"]},PlotRange->Full,FrameLabel->{"Detuning (GHz)","\[CapitalDelta]\[Theta](rad)"}];
	plot=Plot[Legended[nlmFit[\[Delta]]/.nlmFit["BestFitParameters"],"Fit"],
	{\[Delta]}\[Element]Interval[{-45,-10},{10,45}],
	FrameLabel->{"Detuning (GHz)","Rotation Angle (rad)"}
	];
	viz=Show[{plot,lp}];
	lp=ListPlot[{ordPair},PlotRange->Full];
	plot=Plot[nlmFit[\[Delta]]/.nlmFit["BestFitParameters"],{\[Delta],-60,60}];
	graphNoLabel=Show[{plot,lp}];

	AppendTo[return,<|"graph"->viz,"graphNoLabel"->graphNoLabel|>];
	AppendTo[return,"best fit function"->nlmFit["Function"][x]];
	return
];

(*
ListPlot[Transpose[{Map[ToExpression,Keys[a]],Values[a]}],PlotRange\[Rule]Full]
*)


(* ::Section:: *)
(*Polarization *)


(* ::Subsection::Closed:: *)
(*Define Constants*)


(* I desperately want to treat the polarization formula in the same way that I have treated the number density formula. That is, to 
follow along with its derivation and account for missing factors. This will take a good deal of time though, and for now I will simply
use what I have been doing that has been giving me reasonable results. *)
(* Subscript[P, Rb]=polC (\[Delta](Subscript[\[Theta], p]-Subscript[\[Theta], B]))/(Subscript[n, Rb]*L)*)

(*\[CapitalDelta]\[Theta]= (L * Subscript[n, Rb] * Subscript[P, Rb])/polC 1/\[Delta]*)

(*slope=((L * Subscript[n, Rb]*Subscript[P, Rb])/polC) --> Subscript[P, Rb]=polCslope/(Subscript[n, Rb]*L)*)

(* THESE ARE ALMOST THE SAME VALUE*)
(*2019-08-08: Added the negative sign so that the sign of the electron polarization and Rb polarization agree *)
polC=-(16 Pi gHzToHz)/(\[CapitalGamma] (\[Lambda]0*1*^-7)^2);(*MUNIR*)

polC=-(4 gHzToHz)/(re*c*fge);(*KARL*)

lowerLimit=15;
upperLimit=100;
(* Detuning Limits *)
dlmaxp=upperLimit; (* Detuning Limit Maximum for positive values *)
dlminp=lowerLimit; (* Detuning limit Minimum for positive values *)
dlmaxn=-upperLimit; (* Detuning limit Maximum for negative values *)
dlminn=-lowerLimit; (* Detuning limit Minimum for negative values *)


(* ::Subsection:: *)
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
	(* New method *)
	GayLinearFitNoError[fittableData]
	(* Old method *)
	(*** LinearModelFit[fittableData,x,x] ***)
	(*** Adapted for new format ***)
	(* 
	lmf=LinearModelFit[fittableData,x,x]
	Around[lmf["BestFitParameters"][[2]],lmf["ParameterErrors"][[2]]]
	*)
];

(* The function below should be used when the difference between S+ and S- is being used for the angle. *)
Clear[CalculateRubidiumPolarizationTwoPumps];
CalculateRubidiumPolarizationTwoPumps[linearFitSlope_,numberDensity_]:=Module[{pRb},
	(* First we convert the angle shifts to inverse of that quantity so that we can do a linear fit *)
	pRb=polC (linearFitSlope/2)/(numberDensity * cellLength)
];

(* The Function below should be used when the difference between S+ (S-) and no pump is being used for the angle *)
Clear[CalculateRubidiumPolarizationOnePump];
CalculateRubidiumPolarizationOnePump[linearFitSlope_,numberDensity_]:=Module[{pRb},
	(* First we convert the angle shifts to inverse of that quantity so that we can do a linear fit *)
	pRb=polC linearFitSlope/(numberDensity * cellLength)
];


FixWrapping[datainput_,signinput_]:=Module[
{data=datainput,data2,t,fix,div,lm,ignoreme,x,i,sign},
Which[signinput=="+",sign=1,signinput=="-",sign=-1];
FixWrapMakeMonotonic[data,sign]
];
(* This function was developed to account for wrapping in faraday rotation data. *)
FixWrapMakeMonotonic[ordPair_,posOrNegSlope_]:=Module[
{grouped,posDetuningSorted,negDetuningSorted,
t,data,i,sign,mod,fixedData},
t=Transpose[ordPair];
(* Put the data in the form 1/detuning vs. rotation angle, so the variables should be linearly related *)
data=Transpose[{1/t[[1]],t[[2]]}];
(* 
 * We know that when pumping with S+ light, there should be a negative slope, 
 * negative detunings should have a positive rotation and positive detunings 
* should have a negative rotation. The opposite is true for S- light. The
* next bit of code adds wrapping amounts until this condition is satisfied.
*)
For[i=1,i<=Length[data],i++,
(* If the data point is not in the right quadrant,
within some tolerance. *)
If[data[[i]][[1]]*data[[i]][[2]]*posOrNegSlope<0,
(* For negative slopes (positive slopes),
 * If the detuning is negative (positive),
 * numbers must be added (subtracted) 
* to move the data point to the correct quadrant. 
* The detuning being the opposite means
* the opposite operation is needed
 *)
If[data[[i]][[1]]*posOrNegSlope>0,mod=1,mod=-1];
(* Make the corresponding modification *)
data[[i]][[2]]=data[[i]][[2]]+mod*\[Pi]/2;
];
];
(*
 * But it's still possible more modification is needed to the dataset.
* The next piece of code ensures that the rotation amount increases 
* as we move to smaller detunings by adding wrapping amounts until 
* smaller detunings have larger rotation angles. 
*)
(* I start by separating into positive and negative detunings *)
grouped=GroupBy[data,Positive[#1[[1]]]&];
(* Then sort them so that we start at large detunings and work
* towards smaller ones. *)
posDetuningSorted=Sort[grouped[True],Abs[#1[[1]]]<Abs[#2[[1]]]&];
negDetuningSorted=Sort[grouped[False],Abs[#1[[1]]]<Abs[#2[[1]]]&];
(* This function adds the wrapping. The second argument should
* be 1 if rotation is more positive as detuning decreases. It should
* be -1 if rotation is more negative as detuning decreases.
*)
posDetuningSorted=AddWrap[posDetuningSorted,posOrNegSlope];
negDetuningSorted=AddWrap[negDetuningSorted,posOrNegSlope*-1];
(* Package everything together and ship it out in the original format *)
fixedData=Join[posDetuningSorted,negDetuningSorted];
t=Transpose[fixedData];
fixedData=Transpose[{1/t[[1]],t[[2]]}];
fixedData=Sort[fixedData,#1[[1]]<#2[[1]]&]
];

(* This is a helper function for Fixing the rotation. It makes sure that the
* amount of rotation always increases with decreasing detuning.
*)
Clear[AddWrap];
AddWrap[sortedDataset_List,modifier_]:=Module[
{i,next,tol,workingDataset=sortedDataset},
For[i=1,i<Length[workingDataset],i++,
next=i+1;
(* Check that the smaller detuning is a larger value, within some tolerance.*)
tol=1.1;(*tolerance is 50 percent *)
While[Abs[workingDataset[[i,2]]]>Abs[workingDataset[[next,2]]]*tol,
workingDataset[[next,2]]+=(modifier*\[Pi]/2);
];
];
workingDataset
];



(* ::Subsection:: *)
(*Equations to Organize Data*)


CalculateAngleDifferenceOrdPairsFromFileNames[fileNameFinal_,fileNameInitial_,dataColumnName_]:=
	Module[
		{indexToRemove,differenceFirst,differenceLast,d,dd,d1,d2,
		detuningLimit,det,a1,a2,frequencyComponent,return},
		return=<||>;
		
		detuningLimit=\[Delta]Cut;
		frequencyComponent=4;
		det="DETUNING";
		d=ConvertSectionedFaradayDataFileToDatasets[fileNameFinal][[2]];
		d=Dataset[FourierAnalyzeSectionedFaradayDatasets[d,dataColumnName]];
		d=Dataset[GetAnglesFromFourierDataset[d,frequencyComponent]];
		d=Normal[Values[Values[Select[d,Abs[#[det]]>detuningLimit&]]]];
		{d1,a1}=Transpose[d];

		dd=ConvertSectionedFaradayDataFileToDatasets[fileNameInitial][[2]];
		dd=Dataset[FourierAnalyzeSectionedFaradayDatasets[dd,dataColumnName]];
		dd=Dataset[GetAnglesFromFourierDataset[dd,frequencyComponent]];
		dd=Normal[Values[Values[Select[dd,Abs[#[det]]>detuningLimit&]]]];
		{d2,a2}=Transpose[dd];
		
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
		AppendTo[return,"orderedPairs"->Transpose[{d1,(a1-a2)}]];
		AppendTo[return,"rawAnglesGraph"->ListPlot[{d,dd},PlotRange->Full]];
		AppendTo[return,<|"finalAngles"->d,"initialAngles"->dd|>];
		return
		
];

ExportRbPolOrderedPairsFromDataset[dataset_]:=Module[
{export,names,values,fileNamesAndData},
export=dataset[All,All,"S+ Pump","orderedPairs"]//Normal;
names=Flatten[Keys[export]];
values=Flatten[Values[export][[All,1]],0];
fileNamesAndData=Partition[Riffle[names,values],2];
Apply[Export[FileBaseName[#1]<>"_rbPol.csv",#2]&,fileNamesAndData,1]
];


(* ::Subsection:: *)
(*One File, One Function*)


(* The number density should be input in pcc*10^12 *)
Clear[ProcessFaradayRotationPolarizationFilesThreePumps];
ProcessFaradayRotationPolarizationFilesThreePumps[unPumped_,
												sPlusPump_,
												sMinusPump_,
												dataColumnName_,
												numberDensity_Around]:=
Module[{fileResults,return,lm,
		pol,viz,results,polarizationResults,
		detuning,angles,orderedPairs},
	return=<||>;
	polarizationResults=<||>;
	fileResults=CalculateAngleDifferenceOrdPairsFromFileNames[sPlusPump,sMinusPump,dataColumnName];
	AppendTo[polarizationResults,fileResults];
	{detuning,angles}=Transpose[fileResults["orderedPairs"]];
	orderedPairs=Transpose[{detuning,angles/2}];
	orderedPairs=FixWrapping[orderedPairs,"-"];
	AppendTo[polarizationResults,"n_Rb"->numberDensity];
	AppendTo[polarizationResults,"orderedPairsCorrected"->orderedPairs];
	results=CalculateRubidiumPolarizationFromOrderedPairs[orderedPairs,numberDensity];
	AppendTo[polarizationResults,results];
	AppendTo[polarizationResults,"Time"->GetTimeInfoFromFileNameString[unPumped]];
	AppendTo[return,"Both Pumps"->polarizationResults];
	
	fileResults=CalculateAngleDifferenceOrdPairsFromFileNames[sPlusPump,unPumped,dataColumnName];
	AppendTo[polarizationResults,fileResults];
	{detuning,angles}=Transpose[fileResults["orderedPairs"]];
	orderedPairs=Transpose[{detuning,angles}];
	orderedPairs=FixWrapping[orderedPairs,"-"];
	AppendTo[polarizationResults,"n_Rb"->numberDensity];
	AppendTo[polarizationResults,"orderedPairsCorrected"->orderedPairs];
	results=CalculateRubidiumPolarizationFromOrderedPairs[orderedPairs,numberDensity];
	AppendTo[polarizationResults,results];
	AppendTo[polarizationResults,"Time"->GetTimeInfoFromFileNameString[unPumped]];
	AppendTo[polarizationResults,"File"->sPlusPump];
	AppendTo[return,"S+ Pump"->polarizationResults];
	
	fileResults=CalculateAngleDifferenceOrdPairsFromFileNames[sMinusPump,unPumped,dataColumnName];
	AppendTo[polarizationResults,fileResults];
	{detuning,angles}=Transpose[fileResults["orderedPairs"]];
	orderedPairs=Transpose[{detuning,angles}];
	orderedPairs=FixWrapping[orderedPairs,"+"];
	AppendTo[polarizationResults,"n_Rb"->numberDensity];
	AppendTo[polarizationResults,"orderedPairsCorrected"->orderedPairs];
	results=CalculateRubidiumPolarizationFromOrderedPairs[orderedPairs,numberDensity];
	AppendTo[polarizationResults,results];
	AppendTo[polarizationResults,"Time"->GetTimeInfoFromFileNameString[unPumped]];
	AppendTo[polarizationResults,"File"->sMinusPump];
	AppendTo[return,"S- Pump"->polarizationResults];
	
	return
];

(* The number density should be input in pcc*10^12 *)
Clear[ProcessFaradayRotationPolarizationFilesTwoPumps];
ProcessFaradayRotationPolarizationFilesTwoPumps[fileName1_,
												fileName2_,
												dataColumnName_,
												numberDensity_Around]:=
Module[{fileResults,return,lm,
		pol,viz,results,
		detuning,angles,orderedPairs},
	return=<||>;
	fileResults=CalculateAngleDifferenceOrdPairsFromFileNames[fileName1,fileName2,dataColumnName];
	AppendTo[return,fileResults];
	{detuning,angles}=Transpose[fileResults["orderedPairs"]];
	orderedPairs=Transpose[{detuning,angles/2}];
	
	AppendTo[return,"n_Rb"->numberDensity];
	
	results=CalculateRubidiumPolarizationFromOrderedPairs[orderedPairs,numberDensity];
	AppendTo[return,results];
	return
];

(* The number density should be input in pcc*10^12 *)
Clear[ProcessRbPolarizationTwoPumpsDataRuns];
ProcessRbPolarizationTwoPumpsDataRuns[fileNames_,
									dataColumnName_,
									numberDensity_Around]:=
Module[{fileResults,return,lm,
		pol,viz,results,results2,results3,
		detuning,angles,orderedPairs,
		paired,t,pumpedPlus,pumpedMinus,
		unpumped,n},
		If[TrueQ[Head[numberDensity]==List],n=numberDensity,n=ConstantArray[numberDensity,Length[fileNames]/3]];
	(* Make a new list from the long list of file names in which each item is a pair of a pumped and unpumped file.*)
	paired=Partition[fileNames,3];
	(* Restructure the file names so that we have one list of pump files, and a second list of unpumped files *)
	t=Transpose[paired];
	unpumped=t[[1]];
	pumpedPlus=t[[2]];
	pumpedMinus=t[[3]];
	results=MapThread[ProcessFaradayRotationPolarizationFilesTwoPumps,{pumpedPlus,pumpedMinus,ConstantArray[dataColumnName,Length[pumpedPlus]],ConstantArray[numberDensity,Length[pumpedPlus]]}];
	Join[results,results2,results3]
	AssociationThread[unpumped->results]
];

(* The number density should be input in pcc*10^12 *)
Clear[ProcessFaradayRotationPolarizationFilesOnePump];
ProcessFaradayRotationPolarizationFilesOnePump[fileName1_,fileName2_,dataColumnName_,numberDensity_]:=
	Module[{ordPair,return,lm,pol,viz,results,twoFileResults},
	
		return=<||>;
		twoFileResults=CalculateAngleDifferenceOrdPairsFromFileNames[fileName1,fileName2,dataColumnName];
		AppendTo[return,twoFileResults];
		AppendTo[return,"n_Rb"->numberDensity];
		results=CalculateRubidiumPolarizationFromOrderedPairs[twoFileResults["orderedPairs"],numberDensity];
		AppendTo[return,results]
];

ProcessRbPolarizationOnePumpDataRuns[fileNames_,signalColumnTitle_,numberDensity_]:=Module[
{t,paired,pumped,unpumped},
	(* Make a new list from the long list of file names in which each item is a pair of a pumped and unpumped file.*)
	paired=Partition[fileNames,2];
	(* Restructure the file names so that we have one list of pump files, and a second list of unpumped files *)
	t=Transpose[paired];
	pumped=t[[1]];
	unpumped=t[[2]];
	MapThread[ProcessFaradayRotationPolarizationFilesOnePump,{pumped,unpumped,ConstantArray[signalColumnTitle,Length[pumped]],ConstantArray[numberDensity,Length[pumped]]}]
];


CalculateRubidiumPolarizationFromOrderedPairs[ordPair_,numDensity_Around]:=
	Module[{return,lm,pol,viz,polErr},
		return=<||>;
		lm=CalculateRbPolarizationFit[ordPair];
		AppendTo[return,"Model Fit"->lm];
		pol=CalculateRubidiumPolarizationOnePump[lm,numDensity*1*^12];
		(*2019-09-25: pol=CalculateRubidiumPolarizationOnePump[lm["BestFitParameters"][[2]],numDensity*1*^12];*)
		polErr=pol[[2]]/pol[[1]];
		AppendTo[return,"P_Rb"->pol];
		AppendTo[return,"P_RbPercentErr"->polErr];
		viz=Show[ListPlot[{TransformPolDataForPlotting[ordPair]}],Plot[lm["Value"]*x,{x,-.1,.1}],FrameLabel->{"1/\[Delta](\!\(\*SuperscriptBox[\(GHz\), \(-1\)]\))","\[CapitalDelta]\[Theta](rad)"}];
		AppendTo[return,"Graph"->viz];
		viz=Show[ListPlot[{TransformPolDataForPlotting[ordPair]}],Plot[lm["Value"]*x,{x,-.1,.1}]];
		AppendTo[return,"GraphNoLabel"->viz];
		AppendTo[return,"Dataset Notes"->"The ordered pairs are the detuning (GHz) and difference in angle (rad) between the two pump types."];
		return
];


(* ::Subsection:: *)
(*Folder Based Processing*)


ProcessRbPolarizationRun[folder_,numberOfPumps_:3]:=Module[
{fn,density,piPump,sPlus,sMinus,r,
nVals,nAvg,nStdm,
polarizationArgs,r2,return,r3},
SetDirectory[folder];
fn=FileNames[RegularExpression["FDayScan.*_[0-9]*.dat"]];
return = ProcessRbPolarizationSequence[fn,numberOfPumps];
ResetDirectory[];
return
];


ProcessRbPolarizationSequence[fn_,numberOfPumps_:3]:=Module[
{density,piPump,sPlus,sMinus,r,
nVals,nAvg,nStdm,
polarizationArgs,r2,return,r3},
density=Take[fn,{1,-1,numberOfPumps}];

If[numberOfPumps==3,
sPlus=Take[fn,{2,-1,numberOfPumps}];
sMinus=Take[fn,{3,-1,numberOfPumps}];
,
piPump=Take[fn,{2,-1,numberOfPumps}];
sPlus=Take[fn,{3,-1,numberOfPumps}];
sMinus=Take[fn,{4,-1,numberOfPumps}];
];

r=ProcessFaradayRotationNumberDensityFile[density];
nVals=Dataset[r][All,"n_Rbx10^12"]//Normal;
nAvg=Dataset[r][Mean,"n_Rbx10^12"]//Normal;
polarizationArgs=Partition[Flatten[Riffle[Partition[fn,numberOfPumps],nVals]],numberOfPumps+1]; 
If[numberOfPumps==3,
r2=Apply[ProcessFaradayRotationPolarizationFilesThreePumps[#1,#2,#3,"PUMP",nAvg]&,polarizationArgs,1];
r3=Apply[ProcessFaradayRotationPolarizationFilesThreePumps[#1,#2,#3,"PUMP",#4]&,polarizationArgs,1];
,
r2=Apply[ProcessFaradayRotationPolarizationFilesThreePumps[#1,#3,#4,"PUMP",#5]&,polarizationArgs,1];
];
return=<|"n_Rb"->r,"P_Rb_nAvg"->r2,"P_Rb"->r3|>;
return
];


(* ::Chapter:: *)
(*Functions For Datasets*)


ReCalculatePRbWithGivenNRb[dataset_Dataset,nRb_]:=
Module[{d=dataset,densities,slopes,allPrb},
densities=d[All,"n_Rb","n_Rbx10^12"]//Normal;
slopes=d[All,"P_Rb",{"S+ Pump","S- Pump"},"Model Fit"]//Values//Normal;
allPrb=CalculateRubidiumPolarizationOnePump[slopes,nRb]
];

GeneratePTEFromPTEAnalysisFolder[folderName_]:=Module[{dr,de,pe,prb,prbNoPump,avgNRb},
dr=Import[FileNameJoin[{folderName,"P_Rb.wdx"}]];
de=Import[FileNameJoin[{folderName,"P_e.wdx"}]];
pe=de[All,Values,"processedSignal","P_e"]//Normal//Transpose;
avgNRb=(dr[All,"n_Rb","n_Rbx10^12"]//Normal//Mean)*1*^12;
prb=dr[All,"P_Rb_2",{"S+ Pump","S- Pump"},"P_Rb"]//Values//Normal//Transpose;
prb=ReCalculatePRbWithGivenNRb[dr,avgNRb["Value"]]//Transpose;
prbNoPump=ConstantArray[Around[0,.001],Length[prb[[1]]]];
prb=Join[{prbNoPump},prb];
<|"data"->Transpose[{Flatten[prb],Flatten[pe]}],"density"->avgNRb|>
];
