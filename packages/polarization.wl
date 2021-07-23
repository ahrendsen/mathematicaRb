(* ::Package:: *)

(* ::Title:: *)
(*Electron Polarization Functions*)


(* ::Chapter::Closed:: *)
(*Import other needed functions*)


Import["fileManipulation.wl"];
Import["dataManipulation.wl"];
Import["generalPlottingFunctions.wl"];
Import["apparatus.wl"];
Import["physicalConstants.wl"];


(* ::Chapter::Closed:: *)
(*Plotting Functions*)


(* ::Subchapter:: *)
(*DFT*)


(* The fourier coefficients that are input should be the output of the DFT I created *)
FourierCoefficientBarChart[fc_]:=Module[
{bcLabels,blanks,bc,
sinValues,sinValuesNoZero,
cosValues,cosValuesNoZero,
bcLabelsNoZero,bcLabelsSparse,
results},

results=<||>;
bcLabels=Drop[Keys[fc[[1]]],{2,-1,2}];
sinValues=Values[fc[[1]]]/Values[fc[[2]]][[1]];
sinValuesNoZero=Drop[sinValues,{1}];
blanks=ConstantArray["",Length[bcLabels]];
bcLabelsNoZero=Drop[bcLabels,{1}];
bcLabels=Riffle[blanks,bcLabelsNoZero];
bc=BarChart[sinValuesNoZero,ChartStyle->colorBlindPalette[[2]],ChartLabels->bcLabels,Frame->True];
AppendTo[results,"SinCoefficients"->bc];

bcLabels=Drop[Keys[fc[[2]]],{2,-1,2}];
cosValues=Values[fc[[2]]]/Values[fc[[2]]][[1]];
cosValuesNoZero=Drop[cosValues,{1}];
bcLabelsNoZero=Drop[bcLabels,{1}];
blanks=ConstantArray["",Length[bcLabelsNoZero]];
bcLabelsSparse=Riffle[blanks,bcLabelsNoZero];
bc=BarChart[cosValuesNoZero,ChartStyle->colorBlindPalette[[3]],ChartLabels->bcLabelsSparse,Frame->True];
AppendTo[results,"CosCoefficients"->bc];
results
];

(* The fourier coefficients that are input should be the output of the DFT I created *)
FourierCoefficientBarChartUnNormalized[fc_]:=Module[
{bcLabels,blanks,bc,
sinValues,sinValuesNoZero,
cosValues,cosValuesNoZero,
bcLabelsNoZero,bcLabelsSparse,
results},

results=<||>;
bcLabels=Drop[Keys[fc[[1]]],{2,-1,2}];
sinValues=Values[fc[[1]]];
sinValuesNoZero=Drop[sinValues,{1}];
blanks=ConstantArray["",Length[bcLabels]];
bcLabelsNoZero=Drop[bcLabels,{1}];
bcLabels=Riffle[blanks,bcLabelsNoZero];
bc=BarChart[sinValuesNoZero,ChartStyle->colorBlindPalette[[2]],ChartLabels->bcLabels,ImageSize->{720,450},Frame->True];
AppendTo[results,"SinCoefficients"->bc];

bcLabels=Drop[Keys[fc[[2]]],{2,-1,2}];
cosValues=Values[fc[[2]]];
bcLabelsNoZero=Drop[bcLabels,{1}];
blanks=ConstantArray["",Length[bcLabels]];
bcLabelsSparse=Riffle[blanks,bcLabels];
bc=BarChart[cosValues,ChartStyle->colorBlindPalette[[3]],ChartLabels->bcLabelsSparse,ImageSize->{720,450},Frame->True];
AppendTo[results,"CosCoefficients"->bc];
results
];

FourierFitFunction[fc_,var_]:=fc["Cos Coefficients"][0];


(* ::Subchapter::Closed:: *)
(*Raw Data*)


RepeatedMeasurementSummary[dataPoints_List]:=Module[
{dataPointsPlot,dataPointsAvg,dataPointsStd},

	dataPointsAvg=Mean[dataPoints];
	dataPointsStd=StandardDeviation[dataPoints];
	dataPointsPlot=ListPlot[dataPoints,PlotRange->{.8*dataPointsAvg,1.2*dataPointsAvg}];
	<|"currentPlot"->dataPointsPlot,"currentAvg"->dataPointsAvg,"currentStd"->dataPointsStd|>
];


(* ::Chapter:: *)
(*Polarization Functions*)


(* ::Section::Closed:: *)
(*Constants*)


(* ::Input::Initialization:: *)
(*stokesNames={"p0","p1","p2","p3_mag","p3_c2","p3_s2"};*)
stokesNames={"p0","p1","p2","p3"};
(*stokesErrNames={"p0err","p1err","p2err","p3_magerr","p3_c2err","p3_s2err"};*)
stokesErrNames={"p0_err","p1_err","p2_err","p3_err"};
nA=1*^-9;
mTorr=1*^-3;
(* These are defined in the physicalConstansts.wl file now.
alpha=20.4*\[Pi]/180;
beta=68.7*\[Pi]/180;
deltaAB=alpha-beta;
delta=94.54*\[Pi]/180 (*Munir's reported 1.66 rad plus or minus .01*);
*)


(* These are defined in the physicalConstants.wl file now.
alphaOld=3.2*\[Pi]/180;
betaOld=-6.95*\[Pi]/180;
rev=1;
*)


(* ::Section::Closed:: *)
(*Calculate Polarization*)


FourierFit[list_,rev_]:=Module[{fit,function,dataPts,stepSize,dataPtsPerRev,angles,data(*,C0,C2,C4,S2,S4*)},
	function[\[Theta]_]= C0+C2 Cos[2(\[Theta]+\[Theta]o)]+C4 Cos[4(\[Theta]+\[Theta]o)]+S2 Sin[2(\[Theta]+\[Theta]o)]+S4 Sin[4(\[Theta]+\[Theta]o)];
	dataPts=Length[list];
	dataPtsPerRev=dataPts/rev;
	angles=Range[0,(2\[Pi]*rev)-2\[Pi]/dataPtsPerRev,2\[Pi]/dataPtsPerRev];
	data=Transpose[{list,angles}];
	fit=NonlinearModelFit[data,function[\[Theta]],{{C0,1000},C2,C4,S2,S4,\[Theta]o},\[Theta]]
];


CalculateElectronPolarization[stokes_Association]:=Module[{p1,p3},
p3=stokes["p3"];
p1=stokes["p1"];
CalculateElectronPolarization[p1,p3]
];

CalculateElectronPolarization[p1_,p3_]:=p3*(polConstxxp3Coeff/(polConstxxp1Const+polConstxxp1Coeff*p1))


CalculateElectronPolarizationError[stokes_,averageStokes_]:=Module[{pe,peError,p1Coeff,p3Coeff,p1Const,p1Values,p3Values},
pe=CalculateElectronPolarization[averageStokes];
p1Const=polConstxxp1Const;
p1Coeff=polConstxxp1Coeff;
p3Coeff=polConstxxp3Coeff;
p1Values=Transpose[stokes][[2]]; (*We take the second item in the transpose because the first item is p0 *)
p3Values=Transpose[stokes][[4]]; (*The 4th item is p3 *)
peError=Sqrt[pe^2*((averageStokes["p3_err"]/averageStokes["p3"])^2+((averageStokes["p1_err"]*p1Coeff)/(p1Const+p1Coeff*averageStokes["p1"]))^2-Covariance[p1Values*p1Coeff+p1Const,p3Values*p3Coeff]^2/(p3Coeff*averageStokes["p3"]*(p1Const+p1Coeff*averageStokes["p1"])))]
];


CalculateAverageStokesFromFilenames[polFileNames_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{signal,stokes,stokesStdDev,file,rev,j,stokesValues,transposedValues,singleRunInfo,allSingleRunInfo,current,counts,countsPlot,combinedPlot,results,averageStokes},
stokes=<|"p0"->0,"p1"->0,"p2"->0,"p3_mag"->0,"p3_c2"->0,"p3_s2"->0|>;
stokesStdDev={0,0,0,0,0,0};
stokesValues={};
allSingleRunInfo=<||>;
For[j=1,j<=Length[polFileNames],j++,
	stokes=CalculateOneFileStokes[polFileNames[[j]],alpha,beta,delta];
	current=CalculateOneFileCurrent[polFileNames[[j]]];
	counts=GetIntensityArrayFromFileName[polFileNames[[j]]];
	countsPlot=ListPlot[counts];
	AppendTo[stokesValues,stokes];
	singleRunInfo=<|"stokes"->stokes|>;
	AppendTo[singleRunInfo,current];
	AppendTo[singleRunInfo,"counts"->countsPlot];
	AppendTo[singleRunInfo,"rawCounts"->counts];
	AppendTo[allSingleRunInfo, polFileNames[[j]]->singleRunInfo];
];
combinedPlot=ListPlot[GetCurrentNormalizedAverageCountRateFromFileNames[polFileNames][[1]],PlotLabel->"Multiple Runs Average Counts",
																			FrameLabel->{"Stepper Motor Position Number","Current Normalized Counts (Hz/nA)"}];
combinedPlot=ListPlot[GetCurrentNormalizedAverageCountRateFromFileNames[polFileNames][[1]]];
averageStokes=AverageStokes[stokesValues];

results=<|"P_e"->CalculateElectronPolarization[averageStokes],"P_eerr"->CalculateElectronPolarizationError[stokesValues,averageStokes]|>;
AppendTo[results,"timeStart"->ExtractTimeInfoFromFileNameString[polFileNames[[1]]]];
AppendTo[results,"timeEnd"->ExtractTimeInfoFromFileNameString[polFileNames[[-1]]]];
AppendTo[results,averageStokes];
AppendTo[results,ImportFile[polFileNames[[1]]][[1]]]; (* The header of the first file. *)
AppendTo[results,<|"IndividualRunInfo"->allSingleRunInfo,"AllRunsPlot"->combinedPlot|>]
]

CalculateNoBeamAverageStokesFromFilenames[polFileNames_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{signal,stokes,stokesStdDev,file,rev,j,stokesValues,transposedValues,singleRunInfo,allSingleRunInfo,current,counts,countsPlot,combinedPlot,results,averageStokes},
stokes=<|"p0"->0,"p1"->0,"p2"->0,"p3_mag"->0,"p3_c2"->0,"p3_s2"->0|>;
stokesStdDev={0,0,0,0,0,0};
stokesValues={};
allSingleRunInfo=<||>;
For[j=1,j<=Length[polFileNames],j++,
	stokes=CalculateOneFileStokes[polFileNames[[j]],alpha,beta,delta];
	current=CalculateOneFileCurrent[polFileNames[[j]]];
	counts=GetIntensityArrayFromFileName[polFileNames[[j]]];
	countsPlot=ListPlot[counts];
	AppendTo[stokesValues,stokes];
	singleRunInfo=<|"stokes"->stokes|>;
	AppendTo[singleRunInfo,current];
	AppendTo[singleRunInfo,"counts"->countsPlot];
	AppendTo[singleRunInfo,"rawCounts"->counts];
	AppendTo[allSingleRunInfo, polFileNames[[j]]->singleRunInfo];
];
combinedPlot=ListPlot[GetAverageCountRateFromFileNames[polFileNames][[1]],PlotLabel->"Multiple Runs Average Counts",
																			FrameLabel->{"Stepper Motor Position Number","Current Normalized Counts (Hz/nA)"}];
combinedPlot=ListPlot[GetAverageCountRateFromFileNames[polFileNames][[1]]];
averageStokes=AverageStokes[stokesValues];

results=<|"P_e"->CalculateElectronPolarization[averageStokes],"P_eerr"->CalculateElectronPolarizationError[stokesValues,averageStokes]|>;
AppendTo[results,"timeStart"->ExtractTimeInfoFromFileNameString[polFileNames[[1]]]];
AppendTo[results,"timeEnd"->ExtractTimeInfoFromFileNameString[polFileNames[[-1]]]];
AppendTo[results,averageStokes];
AppendTo[results,ImportFile[polFileNames[[1]]][[1]]]; (* The header of the first file. *)
AppendTo[results,<|"IndividualRunInfo"->allSingleRunInfo,"AllRunsPlot"->combinedPlot|>]
]

CalculateAverageCurrentFromFiles[polFileNames_]:=Module[{signal,stokes,stokesStdDev,file,rev,j,stokesValues,transposedValues},
stokes=<|"p0"->0,"p1"->0,"p2"->0,"p3_mag"->0,"p3_c2"->0,"p3_s2"->0|>;
stokesStdDev={0,0,0,0,0,0};
stokesValues={};
For[j=1,j<=Length[polFileNames],j++,
stokes=CalculateOneFileCurrent[polFileNames[[j]]];
AppendTo[stokesValues,stokes]
];
AverageStokes[stokesValues]
]


StokesNormalize[stokesVector_]:=Module[{intensity,newVector,i},
newVector=stokesVector;
intensity=stokesVector[[1]];
For[i=2,i<=Length[stokesVector],i++,
newVector[[i]]=stokesVector[[i]]/intensity;
];
newVector
];

StokesUnNormalize[stokesVector_]:=Module[{intensity,newVector,i},
newVector=stokesVector;
intensity=stokesVector[[1]];
For[i=2,i<=Length[stokesVector],i++,
newVector[[i]]=stokesVector[[i]]*intensity;
];
newVector
];

StokesSubtract[addedStokes_,subtractedStokes_]:=Module[
{unNormAdd,unNormSub},
unNormAdd=UnNormalizeStokes[addedStokes];
unNormSub=UnNormalizeStokes[subtractedStokes];
NormalizeStokes[unNormAdd-unNormSub]
];

StokesDummyCheck[stokesVector_]:=Module[
{polarizedFraction,
linearFraction,
circularFraction},
polarizedFraction=Norm[Values[stokesVector][[2;;]]];
linearFraction=Norm[Values[stokesVector][[2;;3]]];
circularFraction=Norm[Values[stokesVector][[4;;4]]];
<|"polarizedFraction"->polarizedFraction,"linearFraction"->linearFraction,"circularFraction"->circularFraction|>
];

GetIntensityArrayFromFileName[fileName_]:=Module[{f,intensityArray},
f=ImportFile[fileName];
intensityArray=Normal[f[[2]][All,{"ANGLE","COUNT"}]]
];

GetCurrentArrayFromFileName[fileName_]:=Module[{f,intensityArray},
f=ImportFile[fileName];
intensityArray=Normal[f[[2]][All,{"ANGLE","CURRENT"}]]
];

PlotRawPolFromFileName[fileName_]:=Module[
{f,intensityArray},
f=ImportFile[fileName];
intensityArray=Normal[f[[2]][All,"COUNT"]];
ListPlot[intensityArray]
];
SetAttributes[PlotRawPolFromFileName,Listable];

PlotListOfPolFileNames[fileNames_,plotRange_:{Automatic,Automatic}]:=Module[{intensityArrays,i},
intensityArrays={};
For[i=1,i<=Length[fileNames],i++,
AppendTo[intensityArrays,Legended[GetIntensityArrayFromFileName[fileNames[[i]]],fileNames[[i]]]];
];
ListPlot[intensityArrays,PlotRange->plotRange]
];

GetStokesFromFileNames[fileNames_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{i,stokes,signal,stokesData},
stokesData={};
signal=ImportFile[fileNames];
For[i=1,i<=Length[fileNames],i++,
stokes=StokesParametersFromPolFile[signal[[i]],alpha,beta,delta,signal[[i]][[1]]["REV"]];
AppendTo[stokesData,stokes];
];
stokesData
];

GetAverageFourierCoefficientsFromFileNames[fileNames_]:=Module[{fc},
fc=GetFourierCoefficientsFromFileNames[fileNames]];

GetFourierCoefficientsFromFileNames[fileNames_]:=Module[{files,i,coefficientsList,intensityArray,revs},
files=ImportFile[fileNames];
coefficientsList={};
For[i=1,i<=Length[files],i++,
revs=files[[i]][[1]]["REV"];
intensityArray=Normal[files[[i]][[2]][All,"COUNT"]];
AppendTo[coefficientsList,DFT[intensityArray,revs]];
];
coefficientsList
];

GetAoutToEnergyModelFromExcitationFunctionTimeStamp[excitationFileNames_,timeStamp_]:=Module[{exFileIndex,exFile,aoutToElectronEnergyData,x},
exFileIndex=GetIndexOfFilenameFromTimestamp[excitationFileNames,timeStamp];
exFile=ImportFile[excitationFileNames[[exFileIndex]]];
aoutToElectronEnergyData=Normal[exFile[[2]][All,{"Aout","SecondaryElectronEnergy"}][[All,Values]]];
LinearModelFit[aoutToElectronEnergyData,x,x]
];

GetListOfPropertiesFromHeadersOfListOfFileNames[fileName_,property_]:=Module[{f},
f=ImportFile[fileName];
Normal[f[[1]][property]]
];

SetAttributes[GetListOfPropertiesFromHeadersOfListOfFileNames,Listable];

OrganizeFilesByCategories[fileNames_,startFile_,categories_,runs_]:=Module[{allFiles,j,i,listOfNames,subsetFiles,firstFile,lastFile,files},
allFiles=<||>;
subsetFiles=<||>;
files={};

For[j=1,j<=Length[categories[[1]]],j++,
For[i=1,i<=Length[categories[[2]]],i++,
firstFile=((startFile-1)+i)+(j-1)*Length[categories[[2]]]*runs;
lastFile=((startFile-1)+i+(j)*Length[categories[[2]]]*runs)-1;
listOfNames=Take[fileNames,{firstFile,lastFile,Length[categories[[2]]]}];
AppendTo[files,<|"cat1"-> categories[[1]][[j]],"cat2"->categories[[2]][[i]],"fileNames"-> listOfNames|>];
];
];
files
];

GetStokesFromDatasetLine[datasetLine_]:=Module[{stokes,stokesErr,v},
stokes=<||>;
stokesErr=<||>;
For[v=1,v<=Length[stokesNames],v++,
AppendTo[stokes,stokesNames[[v]]->Normal[datasetLine[All,stokesNames[[v]]]][[1]]];
AppendTo[stokesErr,stokesErrNames[[v]]->Normal[datasetLine[All,stokesErrNames[[v]]]][[1]]];
];
{stokes,stokesErr}
];


(* ::Input::Initialization:: *)
OrganizeFscanFilesByCategories[fileNames_,startFile_,categories_,runs_]:=Module[{allFiles,j,i,listOfNames,subsetFiles,firstFile,lastFile,files},
allFiles=<||>;
subsetFiles=<||>;
files={};

For[j=1,j<=Length[categories[[1]]],j++,
For[i=1,i<=Length[categories[[2]]],i++,
firstFile=(startFile+(i-1)*runs)+(j-1)*Length[categories[[2]]]*runs;
AppendTo[files,<|Keys[categories][[1]]-> categories[[1]][[j]],Keys[categories][[2]]->categories[[2]][[i]],"fileNames"->fileNames[[firstFile]]|>];
];
];
files
];

ExtractTimeInfoFromFileNameString[fileNameString_]:=Module[{lastChars,date,year,month,day,hour,min,sec,do},
lastChars=StringTake[fileNameString,-21];
date=StringTake[lastChars,17];
year=StringTake[date,{1,4}];
month=StringTake[date,{6,7}];
day=StringTake[date,{9,10}];
hour=StringTake[date,{12,13}];
min=StringTake[date,{14,15}];
sec=StringTake[date,{16,17}];
do=DateObject[{Interpreter["Number"][year],Interpreter["Number"][month],Interpreter["Number"][day],Interpreter["Number"][hour],Interpreter["Number"][min],Interpreter["Number"][sec]}]
];


(* ::Section::Closed:: *)
(*Quick Polarization*)


QPOLProcessFileDataNoBackground[fileName_]:=Module[
{f,avg,stdDev,
t,p,
normCounts,
p0,p3,pe,peErr,results},
results=<||>;
f=ImportFile[fileName];
t=ToExpression[f[[1]]["DWELL(s)"]];
normCounts=f[[2]][All,<|"COUNT+45"->#["COUNT+45"]/#["CURRENT+45"]/(t),"COUNT-45"->#["COUNT-45"]/#["CURRENT-45"]/(t)|>&];
p3=Normal[Flatten[normCounts[All,{(#["COUNT+45"]-#["COUNT-45"])/(#["COUNT+45"]+#["COUNT-45"])}&]]];
p0=Normal[Flatten[normCounts[All,{(#["COUNT+45"]-#["COUNT-45"])/(#["COUNT+45"]+#["COUNT-45"])}&]]];
pe=CalculateElectronPolarization[polConstxxunPolP1,Mean[p3]];
peErr=CalculateElectronPolarization[polConstxxunPolP1,StandardDeviation[p3]/Sqrt[Length[p3]]];
AppendTo[results,<|"QPolFile"->fileName,"P_e"->pe,"P_eerr"->peErr,"rawP3"->p3,"averageP3"->Mean[p3],"stdDevP3"->StandardDeviation[p3],"p1"->polConstxxunPolP1|>];
AppendTo[results,<|"AvgCurrent+45"->f[[2]][Mean,#["CURRENT+45"]&],
					"AvgCurrent-45"->f[[2]][Mean,#["CURRENT-45"]&],
					"AvgCount+45"->N[f[[2]][Mean,#["COUNT+45"]&]],
					"AvgCount-45"->N[f[[2]][Mean,#["COUNT-45"]&]]|>];
AppendTo[results,f[[1]]];
results
]
SetAttributes[QPOLProcessFileDataNoBackground,Listable]

QPOLProcessFileDataNoBeamBackground[fileName_,noBeamFileName_]:=Module[
{t,f,avgNoBeam,avg,stdDev,i45,i135,i45err,i135err,p3,pe,peErr,results,fNoBeam,p,bnc(*background Normalized Counts*)},
results=<||>;
fNoBeam=ImportFile[noBeamFileName];
t=ToExpression[fNoBeam[[1]]["DWELL(s)"]];
avgNoBeam=fNoBeam[[2]][Mean,<|"COUNT+45"->#["COUNT+45"]/t,"COUNT-45"->#["COUNT-45"]/t|>&];
f=ImportFile[fileName];
t=ToExpression[f[[1]]["DWELL(s)"]];
bnc=f[[2]][All,<|"COUNT+45"->(#["COUNT+45"]/t-avgNoBeam["COUNT+45"])/(#["CURRENT+45"]),"COUNT-45"->(#["COUNT-45"]/t-avgNoBeam["COUNT-45"])/(#["CURRENT-45"])|>&];
p3=Normal[Flatten[bnc[All,{(#["COUNT+45"]-#["COUNT-45"])/(#["COUNT+45"]+#["COUNT-45"])}&]]];
pe=CalculateElectronPolarization[polConstxxunPolP1,Mean[p3]];
peErr=CalculateElectronPolarization[polConstxxunPolP1,StandardDeviation[p3]/Sqrt[Length[p3]]];
AppendTo[results,<|"P_e"->pe,"P_eerr"->peErr,"rawP3"->p3,"averageP3"->Mean[p3],"stdDevP3"->StandardDeviation[p3],"p1"->polConstxxunPolP1|>];
results
]

QPOLProcessFileDataFullBackground[fileName_,noBeamFileName_,belowThresholdFileName_]:=Module[
{t,p,f,avgNoBeam,avgBT,avg,stdDev,i45,i135,i45err,i135err,p3,results,fNoBeam,btf(*Below threshold File*),bnc(*background Normalized Counts*),pe,peErr},
results=<||>;
fNoBeam=ImportFile[noBeamFileName];
t=ToExpression[fNoBeam[[1]]["DWELL(s)"]];
avgNoBeam=fNoBeam[[2]][Mean,<|"COUNT+45"->#["COUNT+45"]/t,"COUNT-45"->#["COUNT-45"]/t|>&];
btf=ImportFile[belowThresholdFileName];
t=ToExpression[btf[[1]]["DWELL(s)"]];
avgBT=btf[[2]][Mean,<|"COUNT+45"->(#["COUNT+45"]/t-avgNoBeam["COUNT+45"])/(#["CURRENT+45"]),"COUNT-45"->(#["COUNT-45"]/t-avgNoBeam["COUNT-45"])/(#["CURRENT-45"])|>&];
f=ImportFile[fileName];
t=ToExpression[f[[1]]["DWELL(s)"]];
bnc=f[[2]][All,<|"COUNT+45"->(#["COUNT+45"]/t-avgNoBeam["COUNT+45"])/(#["CURRENT+45"])-avgBT["COUNT+45"],"COUNT-45"->(#["COUNT-45"]/t-avgNoBeam["COUNT-45"])/(#["CURRENT-45"])-avgBT["COUNT+45"]|>&];
avg=bnc[Mean,{"COUNT+45","COUNT-45"}];
p3=Normal[Flatten[bnc[All,{(#["COUNT+45"]-#["COUNT-45"])/(#["COUNT+45"]+#["COUNT-45"])}&]]];
pe=CalculateElectronPolarization[polConstxxunPolP1,Mean[p3]];
peErr=CalculateElectronPolarization[polConstxxunPolP1,StandardDeviation[p3]/Sqrt[Length[p3]]];
AppendTo[results,<|"P_e"->pe,"P_eerr"->peErr,"rawP3"->p3,"averageP3"->Mean[p3],"stdDevP3"->StandardDeviation[p3],"p1"->polConstxxunPolP1|>];
AppendTo[results,<|"AvgCurrent+45"->f[[2]][Mean,#["CURRENT+45"]&],
					"AvgCurrent-45"->f[[2]][Mean,#["CURRENT-45"]&],
					"AvgCount+45"->N[f[[2]][Mean,#["COUNT+45"]&]],
					"AvgCount-45"->N[f[[2]][Mean,#["COUNT-45"]&]]|>];
results
]


(* ::Section::Closed:: *)
(*Public Functions *)


(* ::Subsection:: *)
(*Single File*)


(* The functions that will be regularly called from the Data Processing notebooks. *)
ProcessElectronPolarizationFile[fn_String,subtractDarkCounts_:True,cn_:True(*Current Normalization*)]:=Module[
{counts,current,dwellTime,header,signal,
f,dataset,results,error},
	results=<||>;
	f=ImportFile[fn];
	header=f[[1]];
	AppendTo[results,header];
	dataset=f[[2]];
	AppendTo[results,"dataset"->Normal[dataset]];

	counts=Normal[dataset[All,"COUNT"]];
	current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
	dwellTime=header["DWELL(s)"];
	If[cn==True,
		signal=GetCurrentNormalizedCountRate[counts,current,dwellTime,subtractDarkCounts];
		error=GetCurrentNormalizedCountRateError[counts,current,dwellTime,subtractDarkCounts];,
		signal=GetCountRate[counts,dwellTime,subtractDarkCounts];
		error=GetCountRateError[counts,dwellTime,subtractDarkCounts];
	];
	
	AppendTo[results,"signal"->signal];
	AppendTo[results,"error"->error];
	AppendTo[results,ProcessElectronPolarizationFromSignal[signal]](* See Private Functions for the rest of processing. *)
];

ProcessElectronPolarizationFileSubtractBackground[fn_String,fnBack_]:=Module[
{counts,current,dwellTime,header,signal,
f,dataset,results,background,d},
	results=<||>;
	
	f=ImportFile[fn];
	header=f[[1]];
	AppendTo[results,header];
	dataset=f[[2]];
	AppendTo[results,"dataset"->Normal[dataset]];

	counts=Normal[dataset[All,"COUNT"]];
	current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
	dwellTime=header["DWELL(s)"];
	signal=GetCurrentNormalizedCountRate[counts,current,dwellTime];
	AppendTo[results,"signal"->signal];

	If[Head[fnBack]!= List,
		f=ImportFile[fnBack];
		header=f[[1]];
		AppendTo[results,"backgroundHeader"->header];
		dataset=f[[2]];
		AppendTo[results,"datasetBack"->Normal[dataset]];

		counts=Normal[dataset[All,"COUNT"]];
		current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
		dwellTime=header["DWELL(s)"];
		background=GetCurrentNormalizedCountRate[counts,current,dwellTime];
		,
		d=ProcessElectronPolarizationFileAverage[fnBack,True,True];
		background=Normal[Dataset[d][Values,"avgSignal"]][[1]];
	];

	AppendTo[results,"background"->background];
	AppendTo[results,"trueSignal"->signal-background];
	AppendTo[results,"MeanCurrent"->Mean[current]];
	AppendTo[results,"MeanCounts"->N[Mean[counts],3]];
	AppendTo[results,"MeanCountRate"->Mean[signal]];
	AppendTo[results,ProcessElectronPolarizationFromSignal[signal-background]](* See Private Functions for the rest of processing. *)
];


(* ::Subsection:: *)
(*Multiple Files*)


ProcessElectronPolarizationFileAverage[fn_List,darkCountsSubtract_:True,cn_:True]:=Module[
{},
ProcessElectronPolarizationFileAverageSUM[fn,darkCountsSubtract,cn]
];

ProcessElectronPolarizationFileSubtractBackgroundAverage[fn_List,fnBack_List,darkCountsSubtract_:True,cn_:True]:=Module[
{},
GetAverageElectronPolarizationWithBackgroudSubtractionSUM[fn,fnBack,cn]
];


(* ::Section:: *)
(*Private Functions*)


(* ::Subsection:: *)
(*Single File*)


(* The "helper functions" that the public functions will access to make their job easier. Won't typically be called from the notebook.*)

(* The old way. Abandoned in favor of adding in an uncertainty in the counts. 
GetCurrentNormalizedCountRate[counts_,current_,dwellTime_,darkCountsSubtract_:True]:=Module[
{countRate,r},
	If[darkCountsSubtract==True,
		countRate=(counts/dwellTime)-darkCountRate;,
	
		countRate=(counts/dwellTime);
	];
	r=countRate/current; (* r= Count rate *)
	Normal[r]
];
*)

GetCurrentNormalizedCountRate[counts_,current_,dwellTime_,darkCountsSubtract_:True]:=Module[
{countRate,r,cts,nA=1*^-9},
	cts=Apply[Around,{#,Sqrt[#]}&/@counts,{1}];
	If[darkCountsSubtract==True,
		countRate=(cts/dwellTime)-darkCountRate;,
	
		countRate=(cts/dwellTime);
	];
	r=countRate/current; (* r= Count rate *)
	Normal[r]
];


GetCurrentNormalizedCountRateError[counts_,current_,dwellTime_,darkCountsSubtract_:True]:=Module[
{countRateError,r},
	If[darkCountsSubtract==True,
		countRateError=Sqrt[(Sqrt[counts]/dwellTime)^2+darkCountRate];,
		countRateError=Sqrt[(Sqrt[counts]/dwellTime)^2];
	];
	r=countRateError/current; (* r= Count rate *)
	Normal[r]/Sqrt[Length[r]]
];

GetCountRate[counts_,dwellTime_,darkCountsSubtract_:True]:=Module[
{countRate,cts,r},
	cts=Apply[Around,{#,Sqrt[#]}&/@counts,{1}];
	If[darkCountsSubtract==True,
		countRate=(cts/dwellTime)-darkCountRate;,
		countRate=(cts/dwellTime);
	];
	r=countRate; (* r= Count rate *)
	N[Normal[r]]
];

GetCountRateError[counts_,dwellTime_,darkCountsSubtract_:True]:=Module[
{countRateError,r},
	If[darkCountsSubtract==True,
		countRateError=Sqrt[(Sqrt[counts]/dwellTime)^2+darkCountRate];,
		countRateError=Sqrt[(Sqrt[counts]/dwellTime)^2];
	];
	r=countRateError; (* r= Count rate *)
	Normal[r]/Sqrt[Length[r]]
];

CalculateStokesFromDFT[dftOutput_,rev_:rev,alpha_:alpha,beta0_:beta,delta_:delta]:=Module[{c,s},
c=dftOutput["Cos Coefficients"];
s=dftOutput["Sin Coefficients"];
CalculateStokesFromFourierCoefficients[c[0],c[2],c[4],s[2],s[4],rev,alpha,beta0,delta]
];

(* 2021-04-15: After working so hard to determine the positions of the 
optical elements, I only get a value of P2=0 when my beta0 is added rather than subtracted, 
though Keith and Will B. seem to indicate that beta0 should be subtracted. 
I'm going to go with the physics here and use my new values and beta0 added. *)
Clear[CalculateStokesFromFourierCoefficients];
CalculateStokesFromFourierCoefficients[c0_,c2_,c4_,s2_,s4_,rev_:rev,alpha_:alpha,beta0_:beta,delta_:delta]:=Module[
{stokes,sign},
	stokes=<||>;
	sign=Sign[c2["Value"]];
	AppendTo[stokes,stokesNames[[1]]->c0-(1+Cos[delta])/(1-Cos[delta])*(c4*Cos[4 alpha+4 beta0]+s4*Sin[4 alpha+4 beta0])];
	AppendTo[stokes,stokesNames[[2]]->2/(1-Cos[delta])*(c4*Cos[2*alpha + 4 *beta0]+s4*Sin[2*alpha + 4 *beta0])/stokes[[1]]];
	AppendTo[stokes,stokesNames[[3]]->2/(1-Cos[delta])*(s4*Cos[2*alpha + 4 *beta0]-c4*Sin[2*alpha + 4 *beta0])/stokes[[1]]];
	AppendTo[stokes,stokesNames[[4]]-> sign*Sqrt[c2^2+s2^2]/(Sin[delta]^2*stokes[[1]])];
	(*AppendTo[stokes,stokesNames[[5]]->c2/(Sin[delta]Sin[2 alpha + 2 beta0]*stokes[[1]])]; 
	AppendTo[stokes,stokesNames[[4]]->-s2/(Sin[delta]Cos[2 alpha + 2beta0]*stokes[[1]])];*)
	stokes
];

ProcessElectronPolarizationFromSignal[signal_List]:=Module[
{fc,s,c,plot,ePlot,fourierFit,fourierFitFunction,
stokes,electronPolarization,dataPts,
results,err},
results=<||>;
fc=DFT[signal];
dataPts=Length[signal];
s=fc["Sin Coefficients"];
c=fc["Cos Coefficients"];
ePlot=ListPlot[{Partition[Riffle[Range[0,2\[Pi]-2\[Pi]/dataPts,2\[Pi]/dataPts],signal],2]},
{PlotLegends->Placed[{"Data"},{.5,.5}],
FrameTicks->{{Automatic,None},{{0,\[Pi]/2,\[Pi],3\[Pi]/2,2\[Pi]},None}},
GridLines->{{0,\[Pi]/2,\[Pi],3\[Pi]/2,2\[Pi]},Automatic}}];
fourierFitFunction=c[0]+c[2]*Cos[2#]+s[2]*Sin[2#]+c[4]*Cos[4#]+s[4]*Sin[4#]&;
plot=Plot[{
fourierFitFunction[\[Theta]]["Value"],
fourierFitFunction[\[Theta]]["Value"]+fourierFitFunction[\[Theta]]["Uncertainty"],
fourierFitFunction[\[Theta]]["Value"]-fourierFitFunction[\[Theta]]["Uncertainty"]},{\[Theta],0,2\[Pi]},
PlotRange->Full,
PlotLegends->Placed[{"Fit"},{.5,.5}],
FrameTicks->{{Automatic,None},{{0,\[Pi]/2,\[Pi],3\[Pi]/2,2\[Pi]},None}},
GridLines->{Automatic,{0,\[Pi]/2,\[Pi],3\[Pi]/2,2\[Pi]}},
Filling->{1->{2},1->{3}},
FillingStyle->RGBColor["#6699CC"],
PlotStyle->{RGBColor["#004488"],{Opacity[0],White},{Opacity[0],White}}];
fourierFit=Show[plot,ePlot,ImageSize->{UpTo[3*72],UpTo[2*72]}];
AppendTo[results,"rawAverageSignal"->ListPlot[{signal},ImageSize->{UpTo[3*72],UpTo[2*72]}]];
AppendTo[results,"fourierFit"->fourierFit];
AppendTo[results,"fourierFitFunction"->fourierFitFunction];
AppendTo[results,"fc"->fc];
AppendTo[results,"fcBarChart"->FourierCoefficientBarChart[fc]];
AppendTo[results,ProcessElectronPolarizationFromFourierCoefficients[fc]]
];

ProcessElectronPolarizationFromFourierCoefficients[fc_]:=Module[
{s,c,results,stokes},
	results=<||>;
	s=fc["Sin Coefficients"];
	c=fc["Cos Coefficients"];
	stokes=CalculateStokesFromFourierCoefficients[c[0],c[2],c[4],s[2],s[4]];
	AppendTo[results,"stokes"->stokes];
	AppendTo[results,ProcessElectronPolarizationFromStokesParameters[stokes]]
];

ProcessElectronPolarizationFromStokesParameters[stokes_Association]:=Module[
{results,electronPolarization},
	results=<||>;
	electronPolarization=CalculateElectronPolarization[stokes["p1"],stokes["p3"]];
	AppendTo[results,"P_e"->electronPolarization]
];


(* ::Subsection::Closed:: *)
(*Multiple Files*)


(* ::Subsubsection::Closed:: *)
(*Obtaining Count Rates*)


GetAverageCurrentNormalizedCountRateFromFilenames[fn_,darkCountsSubtract_:True]:=Module[
	{signals,i,f,header,dataset,
	datasets,currents,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,singleRunInfo},
	results=<||>;
	f=ImportFile[fn];

	signals={};
	datasets={};
	currents={};
	fileNames={};
	For[i=1,i<=Length[f],i++,
		header=f[[i]][[1]];
		dataset=f[[i]][[2]];
		AppendTo[datasets,Normal[dataset]];
		counts=Normal[dataset[All,"COUNT"]];
		current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
		AppendTo[currents,RepeatedMeasurementSummary[current]];
		dwellTime=header["DWELL(s)"];
		signal=GetCurrentNormalizedCountRate[counts,current,dwellTime,darkCountsSubtract];
		AppendTo[signals,signal];
		AppendTo[fileNames,header["File"]];
	];
	AppendTo[results,"averagedFiles"->fn];
	
	singleRunInfo=Transpose[{signals,datasets}];
	threader[arg_]:=AssociationThread[{"signal","dataset"},arg];
	singleRunInfo=AssociationThread[fn,Map[threader,singleRunInfo]];
	
	For[i=1,i<=Length[singleRunInfo],i++,
	AppendTo[singleRunInfo[[i]],currents[[i]]];
	AppendTo[singleRunInfo[[i]],"Time"->GetTimeInfoFromFileNameString[fn[[i]]]];
	AppendTo[singleRunInfo[[i]],header];
	];

	signal=Mean[signals];
	signalStdm=StandardDeviation[signals]/Sqrt[Length[signals]];
	signal=Apply[Around[#1,#2]&,Partition[Riffle[signal,signalStdm],2],{1}];
	AppendTo[results,"avgSignal"->signal];
	AppendTo[results,"Time"->GetTimeInfoFromFileNameString[fn[[1]]]];
	
	AppendTo[results,"IndividualRunInfo"->Dataset[singleRunInfo]]
];

Clear[GetAverageCountRateFromFilenames]
GetAverageCountRateFromFilenames[fn_,darkCountsSubtract_:True]:=Module[{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,datasets,
	currents,singleRunInfo},
	results=<||>;
	f=ImportFile[fn];

	signals={};
	fileNames={};
	datasets={};
	currents={};
	
	(* Extract all the information from each of the files *)
	For[i=1,i<=Length[f],i++,
		header=f[[i]][[1]];
		dataset=f[[i]][[2]];
		AppendTo[datasets,Normal[dataset]];
		counts=Normal[dataset[All,"COUNT"]];
		dwellTime=header["DWELL(s)"];
		signal=GetCountRate[counts,dwellTime,darkCountsSubtract];
		AppendTo[signals,signal];
		AppendTo[fileNames,header["File"]];
	];
	
	(* Save the information so that it can be investigated later *)
	AppendTo[results,"averagedFiles"->Map[FileNameTake[#]&,fn]];
	singleRunInfo=Transpose[{signals,datasets}];
	threader[arg_]:=AssociationThread[{"signal","dataset"},arg];
	singleRunInfo=AssociationThread[fn,Map[threader,singleRunInfo]];
	
	(* Average the files *)
	signal=Mean[signals];
	signalStdm=StandardDeviation[signals]/Sqrt[Length[signals]];
	(* Use the Around[] object to let Mathematica take care of error propogation *)
	signal=Apply[Around[#1,#2]&,Transpose[{signal,signalStdm}],{1}];
	AppendTo[results,"avgSignal"->signal];
	AppendTo[results,"IndividualRunInfo"->singleRunInfo]
];


(* ::Subsubsection::Closed:: *)
(*Averaging Methods of Files*)


Clear[ProcessElectronPolarizationFileAverageSUM];
ProcessElectronPolarizationFileAverageSUM[fn_List,darkCountsSubtract_:True,cn_:True(*Current Normalization*)]:=Module[
	{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,r},
	results=<||>;
	
	If[cn==True,
	AppendTo[results,r=GetAverageCurrentNormalizedCountRateFromFilenames[fn,darkCountsSubtract]];,
	AppendTo[results,r=GetAverageCountRateFromFilenames[fn,darkCountsSubtract]];
	];

	AppendTo[results,ProcessElectronPolarizationFromSignal[r["avgSignal"]]];
	AppendTo[results,"IndividualRunInfo"->Dataset[r]["IndividualRunInfo"]];
	<|FileNameTake[fn[[1]]]->results|>
];

(* This takes as input a list of file names, and outputs a single electron polarization with the
errors calculated by computing the Fourier Coefficients for each dataset and taking the average of 
these values as well as the standard deviation of the mean, which will then be used as the error of
the Fourier Coefficients and will be propogated to calculate the error in the polarization. *)
Clear[ProcessElectronPolarizationFileAverageFOURIER];
ProcessElectronPolarizationFileAverageFOURIER[fn_List,cn_:True,darkCountsSubtract_:True]:=Module[
	{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,dftResults,dfts,
	fcTitles,res,errors,avgFourierComponents,
	irResults (*individual run results*),iri(* invidual run info *),r},
	results=<||>;
	irResults=<||>;

	r=ProcessElectronPolarizationFileAverageSUM[fn,darkCountsSubtract,cn];

	iri=Dataset[r][Values,"IndividualRunInfo"][[1]];

	signals=Normal[Dataset[r][Values,"IndividualRunInfo",Values,"signal"][[1]]];
	irResultsMapper[results_]:=AssociationThread[{"fcCoefficients"},results];
	
	dfts=DFT/@signals;
	fcTitles={"Sin Coefficients","Cos Coefficients"};
	res=Map[StandardDeviation[Values[dfts[[All,#]]]]/Length[dfts]&,fcTitles];
	fcThreader[values_]:=AssociationThread[Keys[dfts[[1,1]]],values];
	res=Map[fcThreader,res];
	errors=AssociationThread[StringInsert[fcTitles," Error",-1],res];
	avgFourierComponents=Join[Mean[dfts],errors];
	AppendTo[irResults,AssociationThread[fn,Map[irResultsMapper,Transpose[{dfts}]]]];
	AppendTo[results,"averagedFiles"->fn];
	
	For[i=1,i<=Length[iri],i++,
		AppendTo[iri[[i]],irResults[[i]]];
	];
	
	AppendTo[results,ProcessElectronPolarizationFromFourierCoefficients[avgFourierComponents]];
	AppendTo[results,avgFourierComponents];
	AppendTo[results,<|"IndividualRunInfo"->iri|>];
	<|fn[[1]]->results|>
];

Clear[ProcessElectronPolarizationFileAverageSTOKES];
ProcessElectronPolarizationFileAverageSTOKES[fn_List,darkCountsSubtract_:True,cn_:True]:=Module[
	{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,dftResults,dfts,
	fcTitles,res,errors,avgFourierComponents,
	d,sin,cos,c0,c2,c4,s2,s4,stokes,
	stokes2,stokesStdDev,avg,
	transposedValues,
	irResults (*individual run results*),iri(* invidual run info *)},
	
	results=<||>;
	irResults=<||>;
	fcTitles={"Sin Coefficients","Cos Coefficients"};

	d=ProcessElectronPolarizationFileAverageFOURIER[fn,darkCountsSubtract,cn];
	
	iri=Normal[Dataset[d][Values,"IndividualRunInfo"][[1]]];
	dfts=Normal[Dataset[d][Values,"IndividualRunInfo",Values,"fcCoefficients"][[1]]];

	sin=Through[dfts["Sin Coefficients"]];
	cos=Through[dfts["Cos Coefficients"]];
	c0=Transpose[Values[cos]][[0+1]];
	c2=Transpose[Values[cos]][[2+1]];
	c4=Transpose[Values[cos]][[4+1]];
	s2=Transpose[Values[sin]][[2+1]];
	s4=Transpose[Values[sin]][[4+1]];
	stokes=Apply[CalculateStokesFromFourierCoefficients,Transpose[{c0,c2,c4,s2,s4}],{1}];

	irResultsMapper[results_]:=AssociationThread[{"stokes"},results];
	fcThreader[values_]:=AssociationThread[Keys[dfts[[1,1]]],values];

	AppendTo[irResults,AssociationThread[fn,Map[irResultsMapper,Transpose[{stokes}]]]];

	avg=AverageStokes[stokes];

	For[i=1,i<=Length[iri],i++,
		AppendTo[iri[[i]],irResults[[i]]];
	];
	
	PrependTo[results,<|"avgStokes"->avg|>];
	AppendTo[results,<|"IndividualRunInfo"->iri|>];
	PrependTo[results,<|ProcessElectronPolarizationFromStokesParameters[avg]|>];
	<|fn[[1]]->results|>
];




(* ::Subsubsection::Closed:: *)
(*Vector Averaging*)


AverageFourierCoefficients[fcList_Association]:=Module[
{},
fcList
];

AverageStokes[stokes_]:=Module[
{i,transposedValues,stokesAvg,stokesStdDevMean},
	stokesAvg=<||>;
	stokesStdDevMean=<||>; (*Just creating another object with the same length as stokes *)
	transposedValues=Transpose[Values[stokes]];
	For[i=2,i<=Length[stokes[[1]]],i++,
		AppendTo[stokesAvg,stokesNames[[i]]->Mean[transposedValues[[i]]*transposedValues[[1]]]];
		AppendTo[stokesStdDevMean,stokesErrNames[[i]]->StandardDeviation[transposedValues[[i]]*transposedValues[[1]]]/Sqrt[Length[transposedValues[[1]]]]];
	];
	PrependTo[stokesAvg,stokesNames[[1]]->Mean[transposedValues[[1]]]];
	PrependTo[stokesStdDevMean,stokesErrNames[[1]]->StandardDeviation[transposedValues[[1]]]/Sqrt[Length[transposedValues[[1]]]]];
	Join[Association[stokesAvg],Association[stokesStdDevMean]]
];

StokesAverage[stokesList_]:=Module[
{},
	NormalizeStokes[Total[Map[UnNormalizeStokes,stokesList]]]
];


StokesAdd[s1_,s2_]:=Module[
{},
	NormalizeStokes[Total[Map[UnNormalizeStokes,{s1,s2}]]]
];

Clear[NormalizeStokes];
NormalizeStokes[stokesVector_]:=Module[{intensity,newVector,i},
	newVector=stokesVector;
	intensity=stokesVector[[1]];
	For[i=2,i<=Length[stokesVector],i++,
		newVector[[i]]=stokesVector[[i]]/intensity;
	];
	newVector
];

UnNormalizeStokes[stokesVector_]:=Module[{intensity,newVector,i},
newVector=stokesVector;
intensity=stokesVector[[1]];
For[i=2,i<=Length[stokesVector],i++,
newVector[[i]]=stokesVector[[i]]*intensity;
];
newVector
];


(* ::Subsection::Closed:: *)
(*Background Subtraction Methods*)


GetAverageElectronPolarizationWithBackgroudSubtractionSUM[fn_,fnBack_,cn_:True]:=Module[
{rawSignal,beamBackgroundSignal,processedSignal,return=<||>},
	rawSignal=GetAverageCurrentNormalizedCountRateFromFilenames[fn];
	AppendTo[return,"signalCurrentNormalized"->rawSignal];
	If[cn==True,
		beamBackgroundSignal=GetAverageCurrentNormalizedCountRateFromFilenames[fnBack];,
		beamBackgroundSignal=GetAverageCountRateFromFilenames[fnBack];
	];
	AppendTo[return,"beamBackground"->beamBackgroundSignal];
	processedSignal=ProcessElectronPolarizationFromSignal[rawSignal["avgSignal"]-beamBackgroundSignal["avgSignal"]];
	AppendTo[return,"processedSignal"->processedSignal]
];

Clear[GetAverageElectronPolarizationWithBackgroudSubtractionELECTRON];
GetAverageElectronPolarizationWithBackgroudSubtractionELECTRON[fn_,fnBack_,cn_:True]:=Module[
{rawSignal,beamBackgroundSignal,processedSignal=<||>,return=<||>,individualRunInfo,
pe,peerr},
	individualRunInfo=Map[ProcessElectronPolarizationFileSubtractBackground[#,fnBack]&,fn];
	AppendTo[processedSignal,"signal"->Normal[Dataset[individualRunInfo][Mean,"signal"]]];
	pe=Dataset[individualRunInfo][Mean,"P_e"];
	peerr=Dataset[individualRunInfo][StandardDeviation[#]/Length[#]&,"P_e"];
	AppendTo[processedSignal,"P_e"->Around[pe,peerr]];
	AppendTo[processedSignal,"fc"->Normal[Dataset[individualRunInfo][Mean,"fc"]]];
	AppendTo[processedSignal,"stokes"->Normal[Dataset[individualRunInfo][Mean,"stokes"]]];
	AppendTo[processedSignal,"individualRunInfo"->individualRunInfo];
	AppendTo[return,"processedSignal"->processedSignal];
	return
];

SubtractElectronPolarizationBackgroundFOURIER[signal_,background_]:=Module[
{},
	signal
];

SubtractElectronPolarizationBackgroundSTOKES[signal_,background_]:=Module[
{},
	signal
];
SubtractElectronPolarizationBackgroundELECTRON[signal_,background_]:=Module[
{},
	signal
];


(* ::Chapter:: *)
(*DataRun Processing*)


(* Place all files from a run in a folder. Provide that folder as an argument to this function, 
and it will return the electron polarization for the no pump, S+ pump and S- pump.

noBeamBackground is an association of pumpTypes to fileNames for the noBeam Background. beamOnBackground  is the same *)

Clear[ProcessElectronPolarizationRun];
ProcessElectronPolarizationRun[folder_,noBeamBackground_,beamOnBackground_,pumpTypes_:{"noPump","s+Pump","s-Pump"}]:=
Module[
{fn,fnCategorized,
i,
pumpType,
return},
SetDirectory[folder];
return=<||>;
fn=FileNames[RegularExpression["POL.*_[0-9]*.dat"]];
fnCategorized=<||>;
For[i=1,i<=Length[pumpTypes],i++,
	pumpType=pumpTypes[[i]];
	AppendTo[fnCategorized,pumpType->Take[fn,{i,-1,Length[pumpTypes]}]];
	darkCountRate=GetAverageCountRateFromFilenames[noBeamBackground[pumpType],False]["avgSignal"];
	AppendTo[return,pumpType->GetAverageElectronPolarizationWithBackgroudSubtractionELECTRON[fnCategorized[pumpType],beamOnBackground[pumpType]]];
];
ResetDirectory[];
return
];

Clear[POLProcessElectronPolarizationNoBackground];
POLProcessElectronPolarizationNoBackground[folder_,pumpTypes_:{"noPump","s+Pump","s-Pump"}]:=
Module[
{fn,fnCategorized,
i,
pumpType,
return},
SetDirectory[folder];
return=<||>;
fn=FileNames[RegularExpression["POL.*_[0-9]*.dat"]];
fnCategorized=<||>;
For[i=1,i<=Length[pumpTypes],i++,
	pumpType=pumpTypes[[i]];
	AppendTo[fnCategorized,pumpType->Take[fn,{i,-1,Length[pumpTypes]}]];
	darkCountRate=100;
	AppendTo[return,pumpType->ProcessElectronPolarizationFileAverage[fnCategorized[pumpType],True,True]];
];
ResetDirectory[];
return
];


(* ::Chapter::Closed:: *)
(*Monitor Counts Functions*)


AnalyzeMonitorCountsFile[fileName_]:=Module[
{f,comments,val,sdm,percentVar,
lpData,plot,results=<||>},

f=ImportFile[fileName];
comments=f[[1]]["Comments"];
val=N[f[[2]][Mean,"Count"]];
sdm=f[[2]][SDM,"Count"];
percentVar=sdm/val;
AppendTo[results,<|"comments"->comments,"avgCount"->val,"sdmCount"->sdm,"percentVarCount"->percentVar|>];

val=N[f[[2]][Mean,"Current"]];
sdm=f[[2]][SDM,"Current"];
lpData=N[f[[2]][All,"Current"]];
percentVar=sdm/val;
plot=ListPlot[lpData,PlotRange->{.8*val,1.2*val}];
AppendTo[results,<|"avgCurrent"->val,"sdmCurrent"->sdm,"percentVarCurrent"->percentVar,"plotCurrent"->plot,"lpDataCurrent"->lpData|>];
results
];


(* ::Chapter::Closed:: *)
(*Operations On Completed Datasets*)


GeneratePolarizationTransferEfficiencyGraph[folder_]:=Module[
{runAnalysisFolder=folder,
dr,de,
pe,allValues,allErrors,
prb,slopes,densities,
a,sa,
data,dataNoError,
uncertaintyPRb,uncertaintyPe,
syf,af,saf
},
dr=Import[FileNameJoin[{runAnalysisFolder,"P_Rb.wdx"}]];
de=Import[FileNameJoin[{runAnalysisFolder,"P_e.wdx"}]];
pe=GetPeErrorFromDataset[de];
allValues=Map[#1[[1]]&,pe];
allErrors=Map[#1[[2]]&,pe];
(* This first method of calculating Subscript[P, Rb] uses the density calculated in each run 
in the calculating of Subscript[P, Rb]*)
prb=Dataset[dr][All,"P_Rb_2",{"S+ Pump","S- Pump"},"P_Rb"]//Values//Normal//Flatten;
(* This second method of calculating Subscript[P, Rb] uses the average of all of the Subscript[n, Rb]
measurements in the set to use in the calculation of Subscript[P, Rb]. *)
slopes=Dataset[dr][All,"P_Rb_2",{"S+ Pump","S- Pump"},"Model Fit"]//Values//Normal;
densities=Dataset[dr][All,"n_Rb","n_Rbx10^12"]//Normal//Mean;
prb=CalculateRubidiumPolarizationOnePump[slopes,Mean[densities]*1*^12];
a=Total[(prb*allValues)/(allErrors)^2]/Total[prb^2/allErrors^2];
sa=1/Sqrt[Total[prb^2/allErrors^2]];
data=Transpose[{prb,pe}];
dataNoError=Apply[{#1["Value"],#2["Value"]}&,data,{1}];
uncertaintyPRb=Apply[#1["Uncertainty"]&,data,{1}];
uncertaintyPe=Apply[#2["Uncertainty"]&,data,{1}];
syf=Sqrt[(a["Value"]*uncertaintyPRb)^2+uncertaintyPe^2];
af=Total[(prb*allValues)/(syf)^2]/Total[prb^2/syf^2];
saf=1/Sqrt[Total[prb^2/syf^2]];
];

GetPeErrorFromDataset[dataset_,beamBackground_]:=Module[
{sig,peValSP,peErrSP,peValSM,peErrSM,ret,
allValues,allErrors,
plusBack,minusBack},
(* First let's get those uncertainties in P_e *)
plusBack=GetAverageCurrentNormalizedCountRateFromFilenames[beamBackground["s+Pump"]]["avgSignal"];
minusBack=GetAverageCurrentNormalizedCountRateFromFilenames[beamBackground["s-Pump"]]["avgSignal"];
sig=Dataset[dataset][All,{"s+Pump"},"signalCurrentNormalized","IndividualRunInfo",Values,ProcessElectronPolarizationFromSignal[#["signal"]-plusBack]["P_e"]&][All]//Normal;
peValSP=Apply[Mean,sig//Values,{1}];
peErrSP=Apply[StandardDeviation[#]/Sqrt[Length[#]]&,sig//Values,{1}];
(* First let's get those uncertainties in P_e *)
sig=Dataset[dataset][All,{"s-Pump"},"signalCurrentNormalized","IndividualRunInfo",Values,ProcessElectronPolarizationFromSignal[#["signal"]-minusBack]["P_e"]&][All]//Normal;
peValSM=Apply[Mean,sig//Values,{1}];
peErrSM=Apply[StandardDeviation[#]/Sqrt[Length[#]]&,sig//Values,{1}];
allValues=Riffle[peValSP,peValSM];
allErrors=Riffle[peErrSP,peErrSM];
Apply[Around[#1,#2]&,Transpose[{allValues,allErrors}],{1}]
];


(* ::Chapter::Closed:: *)
(*Optics Polarization*)


AnalyzeRetroReflection[signalfn_,backgroundfn_]:=Module[
{fn,raw,data,threader,dSig,signal,dBack,background,bss,plotLabel,nlm,dataPlot,fitPlot,fullPlot,
A,B,F,G,offset,roughAmplitude},
fn=signalfn;
raw=Import[fn];
data=raw[[21;;]];
threader[line_]:=AssociationThread[{"th","cts","curr","curr_err","angle"},line];
dSig=Dataset[Map[threader,data]];
signal=dSig[All,{"angle","curr"}][Values]//Normal;
fn=backgroundfn;
raw=Import[fn];
data=raw[[21;;]];
dBack=Dataset[Map[threader,data]];
background=dBack[All,{"angle","curr"}][Values]//Normal;
If[Length[background]!= Length[signal],
background=Mean[background[[All,2]]],
background=background[[All,2]]];
bss=signal;
bss[[All,2]]=signal[[All,2]]-background;
bss=Abs[bss];

plotLabel="Retro Reflection Fit:\n" <> FileNameTake[signalfn];
roughAmplitude=Differences[MinMax[bss[[All,2]]]][[1]];
nlm=NonlinearModelFit[bss,{A+B*Cos[4(F-G)],Pi/4>G>-Pi/4},{A,B,G},F];
dataPlot=ListPlot[{bss},PlotLabel->plotLabel,FrameLabel->{"Polarimeter Position (rad)","Intensity (arb units)"}];
fitPlot=Plot[A+B*Cos[4(F-G)]/.nlm["BestFitParameters"],{F,0,2Pi}];
fullPlot=Show[{dataPlot,fitPlot}];
offset=(G/.nlm["BestFitParameters"])/\[Pi]*180;
If[offset>45, offset=offset-90];
<|"model"->nlm,"plot"->fullPlot,"offset"->Around[offset,nlm["ParameterErrors"][[3]]]|>
]
