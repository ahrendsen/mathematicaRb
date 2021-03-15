(* ::Package:: *)

(* ::Title:: *)
(*Electron Polarization Functions*)


(* ::Chapter:: *)
(*Import other needed functions*)


Needs["ErrorBarPlots`"]
SetDirectory[$PACKAGES];
Import["fileManipulation.wl"];
Import["dataManipulation.wl"];
Import["apparatus.wl"];
Import["physicalConstants.wl"];
ResetDirectory[];


(* ::Chapter:: *)
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
bc=BarChart[sinValuesNoZero,ChartStyle->colorBlindPalette[[2]],ChartLabels->bcLabels,ImageSize->{720,450},Frame->True];
AppendTo[results,"SinCoefficients"->bc];

bcLabels=Drop[Keys[fc[[2]]],{2,-1,2}];
cosValues=Values[fc[[2]]]/Values[fc[[2]]][[1]];
cosValuesNoZero=Drop[cosValues,{1}];
bcLabelsNoZero=Drop[bcLabels,{1}];
blanks=ConstantArray["",Length[bcLabelsNoZero]];
bcLabelsSparse=Riffle[blanks,bcLabelsNoZero];
bc=BarChart[cosValuesNoZero,ChartStyle->colorBlindPalette[[3]],ChartLabels->bcLabelsSparse,ImageSize->{720,450},Frame->True];
AppendTo[results,"CosCoefficients"->bc];
results
];


(* ::Subchapter:: *)
(*Raw Data*)


(* ::Subchapter:: *)
(*Stokes*)


(* ::Chapter:: *)
(*Polarization Functions*)


(* ::Subchapter:: *)
(*Constants*)


(* ::Input::Initialization:: *)
stokesNames={"p0","p1","p2","p3_mag","p3_c2","p3_s2"};
stokesErrNames={"p0err","p1err","p2err","p3_magerr","p3_c2err","p3_s2err"};
nA=1*^-9;
mTorr=1*^-3;
alpha=20.4*\[Pi]/180;
beta=68.7*\[Pi]/180;
deltaAB=alpha-beta;
delta=94.54*\[Pi]/180 (*Munir's reported 1.66 plus or minus .01*);
darkCountRate=14;


alphaOld=3.2*\[Pi]/180;
betaOld=-6.95*\[Pi]/180;
rev=1;


GenerateListOfFileNames[prefix_,date_,runtimes_]:=Module[{j,fileNames},
fileNames={};
For[j=1,j<=Length[runtimes],j++,
AppendTo[fileNames,StringJoin[prefix,date,"_",ToString[runtimes[[j]]],".dat"]]
];
fileNames
];
(* If given a list of files*)
GetIndexOfFilenameFromTimestamp[files_,timestamp_]:=Module[{i,index},
For[i=1,i<=Length[files],i++,
If[StringContainsQ[ files[[i]],timestamp],index=i]
];
index
];

FourierFit[list_,rev_]:=Module[{fit,function,dataPts,stepSize,dataPtsPerRev,angles,data(*,C0,C2,C4,S2,S4*)},

function[\[Theta]_]= C0+C2 Cos[2(\[Theta]+\[Theta]o)]+C4 Cos[4(\[Theta]+\[Theta]o)]+S2 Sin[2(\[Theta]+\[Theta]o)]+S4 Sin[4(\[Theta]+\[Theta]o)];
dataPts=Length[list];
dataPtsPerRev=dataPts/rev;
angles=Range[0,(2\[Pi]*rev)-2\[Pi]/dataPtsPerRev,2\[Pi]/dataPtsPerRev];
data=Transpose[{list,angles}];
fit=NonlinearModelFit[data,function[\[Theta]],{{C0,1000},C2,C4,S2,S4,\[Theta]o},\[Theta]]
];

StokesParametersFromFourierCoefficients[fc_(*Fourrier Coefficients*),rev_:rev,alpha_:alpha,beta0_:beta,delta_:delta]:=Module[{c0,c2,c4,s2,s4,stokes},
stokes={0,0,0,0,0,0};
(*c0=fc[[2]][[1]][[1]];
c2=fc[[2]][[1]][[3]];
c4=fc[[2]][[1]][[5]];
s2=fc[[1]][[1]][[3]];
s4=fc[[1]][[1]][[5]];
*)
c0=fc["Cos Coefficients"][0];
c2=fc["Cos Coefficients"][2];
c4=fc["Cos Coefficients"][4];
s2=fc["Sin Coefficients"][2];
s4=fc["Sin Coefficients"][4];

stokes[[1]]=c0-(1+Cos[delta])/(1-Cos[delta])*(c4*Cos[4 alpha+4 beta0]+s4*Sin[4 alpha+4 beta0]);
stokes[[2]]=2/(1-Cos[delta])*(c4*Cos[2alpha + 4 beta0]+s4*Sin[2alpha + 4 beta0])/stokes[[1]];
stokes[[3]]=2/(1-Cos[delta])*(s4*Cos[2alpha + 4 beta0]-c4*Sin[2alpha + 4 beta0])/stokes[[1]];
stokes[[4]]=Sqrt[c2^2+s2^2]/(Sin[delta]^2 stokes[[1]]);
stokes[[5]]=c2/(Sin[delta]Sin[2 alpha + 2 beta0]stokes[[1]]);
stokes[[6]]=-s2/(Sin[delta]Cos[2 alpha + 2beta0]stokes[[1]]);
stokes];

(* 
Someday when I'm bored, I can try to implement the asoociation version of the stokes vectors. Today is not that day
StokesParametersFromFourierCoefficients[fc_(*Fourrier Coefficients*),rev_:rev,alpha_:alpha,beta0_:beta,delta_:delta]:=Module[{c0,c2,c4,s2,s4,stokes},
stokes=<||>;
c0=fc[[2]][[1]][[1]];
c2=fc[[2]][[1]][[3]];
c4=fc[[2]][[1]][[5]];
s2=fc[[1]][[1]][[3]];
s4=fc[[1]][[1]][[5]];
AppendTo[stokes,stokesNames[[1]]->c0-(1+Cos[delta])/(1-Cos[delta])*(c4*Cos[4 alpha+4 beta0]+s4*Sin[4 alpha+4 beta0])];
AppendTo[stokes,stokesNames[[2]]->2/(1-Cos[delta])*(c4*Cos[2alpha + 4 beta0]+s4*Sin[2alpha + 4 beta0])/stokes[[1]]];
AppendTo[stokes,stokesNames[[3]]->2/(1-Cos[delta])*(s4*Cos[2alpha + 4 beta0]-c4*Sin[2alpha + 4 beta0])/stokes[[1]]];
AppendTo[stokes,stokesNames[[4]]->Sqrt[c2^2+s2^2]/(Sin[delta]^2stokes[[1]])];
AppendTo[stokes,stokesNames[[5]]->c2/(Sin[delta]Sin[2 alpha + 2 beta0]stokes[[1]])];
AppendTo[stokes,stokesNames[[6]]->-s2/(Sin[delta]Cos[2 alpha + 2beta0]stokes[[1]])];
stokes
];
*)

StokesParametersFromRawFourierData[intensityArray_,rev_:rev,alpha_:alpha,beta0_:beta,delta_:delta]:=
Module[{fc(*Fourrier Coefficients*),c0,c2,c4,s2,s4,stokes},
fc=DFT[intensityArray,rev];
StokesParametersFromFourierCoefficients[fc,rev,alpha,beta0,delta]
];

StokesParametersFromPolFile[polFile_,alpha_:alpha,beta0_:beta,delta_:delta,rev_:1]:=
Module[{intensityArray,f},
intensityArray=Normal[polFile[[2]][All,"COUNT"]];
StokesParametersFromRawFourierData[intensityArray,rev,alpha,beta0,delta]
];

StokesParametersFromPolFileName[polFileName_,alpha_:alpha,beta0_:beta,delta_:delta]:=
Module[{intensityArray,f,rev},
f=ImportFile[polFileName];
intensityArray=Normal[f[[2]][All,"COUNT"]];
rev=f[[1]]["REV"];
StokesParametersFromPolFile[f,alpha,beta0,delta,rev_]
];

GetLinearPolarizationFraction[fileName_]:=Module[{data,angle,intensityArray,angleRAD,fc,intensity,linearPolarization,linearPolarizationFraction},
data=GetFileDataset[fileName];
intensityArray=Transpose[data][[2]];
fc=DFT[intensityArray,1] (*fc=fourier coefficients*);
intensity=fc[[2]][[1]];
linearPolarization=Sqrt[fc[[2]][[3]]^2+fc[[1]][[3]]^2];
linearPolarizationFraction=linearPolarization/intensity
];
SetAttributes[GetLinearPolarizationFraction,Listable];
GetAngle[fileName_]:=Module[{data,intensityArray,angleRAD,fc,intensity,linearPolarization,linearPolarizationFraction},
data=GetFileDataset[fileName];
intensityArray=Transpose[data][[2]];
fc=DFT[intensityArray,1] (*fc=fourier coefficients*);
.5*ArcTan[fc[[1]][[3]],fc[[2]][[3]]]
];
SetAttributes[GetAngle,Listable];



GetAverageCountRateFromFileNames[fileNames_]:=Module[{files},
files=ImportFile[fileNames];
GetAverageCountRateFromFiles[files]
];

(* See Nate Clayburn's explanation on Background Subtraction p. 129 eq. 93*)
(* Returns a list of two items.
	the first item is a list of the intensity reaching the PMT at each step of the polarimeter.
	the second item is the standard deviation of these counts
*)
GetAverageCountRateFromFiles[files_]:=Module[{i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j},
nDataPts=Length[Normal[files[[1]][[2]][All,"COUNT"]]];
sum={};errorSum={};
For[j=1,j<= nDataPts,j++,
sum=AppendTo[sum,0];
errorSum=AppendTo[errorSum,0];
];
For[i=1,i<=Length[files],i++,
counts=Normal[files[[i]][[2]][All,"COUNT"]];
current=Normal[files[[i]][[2]][All,"CURRENT"]];
dwell=files[[i]][[1]]["DWELL(s)"];
nCycles=Length[files];
For[j=1,j<=nDataPts,j++,
sum[[j]]+=counts[[j]]/nCycles/dwell;
errorSum[[j]]+=counts[[j]]/(nCycles^2*dwell^2);
];
];
For[j=1,j<= nDataPts,j++,
errorSum[[j]]=Sqrt[errorSum[[j]]];
];
{sum,errorSum}
];

GetCurrentNormalizedAverageCountRateFromFiles[files_]:=Module[{t,scale,darkCountsRate,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j,rev},
nDataPts=files[[1]][[1]]["DATAPPR"];
rev=files[[1]][[1]]["REV"];
scale=files[[1]][[1]]["SCALE"];
darkCountsRate=10;
sum={};errorSum={};
For[j=1,j<= nDataPts,j++,
sum=AppendTo[sum,0];
errorSum=AppendTo[errorSum,0];
];
For[i=1,i<=Length[files],i++,
counts=Normal[files[[i]][[2]][All,"COUNT"]];
current=Normal[files[[i]][[2]][All,"CURRENT"]];
dwell=files[[i]][[1]]["DWELL(s)"];
nCycles=Length[files]*rev;
t=dwell*nCycles;
For[j=1,j<=nDataPts,j++,
sum[[j]]+=(counts[[j]]-darkCountsRate*dwell)/t/Abs[current[[j]]*10^(-scale)/(1*^-9)];
errorSum[[j]]+=(counts[[j]]-darkCountsRate*dwell)/(t^2*(current[[j]]*10^(-scale)/(1*^-9))^2);
];
];
For[j=1,j<= nDataPts,j++,
errorSum[[j]]=Sqrt[errorSum[[j]]];
];
{sum,errorSum}
];

GetCurrentNormalizedAverageCountRateFromFileNames[fileNames_]:=Module[{filesPass,darkCountsRate,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j},
filesPass=ImportFile[fileNames];
GetCurrentNormalizedAverageCountRateFromFiles[filesPass]
];

GetCurrentNormalizedMolyCountRateFromFileNames[fileNames_,darkCountsRate_]:=Module[{files},
files=ImportFile[fileNames];
GetCurrentNormalizedMolyCountRateFromFiles[files,darkCountsRate]
];

GetCurrentNormalizedMolyCountRateFromFiles[files_,darkCountsRate_]:=Module[{currentnA,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j,currentScale},
nDataPts=Length[darkCountsRate[[1]]];
sum={};errorSum={};
For[j=1,j<= nDataPts,j++,
sum=AppendTo[sum,0];
errorSum=AppendTo[errorSum,0];
];
For[i=1,i<=Length[files],i++,
currentScale=files[[i]][[1]]["SCALE"];
counts=Normal[files[[i]][[2]][All,"COUNT"]];
current=Normal[files[[i]][[2]][All,"CURRENT"]];
currentnA=current*10^(-currentScale)/nA;
dwell=files[[i]][[1]]["DWELL(s)"];
nCycles=Length[files];
For[j=1,j<=nDataPts,j++,
sum[[j]]+=(counts[[j]]/dwell-darkCountsRate[[1]][[j]])/(nCycles*Abs[currentnA[[j]]]);
errorSum[[j]]+=counts[[j]]/(nCycles^2*dwell^2*currentnA[[j]]^2)-darkCountsRate[[2]][[j]]/(nCycles^2*currentnA[[j]]^2);
];
];
For[j=1,j<= nDataPts,j++,
errorSum[[j]]=Sqrt[errorSum[[j]]];
];
{sum,errorSum}
];

BackgroundRateCurrentNorm[files_,darkCountsFiles_]:=Module[{darkCountsRate,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j},
darkCountsRate=GetAverageCountRateFromFiles[darkCountsFiles];
nDataPts=Length[darkCountsRate[[1]]];
sum={};errorSum={};
For[j=1,j<= nDataPts,j++,
sum=AppendTo[sum,0];
errorSum=AppendTo[errorSum,0];
];
For[i=1,i<=Length[files],i++,
counts=Normal[files[[i]][[2]][All,"COUNT"]];
current=Normal[files[[i]][[2]][All,"CURRENT"]];
dwell=files[[i]][[1]]["DWELL(s)"];
nCycles=Length[files];
For[j=1,j<=nDataPts,j++,
sum[[j]]+=counts[[j]]/nCycles/dwell/Abs[current[[j]]]-darkCountsRate[[1]][[j]]/nCycles/Abs[current[[j]]];
errorSum[[j]]+=counts[[j]]/(nCycles^2*dwell^2*current[[j]]^2)-darkCountsRate[[2]][[j]]/(nCycles^2*current[[j]]^2);
];
];
For[j=1,j<= nDataPts,j++,
errorSum[[j]]=Sqrt[errorSum[[j]]];
];
{sum,errorSum}
];

AverageSignalRateCurrentNorm[files_,molyCountsRate_,darkCountsFiles_]:=Module[{darkCountsRate,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,rev,j},
darkCountsRate=GetAverageCountRateFromFiles[darkCountsFiles];
nDataPts=Length[darkCountsRate[[1]]];
sum={};errorSum={};
For[j=1,j<= nDataPts,j++,
sum=AppendTo[sum,0];
errorSum=AppendTo[errorSum,0];
];
For[i=1,i<=Length[files],i++,
counts=Normal[files[[i]][[2]][All,"COUNT"]];
current=Normal[files[[i]][[2]][All,"CURRENT"]];
dwell=files[[i]][[1]]["DWELL(s)"];
rev=files[[i]][[1]]["REV"];
nCycles=Length[files];
For[j=1,j<=nDataPts,j++,
sum[[j]]+=counts[[j]]/nCycles/rev/dwell/Abs[current[[j]]]-darkCountsRate[[1]][[j]]/rev/nCycles/Abs[current[[j]]]-molyCountsRate[[1]][[j]]/rev/nCycles;
errorSum[[j]]+=counts[[j]]/(nCycles^2*dwell^2*rev^2*current[[j]]^2)-darkCountsRate[[2]][[j]]/(nCycles^2*rev^2*current[[j]]^2)-molyCountsRate[[2]][[j]]/(nCycles^2*rev^2);
];
];
For[j=1,j<= nDataPts,j++,
errorSum[[j]]=Sqrt[errorSum[[j]]];
];
{sum,errorSum}
];

SignalRateCurrentNormFromFiles[file_,molyCountFiles_,darkCountsFiles_]:=Module[{molyCountsRate,darkCountsRate,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j},
darkCountsRate=AverageGetAverageCountRateFromFiles[darkCountsFiles];
molyCountsRate=AverageBackgroundRateCurrentNorm[molyCountFiles,darkCountsFiles];
SignalRateCurrentNormFromRates[file,molyCountFiles,darkCountsFiles]
];
SignalRateCurrentNormFromRates::usage ="SignalRateCurrentNormFromRates takes three arguments. First, the file you wish to background subtract, second, the fourier data for the current dependent background, third, the fourier data for the current independent background."
SignalRateCurrentNormFromRates[fileName_,molyCountsRate_,darkCountsRate_]:=Module[{file,i,sum,counts,dwell,current,nCycles,errorSum, error,nDataPts,j,mTorrHe},
file=ImportFile[fileName];
nDataPts=Length[Normal[file[[2]][All,"COUNT"]]];
sum={};errorSum={};
For[j=1,j<= nDataPts,j++,
sum=AppendTo[sum,0];
errorSum=AppendTo[errorSum,0];
];
counts=Normal[file[[2]][All,"COUNT"]];
current=Normal[file[[2]][All,"CURRENT"]];
dwell=file[[1]]["DWELL(s)"];
nCycles=file[[1]]["REV"];
mTorrHe=file[[1]]["CVGauge(He)(Torr)"]/mTorr;
For[j=1,j<=nDataPts,j++,
sum[[j]]+=((counts[[j]]/dwell-darkCountsRate[[1]][[j]])/Abs[current[[j]]]-molyCountsRate[[1]][[j]])/(nCycles*mTorrHe);
errorSum[[j]]+=counts[[j]]/(nCycles^2*dwell^2*current[[j]]^2*mTorrHe^2)-darkCountsRate[[2]][[j]]^2/(nCycles^2*current[[j]]^2*mTorrHe^2)-molyCountsRate[[2]][[j]]^2/(nCycles^2*mTorrHe^2);
];
For[j=1,j<= nDataPts,j++,
errorSum[[j]]=Sqrt[errorSum[[j]]];
];
{sum,errorSum}
];

CalculateOneFileBackgroundSubtractedStokes[polFileName_,molyCountRate_,darkCountRate_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{signal,stokes,file,rev},
file=ImportFile[polFileName];
rev=file[[1]]["REV"];
signal=SignalRateCurrentNormFromRates[polFileName,molyCountRate,darkCountRate];
StokesParametersFromRawFourierData[signal[[1]],rev,alpha,beta,delta]
];

CalculateOneFileStokes[polFileName_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{signal,stokes,file,rev},
file=ImportFile[polFileName];
rev=file[[1]]["REV"];
signal=GetCurrentNormalizedAverageCountRateFromFiles[{file}];
Sow[signal,"sig"];
StokesParametersFromRawFourierData[signal[[1]],rev,alpha,beta,delta]
];

CalculateNoBeamOneFileStokes[polFileName_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{signal,stokes,file,rev},
file=ImportFile[polFileName];
rev=file[[1]]["REV"];
signal=GetAverageCountRateFromFiles[{file}];
StokesParametersFromRawFourierData[signal[[1]],rev,alpha,beta,delta]
];

CalculateOneFileCurrent[polFileName_]:=Module[{signal,stokes,file,rev,currentScale,currentValues,currentPlot,currentAvg,currentStd},
file=ImportFile[polFileName];
rev=file[[1]]["REV"];
currentScale=file[[1]]["SCALE"];
currentValues=Normal[file[[2]][All,"CURRENT",Abs]]*10^(-currentScale)/nA;
currentAvg=Mean[currentValues];
currentStd=StandardDeviation[currentValues];
currentPlot=ListPlot[currentValues/currentAvg,PlotLabel->"Current, Normalized to Average Current",FrameLabel->{"Stepper Motor Position","Current (nA)"},PlotRange->{.8,1.2}];
currentPlot=ListPlot[currentValues,PlotRange->{.8*currentAvg,1.2*currentAvg}];
<|"currentPlot"->currentPlot,"currentAvg"->currentAvg,"currentStd"->currentStd|>
];


CalculateAverageBackgroundSubtractedStokesFromFilesAndCountRates[polFileNames_,molyCountRate_,darkCountRate_,alpha_:alpha,beta_:beta,delta_:delta]:=Module[{signal,stokes,stokesStdDev,file,rev,j,stokesValues,transposedValues},
stokes=<|"p0"->0,"p1"->0,"p2"->0,"p3_mag"->0,"p3_c2"->0,"p3_s2"->0|>;
stokesStdDev={0,0,0,0,0,0};
stokesValues={};
For[j=1,j<=Length[polFileNames],j++,
stokes=CalculateOneFileStokes[polFileNames[[j]],molyCountRate,darkCountRate,alpha,beta,delta];
AppendTo[stokesValues,stokes]
];

Print[AverageStokes[stokesValues]];
AverageStokes[stokesValues]
]

CalculateElectronPolarization[stokes_Association]:=Module[{p1,p3},
p3=stokes["p3_s2"];
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
p3Values=Transpose[stokes][[6]]; (*The sixth item is p3_s2 *)
peError=Sqrt[pe^2*((averageStokes["p3_s2err"]/averageStokes["p3_s2"])^2+((averageStokes["p1err"]*p1Coeff)/(p1Const+p1Coeff*averageStokes["p1"]))^2-Covariance[p1Values*p1Coeff+p1Const,p3Values*p3Coeff]^2/(p3Coeff*averageStokes["p3_s2"]*(p1Const+p1Coeff*averageStokes["p1"])))]
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

SubtractStokes[addedStokes_,subtractedStokes_]:=Module[
{unNormAdd,unNormSub},
unNormAdd=UnNormalizeStokes[addedStokes];
unNormSub=UnNormalizeStokes[subtractedStokes];
NormalizeStokes[unNormAdd-unNormSub]
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


QPOLProcessFileDataNoBackground[fileName_]:=Module[
{f,avg,stdDev,
t,p,
normCounts,
p3,pe,peErr,results},
results=<||>;
f=ImportFile[fileName];
t=ToExpression[f[[1]]["DWELL(s)"]];
p=ToExpression[f[[1]]["CVGauge(He)(Torr)"]];
normCounts=f[[2]][All,<|"COUNT+45"->#["COUNT+45"]/#["CURRENT+45"]/(t*p),"COUNT-45"->#["COUNT-45"]/#["CURRENT-45"]/(t*p)|>&];
p3=Normal[Flatten[normCounts[All,{(#["COUNT+45"]-#["COUNT-45"])/(#["COUNT+45"]+#["COUNT-45"])}&]]];
pe=CalculateElectronPolarization[polConstxxunPolP1,Mean[p3]];
peErr=CalculateElectronPolarization[polConstxxunPolP1,StandardDeviation[p3]/Sqrt[Length[p3]]];
AppendTo[results,<|"fileName"->fileName,"P_e"->pe,"P_eerr"->peErr,"rawP3"->p3,"averageP3"->Mean[p3],"stdDevP3"->StandardDeviation[p3],"p1"->polConstxxunPolP1|>];
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
p=ToExpression[f[[1]]["CVGauge(He)(Torr)"]];
bnc=f[[2]][All,<|"COUNT+45"->(#["COUNT+45"]/t-avgNoBeam["COUNT+45"])/(#["CURRENT+45"]*p),"COUNT-45"->(#["COUNT-45"]/t-avgNoBeam["COUNT-45"])/(#["CURRENT-45"]*p)|>&];
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
p=ToExpression[btf[[1]]["CVGauge(He)(Torr)"]];
avgBT=btf[[2]][Mean,<|"COUNT+45"->(#["COUNT+45"]/t-avgNoBeam["COUNT+45"])/(#["CURRENT+45"]*p),"COUNT-45"->(#["COUNT-45"]/t-avgNoBeam["COUNT-45"])/(#["CURRENT-45"]*p)|>&];
f=ImportFile[fileName];
t=ToExpression[f[[1]]["DWELL(s)"]];
p=ToExpression[f[[1]]["CVGauge(He)(Torr)"]];
bnc=f[[2]][All,<|"COUNT+45"->(#["COUNT+45"]/t-avgNoBeam["COUNT+45"])/(#["CURRENT+45"]*p)-avgBT["COUNT+45"],"COUNT-45"->(#["COUNT-45"]/t-avgNoBeam["COUNT-45"])/(#["CURRENT-45"]*p)-avgBT["COUNT+45"]|>&];
avg=bnc[Mean,{"COUNT+45","COUNT-45"}];
p3=Normal[Flatten[bnc[All,{(#["COUNT+45"]-#["COUNT-45"])/(#["COUNT+45"]+#["COUNT-45"])}&]]];
pe=CalculateElectronPolarization[polConstxxunPolP1,Mean[p3]];
peErr=CalculateElectronPolarization[polConstxxunPolP1,StandardDeviation[p3]/Sqrt[Length[p3]]];
AppendTo[results,<|"P_e"->pe,"P_eerr"->peErr,"rawP3"->p3,"averageP3"->Mean[p3],"stdDevP3"->StandardDeviation[p3],"p1"->polConstxxunPolP1|>];
results
]


(* ::Subchapter:: *)
(*Calculate Polarization*)


(* ::Section:: *)
(*Single File*)


(* ::Subsection:: *)
(*Public Functions *)


(* The functions that will be regularly called from the Data Processing notebooks. *)
ProcessElectronPolarizationFile[fn_String]:=Module[
{counts,current,dwellTime,header,signal,
f,dataset,results},
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
AppendTo[results,ProcessElectronPolarizationFromSignal[signal,GetCurrentNormalizedCountRateError[counts,current,dwellTime]]](* See Private Functions for the rest of processing. *)
];

ProcessElectronPolarizationNoBeamBackgroundFile[fn_String]:=Module[
{counts,current,dwellTime,header,signal,
f,dataset,results},
results=<||>;
f=ImportFile[fn];
header=f[[1]];
AppendTo[results,header];
dataset=f[[2]];
AppendTo[results,"dataset"->Normal[dataset]];

counts=Normal[dataset[All,"COUNT"]];
dwellTime=header["DWELL(s)"];
signal=NormalizeNoBeamPolarizationData[counts,dwellTime];
AppendTo[results,"signal"->signal];
AppendTo[results,ProcessElectronPolarizationFromSignal[signal,GetCountRateError[counts,dwellTime]]](* See Private Functions for the rest of processing. *)
];

ProcessElectronPolarizationFileSubtractBackground[fn_String,fnBack_String]:=Module[
{counts,current,dwellTime,header,signal,
f,dataset,results,background},
results=<||>;
f=ImportFile[fn];
header=f[[1]];
AppendTo[results,header];
dataset=f[[2]];
AppendTo[results,"dataset"->Normal[dataset]];

counts=Normal[dataset[All,"COUNT"]];
current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
dwellTime=header["DWELL(s)"];
signal=NormalizePolarizationData[counts,current,dwellTime];
AppendTo[results,"signal"->signal];


f=ImportFile[fnBack];
header=f[[1]];
AppendTo[results,"backgroundHeader"->header];
dataset=f[[2]];
AppendTo[results,"datasetBack"->Normal[dataset]];

counts=Normal[dataset[All,"COUNT"]];
dwellTime=header["DWELL(s)"];
background=NormalizePolarizationData[counts,current,dwellTime];
AppendTo[results,"background"->background];

AppendTo[results,ProcessElectronPolarizationFromSignal[signal-background,Sqrt[signal+background]]](* See Private Functions for the rest of processing. *)
];


(* ::Subsection:: *)
(*Private Functions*)


(* The "helper functions" that the public functions will access to make their job easier. Won't typically be called from the notebook.*)

GetCurrentNormalizedCountRate[counts_,current_,dwellTime_]:=Module[
{countRate,r},
countRate=(counts/dwellTime)-darkCountRate;
r=countRate/current; (* r= Count rate *)
Normal[r]
];

GetCurrentNormalizedCountRateError[counts_,current_,dwellTime_]:=Module[
{countRateError,r},
countRateError=Sqrt[(Sqrt[counts]/dwellTime)^2+darkCountRate];
r=countRateError/current; (* r= Count rate *)
Normal[r]/Sqrt[Length[r]]
];

GetCountRate[counts_,dwellTime_]:=Module[
{countRate,r},
countRate=(counts/dwellTime)-darkCountRate;
r=countRate; (* r= Count rate *)
N[Normal[r]]
];

GetCountRateError[counts_,dwellTime_]:=Module[
{countRateError,r},
countRateError=Sqrt[(Sqrt[counts]/dwellTime)^2+darkCountRate];
r=countRateError; (* r= Count rate *)
Normal[r]/Sqrt[Length[r]]
];

CalculateStokesFromFourierCoefficients[c0_,c2_,c4_,s2_,s4_,rev_:rev,alpha_:alpha,beta0_:beta,delta_:delta]:=Module[{stokes},

stokes=<||>;
AppendTo[stokes,stokesNames[[1]]->c0-(1+Cos[delta])/(1-Cos[delta])*(c4*Cos[4 alpha+4 beta0]+s4*Sin[4 alpha+4 beta0])];
AppendTo[stokes,stokesNames[[2]]->2/(1-Cos[delta])*(c4*Cos[2alpha + 4 beta0]+s4*Sin[2alpha + 4 beta0])/stokes[[1]]];
AppendTo[stokes,stokesNames[[3]]->2/(1-Cos[delta])*(s4*Cos[2alpha + 4 beta0]-c4*Sin[2alpha + 4 beta0])/stokes[[1]]];
AppendTo[stokes,stokesNames[[4]]->Sqrt[c2^2+s2^2]/(Sin[delta]^2stokes[[1]])];
AppendTo[stokes,stokesNames[[5]]->c2/(Sin[delta]Sin[2 alpha + 2 beta0]stokes[[1]])];
AppendTo[stokes,stokesNames[[6]]->-s2/(Sin[delta]Cos[2 alpha + 2beta0]stokes[[1]])];
stokes

(*
stokes={0,0,0,0,0,0};

stokes[[1]]=c0-(1+Cos[delta])/(1-Cos[delta])*(c4*Cos[4 alpha+4 beta0]+s4*Sin[4 alpha+4 beta0]);
stokes[[2]]=2/(1-Cos[delta])*(c4*Cos[2alpha + 4 beta0]+s4*Sin[2alpha + 4 beta0])/stokes[[1]];
stokes[[3]]=2/(1-Cos[delta])*(s4*Cos[2alpha + 4 beta0]-c4*Sin[2alpha + 4 beta0])/stokes[[1]];
stokes[[4]]=Sqrt[c2^2+s2^2]/(Sin[delta]^2 stokes[[1]]);
stokes[[5]]=c2/(Sin[delta]Sin[2 alpha + 2 beta0]stokes[[1]]);
stokes[[6]]=-s2/(Sin[delta]Cos[2 alpha + 2beta0]stokes[[1]]);
stokes
*)
];

ProcessElectronPolarizationFromSignal[signal_List,signalErr_List:{}]:=Module[
{fc,s,c,plot,ePlot,fourierFit,fourierFitFunction,
stokes,electronPolarization,
results,err},
If[Length[signalErr]==0,err=ConstantArray[0,Length[signal]],err=signalErr];
results=<||>;
fc=DFT[signal];
s=fc["Sin Coefficients"];
c=fc["Cos Coefficients"];
ePlot=ListEPlot[Range[0,2\[Pi],2\[Pi]/60],signal,err,PlotLegends->{"Averaged Signal Data"}];
fourierFitFunction=c[0]+c[2]*Cos[2#]+s[2]*Sin[2#]+c[4]*Cos[4#]+s[4]*Sin[4#]&;
plot=Plot[fourierFitFunction[\[Theta]],{\[Theta],0,2\[Pi]},PlotRange->Full,PlotLegends->{"Fourier Coefficient Fit"}];
fourierFit=Show[ePlot,plot];
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
electronPolarization=CalculateElectronPolarization[stokes["p1"],stokes["p3_s2"]];
AppendTo[results,"P_e"->electronPolarization]
];

GetFourierCoefficientsFromFileNames


(* ::Subsection:: *)
(*Background Subtraction Methods*)


GetAverageElectronPolarizationWithBackgroudSubtractionSUM[fn_,fnBack_]:=Module[
{rawSignal,beamBackgroundSignal},
rawSignal=GetAverageCurrentNormalizedCountRateFromFilenames[fn];
beamBackgroundSignal=GetAverageCurrentNormalizedCountRateFromFilenames[fnBack];
ProcessElectronPolarizationFromSignal[rawSignal["signal"]-beamBackgroundSignal["signal"],Sqrt[rawSignal["error"]^2+beamBackgroundSignal["error"]^2]]
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


(* ::Section:: *)
(*Multiple Files*)


ProcessElectronPolarizationFileAverage[fn_List]:=Module[
{},
ProcessElectronPolarizationFileAverageSUM[fn]
];

ProcessElectronPolarizationBeamBackgroundFileAverage[fn_List]:=Module[
{},
fn
];

ProcessElectronPolarizationNoBeamBackgroundFileAverage[fn_List]:=Module[
{},
fn
];

SubtractElectronPolarizationBackground[signal_,background_]:=Module[
{},
signal
];


(* ::Subsection:: *)
(*Averaging Methods*)


ProcessElectronPolarizationFileAverageSUM[fn_List]:=Module[
	{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,datasets,currents,
	allCounts},
	results=<||>;
	f=ImportFile[fn];

	signals={};
	allCounts={};
	datasets={};
	currents={};
	fileNames={};
	For[i=1,i<=Length[f],i++,
		header=f[[i]][[1]];
		dataset=f[[i]][[2]];
		AppendTo[datasets,Normal[dataset]];
		counts=Normal[dataset[All,"COUNT"]];
		AppendTo[allCounts,ListPlot[counts]];
		current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
		AppendTo[currents,RepeatedMeasurementSummary[current]];
		dwellTime=header["DWELL(s)"];
		signal=GetCurrentNormalizedCountRate[counts,current,dwellTime];
		AppendTo[signals,signal];
		AppendTo[fileNames,header["File"]];
	];
	AppendTo[results,"averagedFiles"->fn];
	
	singleRunInfo=Transpose[{signals,datasets,allCounts}];
	threader[arg_]:=AssociationThread[{"signal","dataset","countsPlot"},arg];
	singleRunInfo=AssociationThread[fn,Map[threader,singleRunInfo]];
	For[i=1,i<=Length[singleRunInfo],i++,
	AppendTo[singleRunInfo[[i]],currents[[i]]];
	];

	signal=Mean[signals];
	AppendTo[results,"avgSignal"->signal];
	signalStdm=StandardDeviation[signals]/Sqrt[Length[signals]];
	AppendTo[results,"error"->signalStdm];
	
	AppendTo[results,"IndividualRunInfo"->singleRunInfo];
	PrependTo[results,ProcessElectronPolarizationFromSignal[signal,signalStdm]]
];

ProcessElectronPolarizationFileAverageFOURIER[fn_List]:=Module[
	{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results,dftResults,dfts,
	fcTitles,res,errors,avgFourierComponents,
	irResults (*individual run results*)},
	results=<||>;
	irResults=<||>;
	f=ImportFile[fn];

	signals={};
	fileNames={};
	For[i=1,i<=Length[f],i++,
		header=f[[i]][[1]];
		dataset=f[[i]][[2]];
		counts=Normal[dataset[All,"COUNT"]];
		current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
		dwellTime=header["DWELL(s)"];
		signal=GetCurrentNormalizedCountRate[counts,current,dwellTime];
		AppendTo[signals,signal];
		AppendTo[fileNames,header["File"]];
	];
	irResultsMapper[results_]:=AssociationThread[{"signal","fcCoefficients"},results];
	dfts=DFT/@signals;
	fcTitles={"Sin Coefficients","Cos Coefficients"};
	res=Map[StandardDeviation[Values[dfts[[All,#]]]]/Length[dfts]&,fcTitles];
	fcThreader[values_]:=AssociationThread[Keys[dfts[[1,1]]],values];
	res=Map[fcThreader,res];
	errors=AssociationThread[StringInsert[fcTitles," Error",-1],res];
	avgFourierComponents=Join[Mean[dfts],errors];
	AppendTo[irResults,AssociationThread[fn,Map[irResultsMapper,Transpose[{signals,dfts}]]]];
	AppendTo[results,"averagedFiles"->fn];
	AppendTo[results,<|ProcessElectronPolarizationFromFourierCoefficients[avgFourierComponents],"individualRunInfo"->irResults|>];
	<|fn[[1]]->results|>
];

ProcessElectronPolarizationFileAverageSTOKES[fn_List]:=Module[
{},
fn
];

ProcessElectronPolarizationFileAverageELECTRON[fn_List]:=Module[
{},
fn
];

AverageFourierCoefficients[fcList_Association]:=Module[
{},
fcList
];

AverageStokes[stokesVectors_]:=Module[{i,transposedValues,stokes,stokesStdDev},
stokes=stokesVectors[[1]];
stokesStdDev=stokes; (*Just creating another object with the same length as stokes *)
transposedValues=Transpose[stokesVectors];
For[i=1,i<=Length[stokesVectors[[1]]],i++,
stokes[[i]]=stokesNames[[i]]->Mean[transposedValues[[i]]];
stokesStdDev[[i]]=stokesErrNames[[i]]->StandardDeviation[transposedValues[[i]]]/Sqrt[Length[transposedValues]];
];
Join[Association[stokes],Association[stokesStdDev]]
];



GetAverageCurrentNormalizedCountRateFromFilenames[fn_]:=Module[
	{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results},
	results=<||>;
	f=ImportFile[fn];

	signals={};
	fileNames={};
	For[i=1,i<=Length[f],i++,
		header=f[[i]][[1]];
		dataset=f[[i]][[2]];
		counts=Normal[dataset[All,"COUNT"]];
		current=Abs[Normal[dataset[All,"CURRENT"]]]*10^(-header["SCALE"])/nA;
		dwellTime=header["DWELL(s)"];
		signal=GetCurrentNormalizedCountRate[counts,current,dwellTime];
		AppendTo[signals,signal];
		AppendTo[fileNames,header["File"]];
	];
	AppendTo[results,"averagedFiles"->fn];

	signal=Mean[signals];
	AppendTo[results,"signal"->signal];
	signalStdm=StandardDeviation[signals]/Sqrt[Length[signals]];
	AppendTo[results,"error"->signalStdm]
];

GetAverageCountRateFromFilenames[fn_]:=Module[{signals,i,f,header,dataset,
	counts,current,dwellTime,signal,signalStdm, (*Stdm stands for standard deviation of the mean *)
	fileNames,results},
	results=<||>;
	f=ImportFile[fn];

	signals={};
	fileNames={};
	For[i=1,i<=Length[f],i++,
		header=f[[i]][[1]];
		dataset=f[[i]][[2]];
		counts=Normal[dataset[All,"COUNT"]];
		dwellTime=header["DWELL(s)"];
		signal=GetCurrentNormalizedCountRate[counts,dwellTime];
		AppendTo[signals,signal];
		AppendTo[fileNames,header["File"]];
	];
	AppendTo[results,"averagedFiles"->fn];

	signal=Mean[signals];
	AppendTo[results,"signal"->signal];
	signalStdm=StandardDeviation[signals]/Sqrt[Length[signals]];
	AppendTo[results,"error"->signalStdm]
];

