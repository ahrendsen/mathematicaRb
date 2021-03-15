(* ::Package:: *)

(* ::Title:: *)
(*Absorption Profile Functions*)


(* ::Chapter:: *)
(*Initialization*)


detuningCT="DET"; (* Detuning (C)olumn (T)itle *)
absWavelengthColumnName="WAV";
signalCT={"PMP","PRB","REF"}; (* The old columnTitles were these. I'm saving them here just in case I need to process an old file. *)
signalCT={"HORIZ","VERT","REF"};

SetDirectory[$PACKAGES];
Import["fileManipulation.wl"];
Import["physicalConstants.wl"];
ResetDirectory[];


(* ::Chapter:: *)
(*Labels*)


absLabels={FrameLabel->{"Detuning (GHz)","Intensity (normalized)"}}


(* ::Chapter:: *)
(*Data Cleaning*)


ApproximateFrequency[data_,columnMissing_,columnReference_,
						missingEntry_]:=
Module[{j,adjacentSpots,referenceSpacing,
		validEntries,fit}(*These are your local variables that aren't stored beyond the scope of the module.*),
		
(*Selects the data that is certainly an appropriate measure of wavelength *)
adjacentSpots={};
validEntries=data[All,{columnReference,columnMissing}][Select[#[columnMissing]>793&]];
(*The "reference column" is the column that is monotonically increasing, we need to know what the size of spaces between these values is *)
referenceSpacing=2;
j=1;
(* This "While" loop accumulates adjacent data points used to approximate the missing frequency, we only want the closest datapoints. *)
While[And[Length[adjacentSpots]<3,Abs[j*referenceSpacing]<117](* We want a few data points to approximate where our point lies, keep going until we have at least 3*),
adjacentSpots=validEntries[Select[Abs[#[columnReference]-missingEntry]<=referenceSpacing*j&]];
j++;
];
fit=LinearModelFit[Normal[adjacentSpots[Values]],x,x];
fit[missingEntry] (* Last line is the "return" statement, all other lines need a semi-colon*)
];

(* This is called on a dataset to fill in all of the missing frequency gaps.*)
FillInMissedFrequencies[rawdata_,columnMissing_,columnReference_]:=Module[{referenceSpacing,missingEntries,validEntries,adjacentSpots,estimatedFrequency,linearModel,dataCopy,dataReturn,position,approxValue,data2,k},
dataCopy=Dataset[rawdata];
dataReturn=Normal[dataCopy];
missingEntries=dataCopy[Select[#[columnMissing]<793&]];
missingEntries=Flatten[Normal[missingEntries[All,{columnReference}][Values]]];
For[k=1,k<=Length[missingEntries],k++,
approxValue=ApproximateFrequency[dataCopy,columnMissing,columnReference,missingEntries[[k]]];
(* Find the index of the missing value in the list *)
position=Position[dataCopy,missingEntries[[k]]][[1,1]];
(* Use the index to insert the approximate value in the data  *)
dataReturn[[position,columnMissing]]=approxValue;
];
dataReturn
];

CleanDataset[dataset_,missingColumn_,referenceColumn_]:=Module[{d,d2,i},
d=FillInMissedFrequencies[Normal[dataset[All,{referenceColumn,missingColumn}]],missingColumn,referenceColumn];
d2=Normal[dataset];
For[i=1,i<=Length[d2],i++,
AppendTo[d2[[i]],absWavelengthColumnName->d[[i]][absWavelengthColumnName]]
];
Dataset[d2]
];


RemoveBackgroundFromSaturatedProfile[absorptionDataset_]:=Module[{dataset,noAbsNeg,noAbsPos,wings,linearModelWings,intensitySlope,intensityIntercept,detuning,intensity,x},
{detuning,intensity}=Transpose[absorptionDataset];
intensity=intensity-Min[intensity];
Transpose[{detuning,intensity}
]];

SlopeInterceptProfileConversion[dataset_,slope_,intercept_]:=
Module[
{},
dataset[All,<|detuningCT->#AOUT*slope+intercept,signalCT[[1]]->#[signalCT[[1]]],signalCT[[2]]->#[signalCT[[2]]],signalCT[[3]]->#[signalCT[[3]]]|>&]
];

TrimDataByDetuning[dataset_]:=Module[{lowCut,highCut,temp},
lowCut=5 ;
highCut=30;
temp=dataset[Select[Abs[#[detuningCT]]>lowCut&]];
temp[Select[Abs[#[detuningCT]]<highCut&]]
];


(* ::Chapter:: *)
(*Data Summary Text Form*)


GetTemperaturesOverTime[importedFiles_]:=Module[
{justHeaderData,resTempVsTime,cellTempVsTime},
justHeaderData=Dataset[Map[#[[1]]&,importedFiles]];
resTempVsTime=justHeaderData[All,{"Time","T_res"}]//Values;
cellTempVsTime=justHeaderData[All,{"Time","T_trg"}]//Values;
<|"CellTemperature"->cellTempVsTime,"ResTemperature"->resTempVsTime|>
];


(* ::Chapter:: *)
(*Plotting*)


transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},{{transGHz[[2]],0},{transGHz[[2]],10}},{{transGHz[[3]],0},{transGHz[[3]],10}},{{transGHz[[4]],0},{transGHz[[4]],10}}}];

PlotSignalProfile[fileName_]:=Module[
{datasetRaw,headerInfo,datasetZeroLineRemoved,
dataset,plotList,transitionLines,rawFile,
headerStrippedData,columnHeaderLineNumber,
j,k,ass,tabularData,ds,rbAbsorptionRef,
rbAbsorptionPro,rbAbsorptionPmp},

datasetRaw=GetFileDataset[fileName];
headerInfo=GetFileHeaderInfo[fileName];
dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-1.06|>&];
plotList={Normal[dataset[All,{#[detuningCT],#[signalCT[[2]]]}&]],Normal[dataset[All,{#[detuningCT],#[signalCT[[1]]]}&]]};
transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},{{transGHz[[2]],0},{transGHz[[2]],10}},{{transGHz[[3]],0},{transGHz[[3]],10}},{{transGHz[[4]],0},{transGHz[[4]],10}}}];
ListPlot[plotList,Joined->True,PlotRange->{{-8,8},Full},Epilog->Style[transitionLines,Thickness[.005],Dotted],PlotLabel->"Absorption Profile",FrameLabel->{"Detuning(GHz)","Intensity"}]
];

PlotSignalProfileComparison[fileNames_:List,plotOptions_:{}]:=Module[

{datasetRaw,headerInfo,datasetZeroLineRemoved,
dataset,plotList,transitionLines,rawFile,
headerStrippedData,columnHeaderLineNumber,
j,ass,tabularData,ds,rbAbsorptionRef,
rbAbsorptionPro,rbAbsorptionPmp,orderedPairs},

plotList={};
For[j=1,j<=Length[fileNames],j++,
datasetRaw=GetFileDataset[fileNames[[j]]];
headerInfo=GetFileHeaderInfo[fileNames[[j]]];
dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-1.0|>&];
orderedPairs=Normal[dataset[All,{#[detuningCT],#[signalCT[[2]]]}&]];
AppendTo[plotList,orderedPairs];
];

transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},{{transGHz[[2]],0},{transGHz[[2]],10}},{{transGHz[[3]],0},{transGHz[[3]],10}},{{transGHz[[4]],0},{transGHz[[4]],10}}}];
ListPlot[plotList,
	Joined->True,PlotRange->{{-5,8},Full},
	Epilog->Style[transitionLines,Thickness[.005],Dotted],
	PlotLabel->"Absorption Profile",FrameLabel->{"Detuning(GHz)","Intensity"},
	plotOptions]
];

PlotSignalProfileDifference[fileName1_,fileName2_,plotOptions_:{}]:=Module[
{datasetRaw,headerInfo,datasetZeroLineRemoved,
dataset,intensity,transitionLines,rawFile,fileNames,
headerStrippedData,columnHeaderLineNumber,
j,ass,tabularData,ds,rbAbsorptionRef,
rbAbsorptionPro,rbAbsorptionPmp,orderedPairs,
detunings,differences,timeStrings},

fileNames={fileName1,fileName2};
intensity={};
For[j=1,j<=Length[fileNames],j++,
datasetRaw=GetFileDataset[fileNames[[j]]];
headerInfo=GetFileHeaderInfo[fileNames[[j]]];
dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-wavemeterOffset|>&];
orderedPairs=Normal[dataset[All,{#[detuningCT],#[signalCT[[2]]]}&]];
AppendTo[intensity,Transpose[orderedPairs][[2]]];
detunings=Transpose[orderedPairs][[1]];
];

differences=intensity[[2]]-intensity[[1]];
timeStrings=Map[GetTimeStringFromFileNameString,fileNames];

transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},{{transGHz[[2]],0},{transGHz[[2]],10}},{{transGHz[[3]],0},{transGHz[[3]],10}},{{transGHz[[4]],0},{transGHz[[4]],10}}}];
ListPlot[Transpose[{detunings,differences}],
	Joined->True,PlotRange->{{-5,8},Full},
	Epilog->Style[transitionLines,Thickness[.005],Dotted],
	PlotLabel->timeStrings[[2]] <> " Minus "<> timeStrings[[1]],FrameLabel->{"Detuning(GHz)","Intensity Difference"},
	plotOptions]
];

PlotSignalProfileNormalized[fileName_]:=Module[{datasetRaw,headerInfo,datasetZeroLineRemoved,dataset,plotList,transitionLines,rawFile,headerStrippedData,columnHeaderLineNumber,j,k,ass,tabularData,ds,rbAbsorptionRef,rbAbsorptionPro,rbAbsorptionPmp},

datasetRaw=GetFileDataset[fileName];
headerInfo=GetFileHeaderInfo[fileName];

(*datasetZeroLineRemoved=Select[datasetRaw,#[absWavelengthColumnName]>793&];*)
	
(*dataset=datasetZeroLineRemoved[All,<|#,detuningCT->c/100/(#[absWavelengthColumnName])-\[Nu]0|>&] ;*)
dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-1.0|>&];
dataset=NormalizeProfileWithWings[dataset,detuningCT,signalCT[[2]],-5,6];
plotList=Normal[dataset[All,{#[detuningCT],#[signalCT[[2]]<>"NORMWW"]}&]](*,Normal[dataset[All,{#[detuningCT],#[signalCT[[1]]<>"NORMWW"]&}]]*)
(*transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},{{transGHz[[2]],0},{transGHz[[2]],10}},{{transGHz[[3]],0},{transGHz[[3]],10}},{{transGHz[[4]],0},{transGHz[[4]],10}}}];*)
ListPlot[plotList,Joined->True,PlotRange->{{-5,8},{0,1.1}},Epilog->Style[transitionLines,Thickness[.005],Dotted],PlotLabel->"Absorption Profile",FrameLabel->{"Detuning(GHz)","Intensity"}]
];

GetReferenceOrderedPairs[fileName_,background_,scale_]:=Module[{datasetRaw,headerInfo,datasetZeroLineRemoved,dataset,plotList,transitionLines,rawFile,headerStrippedData,columnHeaderLineNumber,j,k,ass,tabularData,ds,rbAbsorptionRef,rbAbsorptionPro,rbAbsorptionPmp},

datasetRaw=GetFileDataset[fileName];
headerInfo=GetFileHeaderInfo[fileName];

(*datasetZeroLineRemoved=Select[datasetRaw,#[absWavelengthColumnName]>793&];*)
	
(*dataset=datasetZeroLineRemoved[All,<|#,detuningCT->c/100/(#[absWavelengthColumnName])-\[Nu]0|>&] ;*)
dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-1.0|>&];
Normal[dataset[All,{#[detuningCT],scale*(#[signalCT[[3]]]-background)}&]]
];

PlotReferenceProfile[fileName_,background_,scale_]:=Module[{datasetRaw,headerInfo,datasetZeroLineRemoved,dataset,plotList,transitionLines,rawFile,headerStrippedData,columnHeaderLineNumber,j,k,ass,tabularData,ds,rbAbsorptionRef,rbAbsorptionPro,rbAbsorptionPmp},

datasetRaw=GetFileDataset[fileName];
headerInfo=GetFileHeaderInfo[fileName];

(*datasetZeroLineRemoved=Select[datasetRaw,#[absWavelengthColumnName]>793&];*)
	
(*dataset=datasetZeroLineRemoved[All,<|#,detuningCT->c/100/(#[absWavelengthColumnName])-\[Nu]0|>&] ;*)
dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-1.0|>&];
plotList={Normal[dataset[All,{#[detuningCT],scale*(#[signalCT[[3]]]-background)}&]]};
transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},{{transGHz[[2]],0},{transGHz[[2]],10}},{{transGHz[[3]],0},{transGHz[[3]],10}},{{transGHz[[4]],0},{transGHz[[4]],10}}}];
ListPlot[plotList,Joined->True,PlotRange->{{-5,8},Full},Epilog->Style[transitionLines,Thickness[.005],Dotted],PlotLabel->"Absorption Profile",FrameLabel->{"Detuning(GHz)","Intensity"}]
];

QuickPlotRbProfileWithRange[dataset_,plotRange_]:=Module[{rawFile,headerStrippedData,columnHeaderLineNumber,j,k,ass,tabularData,ds,rbAbsorptionRef,rbAbsorptionPro,rbAbsorptionPmp},
rbAbsorptionRef=Normal[dataset[All,{"VOLT",signalCT[[3]]}][[All,Values]]];
rbAbsorptionPro=Normal[dataset[All,{"VOLT",signalCT[[1]]}][[All,Values]]];
rbAbsorptionPmp=Normal[dataset[All,{"VOLT",signalCT[[2]]}][[All,Values]]];
ListPlot[{rbAbsorptionPro,rbAbsorptionRef,rbAbsorptionPmp},PlotRange->{plotRange,Full},Ticks-> {Range[0,1023,20],Automatic},GridLines->{Range[0,1023,20],None},ImageSize->{Full,Automatic}]
];

CompareProfiles[fileName1_,fileName2_]:=
Module[{rawFile,headerStrippedData,columnHeaderLineNumber,j,k,ass,tabularData,ds,profile1,profile2},
profile1=Normal[ConvertAoutProfileToDetuning[GetRbProfileDataset[fileName1],165,"D","A"][All,{detuningCT,signalCT[[3]]}][[All,Values]]];
profile2=Normal[ConvertAoutProfileToDetuning[GetRbProfileDataset[fileName2],165,"D","A"][All,{detuningCT,signalCT[[3]]}][[All,Values]]];
ListPlot[{profile1,profile2},PlotRange->{Full,Full},Ticks-> {Range[0,1023,20],Automatic},GridLines->{Range[0,1023,20],None},ImageSize->{Full,Automatic}]
];

(* Uses the far detunings of the profile, the "wings", to establish the unabsorbed intensity and normalize the profile *)
NormalizeProfile[absorptionDataset_]:=Module[{dataset,noAbsNeg,noAbsPos,wings,linearModelWings,intensitySlope,intensityIntercept,detuning,intensity,x},
dataset=absorptionDataset;
noAbsNeg=dataset[Select[#[detuningCT]<-4.5&]];(*Selects only data with detunings less than -4.5, where no absorption is present. *)
noAbsNeg=noAbsNeg[All,{detuningCT,signalCT[[1]]}];(*Removes all data from dataset except the detunings and the intensity of the probe laser *)
noAbsNeg=Normal[noAbsNeg[[All,Values]]]; (* Removes the headers from the dataset, giving a matrix of the ordered pairs of the data *)

noAbsPos=dataset[Select[#[detuningCT]>6&]];(*See the above comments*)
noAbsPos=noAbsPos[All,{detuningCT,signalCT[[1]]}];
noAbsPos=Normal[noAbsPos[[All,Values]]];

wings=Join[Normal[noAbsNeg,noAbsPos]];
linearModelWings=LinearModelFit[wings,x,x];
intensitySlope=linearModelWings["BestFitParameters"][[2]];
intensityIntercept=linearModelWings["BestFitParameters"][[1]];
detuning=Normal[absorptionDataset[All,detuningCT]];
intensity=Normal[absorptionDataset[All,signalCT[[1]]]];
intensity=intensity/(intensityIntercept+intensitySlope*detuning);
(*Transpose[{detuning,intensity};*)
absorptionDataset[All,<|#,"PRBNRMWWING"-> #[signalCT[[1]]]/(intensityIntercept+intensitySlope*#[detuningCT])|>&]
];

ProcessFullyAbsorbedProfile[fileName_,columnTitle_:signalCT[[2]]]:=Module[
{datasetRaw,headerInfo,dataset,
plotList2,transitionLines,graph},

datasetRaw=GetFileDataset[fileName];
headerInfo=GetFileHeaderInfo[fileName];

dataset=datasetRaw[All,<|#,detuningCT->#[detuningCT]-1.06|>&];

FindBackgroundIntensity[dataset,detuningCT,columnTitle,transGHz[[2]]];
dataset=RemoveBackgroundIntensity[dataset,detuningCT,columnTitle,transGHz[[2]]];

dataset=NormalizeProfileWithWings[dataset,detuningCT,columnTitle,-5,6.5];

plotList2=Normal[dataset[All,{#[detuningCT],#[columnTitle<>"NORMWW"]}&]];
transitionLines=Line[{{{transGHz[[1]],0},{transGHz[[1]],10}},
						{{transGHz[[2]],0},{transGHz[[2]],10}},
						{{transGHz[[3]],0},{transGHz[[3]],10}},
						{{transGHz[[4]],0},{transGHz[[4]],10}}}];
graph=ListPlot[plotList2,Joined->True,PlotRange->{{-5,8},Full},Epilog->Style[transitionLines,Thickness[.005],Dotted],PlotLabel->"Absorption Profile",FrameLabel->{"Detuning(GHz)","Intensity"}];
<|"graph"->graph,"data"->plotList2|>
];



FindWavemeterOffset[dataset_]:=Module[{CtoDDiff,BtoCDiff,AtoBDiff,aoutSteps,s,startPos,rbAbsorptionRef,scanStepSize,posLower,posUpper,i,j,k,transAout,peak,aoutPos,aoutPrev,aoutStep,aoutNextApprox,aoutConvData,lm,aoutConv,aoutIntercept},

(*Use the position in the list to find the first absorption peak's value*)
scanStepSize=Abs[(rbAbsorptionRef[[startPos]]-rbAbsorptionRef[[startPos+1]])[[1]]];
posLower=startPos-5;
posUpper=startPos+5;
transAout={0,0,0,0};
i=0;

If[transitionLetter=="D",j=4,If[transitionLetter=="C",j=3,If[transitionLetter=="B",j=2,If[transitionLetter=="A",j=1,j=0]]]];
If[lastTransitionLetter=="D",k=4,If[lastTransitionLetter=="C",k=3,If[lastTransitionLetter=="B",k=2,If[lastTransitionLetter=="A",k=1,k=0]]]];
(*Cycle through and find the remaining transisition "peaks"*)
While[And[j>= k],
peak=Min[Take[rbAbsorptionRef,{posLower,posUpper},All]];
(*Use the value of the first absorption peak to find its position in the list*)
aoutPos=Position[rbAbsorptionRef,peak][[1,1]];
(*Use the position in the list to get the Aout value located there*)
transAout[[j]]=rbAbsorptionRef[[aoutPos,1]];
(*Because the profile is consistent,we can guess at where the next min.will be*)
If[j>k,aoutPrev=transAout[[j]];aoutStep=aoutSteps[[j-1]];
aoutNextApprox=aoutPrev+aoutStep;
s=Nearest[Take[rbAbsorptionRef,All,1],aoutNextApprox];
posUpper=aoutPos+Floor[aoutStep/scanStepSize]+2*Floor[10/scanStepSize];
posLower=aoutPos+Floor[aoutStep/scanStepSize]-2*Floor[10/scanStepSize];];
j--;];

(*Remove the transistions that we don't have from the list so that we can estimate a linear function relating the data points that we do have*)
aoutConvData=Dataset[Transpose[{transAout,transGHz}]]; 
aoutConvData=aoutConvData[Select[Abs[#[[1]]]>0&]];
aoutConvData=Normal[aoutConvData];


(*Make a linear fit on the data to obtain an AOUT\[Rule]GHz Conversion*)
lm=LinearModelFit[aoutConvData,x,x];
Print[lm];
dataset[All,<|"AOUT"->Key[aoutColumnName],absWavelengthColumnName-> Key[wavelengthColumnName],detuningCT->Key[aoutColumnName]/* lm,signalCT[[1]]->probeColumnName,signalCT[[2]]->pumpColumnName,signalCT[[3]]->referenceColumnName|>]

];


PlotVoltageVsDetuningFromDataset[file_]:=Module[{d,h,s,fName},
	d=file[[2]]; (* data *)
	h=file[[1]]; (* header *)
	s=h["Filename"]; (* string (for filename) *)
	fName=First[StringCases[s,RegularExpression["RbAbs.*\\.dat"]]];
	ListPlot[d[All,{"VOLT","DET"}],FrameLabel->{"Piezo Setting (V)","Reported Detuning (GHz)"},PlotLabel->fName]
];

PlotVoltageVsDetuning[fileName_]:=Module[{file},
file=ImportFile[fileName];
PlotVoltageVsDetuningFromDataset[file]
];
SetAttributes[PlotVoltageVsDetuning,Listable];



PlotRbAbsSeries[importedFiles_]:=Module[
{ordColName,datasets,dataTable,normalizedDataTables,
fileNames,legendPlace},
legendPlace={.75,.25};
legendPlace=Right;
ordColName="VERT";
datasets=importedFiles;
fileNames=Map[FileNameTake,Normal[Dataset[Map[#[[1]]&,datasets]][All,"File"]]];
dataTable=Map[Normal[#[[2]]]&,datasets];
normalizedDataTables=Map[Normal[NormalizeProfileWithWings[Dataset[#],"DET",ordColName]]&,dataTable];
ListPlot[Dataset[normalizedDataTables][All,All,{#["DET"]-wavemeterOffset&,ordColName<>"NORMWW"}],Joined->True,PlotRange->{{-5,6},{-.1,1.1}},,PlotLegends->Placed[LineLegend[fileNames,LegendFunction->legendFrame],Right],
FrameLabel->{"Detuning (GHz)","Intensity (Normalized with Wings)"}]
];


legendFrame[legend_]:=Framed[legend,FrameStyle->Black,RoundingRadius->0,FrameMargins->0,Background->White]
PlotTemperaturesOverTime[importedFiles_]:=Module[
{justHeaderData,resTempVsTime,cellTempVsTime},
justHeaderData=Dataset[Map[#[[1]]&,importedFiles]];
resTempVsTime=justHeaderData[All,{"Time","T_res"}]//Values;
cellTempVsTime=justHeaderData[All,{"Time","T_trg"}]//Values;
DateListPlot[{cellTempVsTime,resTempVsTime},
FrameLabel->{DateString[resTempVsTime[[1]][[1]],{"ISODate"}]<>", HH:MM",
"Temperature (\[Degree]C)"},
PlotLegends->Placed[LineLegend[{"Collision Cell Temperature","Reservoir Temperature"},LegendFunction->legendFrame],{.75,.25}]
]
];



(* ::Subchapter:: *)
(*Normalization*)


(* Normalizes the data by first substracting away the minValue, which is the intensity 
when the probe light is fully absorbed, then dividing by the max value, which is the intensity
far away from the abosrption profile. *)
NormalizeDataTo1[dataset_,minValue_,maxValue_]:=dataset[All,<|#,"NORMTO1"-> (#[signalCT[[1]]]/#[signalCT[[3]]]-minValue)/(maxValue-minValue)|>&];

(* Selects intensity data far enough away from the transitions so that the 
intensity is still above 95% of far detuned. Then it calculates a linear slope from
this far detuned data and uses the value of the linear fit function to determine 
the number to divide by to normalize the intensity to 1. *)
Clear[NormalizeProfileWithWings];
NormalizeProfileWithWings[dataset_,abscissaColumnName_,ordinateColumnName_,
						  lowCutoff_:-6,highCutoff_:9,lowMax_:-20,highMax_:20]:=
Module[{detuningColumnName,datasetZeroLineRemoved,
		wingUpper,wingLower,wingsDataset,
		nlmFit,a,b,d,x,newList,i},
wingUpper=dataset[Select[#[abscissaColumnName]> highCutoff&&#[abscissaColumnName]< highMax&]];
wingLower=dataset[Select[#[abscissaColumnName]< lowCutoff&&#[abscissaColumnName]> lowMax&]];
wingsDataset=Join[wingUpper,wingLower][SortBy[abscissaColumnName]];

nlmFit=NonlinearModelFit[Normal[wingsDataset[All,{abscissaColumnName,ordinateColumnName}][All,Values]],a x^2+b x+d,{a,b,d},x];
dataset[All,<|#,ordinateColumnName<>"NORMWW"->#[ordinateColumnName]/nlmFit[#[abscissaColumnName]]|>&]
];

RemoveBackgroundIntensity[dataset_,abcissaColumnName_,ordinateColumnName_,abcissaPoint_]:=
Module[{absPoint,positions,i,
		fullAbsorptions,minForNorm},
minForNorm=FindBackgroundIntensity[dataset,abcissaColumnName,ordinateColumnName,abcissaPoint];
dataset[All,<|#,ordinateColumnName->#[ordinateColumnName]-minForNorm|>&]
];

FindBackgroundIntensity[dataset_,abcissaColumnName_,ordinateColumnName_,abcissaPoint_]:=
Module[{absPoint,positions,i,
		fullAbsorptions,minForNorm},
absPoint=Nearest[Normal[dataset[All,abcissaColumnName]],abcissaPoint,3];
positions={};
For[i=1,i<=Length[absPoint],i++,
AppendTo[positions,Position[Normal[dataset[All,abcissaColumnName]],absPoint[[i]]]];
];
fullAbsorptions=Normal[dataset[All,ordinateColumnName]][[Flatten[positions]]];
minForNorm=Mean[fullAbsorptions]
];


(* ::Chapter:: *)
(*Theory*)


(* Returns Pressure in torr*)
RbPressure[t_]:=1*10^(15.88253 - 4529.635/t +.00058663*t-2.99138*Log10[t])
RbNDensity[t_]:=RbPressure[t]*133.323/(kB*t)/1*^6
