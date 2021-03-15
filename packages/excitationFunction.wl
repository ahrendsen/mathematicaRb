(* ::Package:: *)

(* ::Title:: *)
(*Excitation Function*)


(* ::Chapter:: *)
(*Initialization*)


SetDirectory[$PACKAGES];
Import["fileManipulation.wl"];
Import["physicalConstants.wl"];
ResetDirectory[];

logPlotZero=.0001;
potentialsLineThickness=.005;
sweepLineThickness=.0001;
rfaPlotOptions={FrameLabel-> {"He Potential (V)","Transmitted Current (nA)"}};
potentialsEpilog:={Thickness[potentialsLineThickness],Dotted,potentials[[1]],DotDashed,potentials[[2]]};


(* ::Chapter:: *)
(*Equations To Work With Data*)


dCT(*Data column Titles *)=<|
"bias"->"V_fil",
"scale"->"MagnitudeOfCurrent(*10^-X)",
"V_N2"->"V_N2",
"V_he"->"V_he",
"V_sw"->"V_sw",
"V_1D"->"V_1D",
"V_2A"->"V_2A",
"see"->"SecondaryElectronEnergy",
"I_f"->"I_f"|>;

MinMaxAbcissa[ordPair_]:=MinMax[Transpose[ordPair][[1]]];
MinMaxOrdinate[ordPair_]:=MinMax[Transpose[ordPair][[2]]];

PlotExFn[fileName_,primaryAbcissa_:dCT["V_he"],secondaryAbcissa_:dCT["see"]]:=
Module[{rawImport,d,h,divisions,
exp,vPR,newDat,
secondaryGoodDiv,topAxis,
potentials,plot,plotNoLabel},
rawImport=ImportFile[fileName];
d=rawImport[[2]];
h=rawImport[[1]];

exp=ToExpression[h["MagnitudeOfCurrent(*10^-X)"]];
newDat=d[All,<|#,dCT["I_f"]-> Around[#[dCT["I_f"]],Max[#["I_fStDev"],1/1024]]*10^(-exp)|>&];

PlotExFnFromDataset[newDat,primaryAbcissa]
];
SetAttributes[PlotExFn,Listable];

PlotExFnRaw[fileName_,primaryAbcissa_:dCT["V_he"],
			secondaryAbcissa_:dCT["see"]]:=
Module[{rawImport,d,h,divisio,ns,
vPR,exp,newDat,
secondaryGoodDiv,topAxis,
exfnPlotRaw,exfnPlotRawNoLabel},
rawImport=ImportFile[fileName];
	d=rawImport[[2]];
	h=rawImport[[1]];

	exp=ToExpression[h["MagnitudeOfCurrent(*10^-X)"]];
	newDat=d[All,<|#,dCT["I_f"]-> Around[#[dCT["I_f"]],Max[#["I_fStDev"],1/1024]]*10^(-exp)|>&];
	PlotExFnRawFromDataset[newDat,primaryAbcissa]
];
SetAttributes[PlotExFnRaw,Listable];


PlotRFA[fileName_,
		primaryAbcissa_:dCT["V_he"],
		secondaryAbcissa_:dCT["see"]]:=
Module[{rawImport,d,h,divisions,mmSA,vPR,
mmPA,primaryGoodDiv,secondaryGoodDiv,newDat,
topAxis,plot,plotNoLabel,exp,potentials},
rawImport=ImportFile[fileName];
	d=rawImport[[2]];
	h=rawImport[[1]];

	exp=ToExpression[h["MagnitudeOfCurrent(*10^-X)"]];
	newDat=d[All,<|#,dCT["I_f"]-> Around[#[dCT["I_f"]],Max[#["I_fStDev"],1/1024]]*10^(-exp)|>&];

	PlotRFAFromDataset[newDat]
];
SetAttributes[PlotRFA,Listable];

ProcessExFnFile[fileName_]:=Module[{f,h,d,plot,rawPlot,rfaPlot,results},
results=<||>;
f=ImportFile[fileName];
h=f[[1]];
d=f[[2]];
plot=PlotExFn[fileName];
AppendTo[results,plot];
rawPlot=PlotExFnRaw[fileName];
AppendTo[results,rawPlot];
rfaPlot=PlotRFA[fileName];
AppendTo[results,rfaPlot];
AppendTo[results,<|"dataset"->d|>];
AppendTo[results,h];
results
];
SetAttributes[ProcessExFnFile,Listable];

Clear[ProcessExFnFileSequence];
ProcessExFnFileSequence[fileNames_]:=Module[
{combinedData,results,h},
results=<||>;
combinedData=CombineExcitationFunctionData[fileNames];
AppendTo[results,"Date"->ExtractTimeInfoFromFileNameString[fileNames[[1]]]];
h=ImportFile[fileNames[[1]]][[1]];
AppendTo[results,PlotExFnFromDataset[combinedData]];
AppendTo[results,PlotExFnRawFromDataset[combinedData]];
AppendTo[results,PlotRFAFromDataset[combinedData]];
AppendTo[results,"Filenames"->fileNames];
AppendTo[results,"rawData"->Normal[combinedData]];
AppendTo[results,h]
];

TwoAxisListPlot[{f_, g_}] := 
 Module[{fgraph, ggraph, frange, grange, fticks, 
   gticks}, {fgraph, ggraph} = 
   MapIndexed[
    ListPlot[#, Axes -> True, 
      PlotStyle -> ColorData[1][#2[[1]]]] &, {f, g}]; {frange, 
    grange} = (PlotRange /. AbsoluteOptions[#, PlotRange])[[
      2]] & /@ {fgraph, ggraph}; fticks = N@FindDivisions[frange, 5]; 
  gticks = Quiet@
    Transpose@{fticks, 
      ToString[NumberForm[#, 2], StandardForm] & /@ 
       Rescale[fticks, frange, grange]}; 
  Show[fgraph, 
   ggraph /. 
    Graphics[graph_, s___] :> 
     Graphics[
      GeometricTransformation[graph, 
       RescalingTransform[{{0, 1}, grange}, {{0, 1}, frange}]], s], 
   Axes -> False, Frame -> True, 
   FrameStyle -> {ColorData[1] /@ {1, 2}, {Automatic, Automatic}}, 
   FrameTicks -> {{fticks, gticks}, {Automatic, Automatic}}]]


(* ::Chapter:: *)
(*Stitching Data Together*)


(* Take input of list of file names, output stitched data and graphs of excitation functions. *)
CombineExcitationFunctionData[files_]:=Module[
{constantsFromDataToInclude,constantsFromHeaderToInclude,
i,header, magnitude, dat, newDat, 
allDat,currentData,countData,
currentPlot,countPlot,
output,sweep,filBias},

For[i=1,i<=Length[files],i++,
header=GetFileHeaderInfo[files[[i]]];
magnitude=header["MagnitudeOfCurrent(*10^-X)"];
dat=GetFileDataset[files[[i]]];
dat=dat[All,<|#,dCT["I_f"]-> #[dCT["I_f"]]*10^(-magnitude)|>&];
allDat=If[i==1,dat,Join[allDat,dat]];
];
allDat
];


(* ::Chapter:: *)
(*Private Functions*)


PlotExFnFromDataset[dataset_,primaryAbcissa_:dCT["V_he"]]:=Module[{rawImport,h,divisions,
exp,vPR,d,dataPoints,
secondaryGoodDiv,topAxis,
potentials,potentialLines,plot,plotNoLabel,
plotFull,plotFullNoLabel},
d=dataset;

(*potentials=          <|dCT["bias"]->d[[1]][dCT["bias"]],
				        "n2"->d[[1]][dCT["V_N2"]]+d[[1]][dCT["bias"]],
(*Helium Threshhold *) "HeT"->d[[1]][dCT["V_N2"]]+d[[1]][dCT["bias"]]+23,
				       "swp"->d[[1]][dCT["V_sw"]]|>;
				       *)
				       
potentials=          <|dCT["bias"]->d[[1]][dCT["bias"]],
				        "n2"->d[[1]][dCT["V_N2"]]+d[[1]][dCT["bias"]],
(*Helium Threshhold *) "HeT"->d[[1]][dCT["V_N2"]]+d[[1]][dCT["bias"]]+23,
				       "V_1D"->d[[1]][dCT["V_1D"]],
				       "V_2A"->d[[1]][dCT["V_2A"]]|>;
				       
potentialLines={
Line[{
{potentials[dCT["bias"]],Log[logPlotZero]},
{potentials[dCT["bias"]],3000}
}],
Line[{
{potentials["n2"],Log[logPlotZero]},
{potentials["n2"],3000}
}],
Line[{
{potentials["HeT"],Log[logPlotZero]},
{potentials["HeT"],3000}
}],
Line[{
{potentials["HeT"]+potentials["V_2A"],Log[logPlotZero]},
{potentials["HeT"]+potentials["V_2A"],3000}
}],
Line[{
{potentials["HeT"]-potentials["V_1D"],Log[logPlotZero]},
{potentials["HeT"]-potentials["V_1D"],3000}
}]
};
	(* It would be really nice to plot the bottom axis with the total He offset, and the top axis with the secondary electron energy. What follows is attempting to accomplish that *)

	vPR={.1,50};(*Vertical Plot Range *)
	dataPoints=d[All,{primaryAbcissa,(#["CountRate"]-23)/(#[dCT["I_f"]]/nA)&}];
	plotFull=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{{Automatic,Automatic},vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
		FrameLabel->{{"Counts/nA",Automatic},{"Total He Offset (V)",Automatic}},
					Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentialLines[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentialLines[[2]],
		colorBlindPalette[[3]],
		potentialLines[[3]],
		Thickness[sweepLineThickness],
		potentialLines[[4]],
		potentialLines[[5]]
		}
	];
	plotFullNoLabel=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{Automatic,vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
					Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentialLines[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentialLines[[2]],
		colorBlindPalette[[3]],
		potentialLines[[3]],
		Thickness[sweepLineThickness],
		potentialLines[[4]],
		potentialLines[[5]]
		}
	];
	
	
	plot=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]+5,Automatic},vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
		FrameLabel->{{"Counts/nA",Automatic},{"Total He Offset (V)",Automatic}},
					Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentialLines[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentialLines[[2]],
		colorBlindPalette[[3]],
		potentialLines[[3]],
		Thickness[sweepLineThickness],
		potentialLines[[4]],
		potentialLines[[5]]
		}
	];
	plotNoLabel=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]+5,Automatic},vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
					Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentialLines[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentialLines[[2]],
		colorBlindPalette[[3]],
		potentialLines[[3]],
		Thickness[sweepLineThickness],
		potentialLines[[4]],
		potentialLines[[5]]
		}
	];
	<|"exFnFullPlot"->plotFull,"exFnFullPlotNoLabel"->plotFullNoLabel,"exFnPlot"->plot,"exFnPlotNoLabel"->plotNoLabel,"exFnData"->Normal[dataPoints]|>
];


PlotExFnRawFromDataset[dataset_,primaryAbcissa_:dCT["V_he"]]:=Module[
{rawImport,d,h,divisions,
exp,vPR,dataPoints,
secondaryGoodDiv,topAxis,
potentials,potentialLines,
exfnPlotRaw,exfnPlotRawNoLabel},

d=dataset;
potentials=          <|dCT["bias"]->d[[1]][dCT["bias"]],
				        "n2"->d[[1]][dCT["V_N2"]]+d[[1]][dCT["bias"]],
(*Helium Threshhold *) "HeT"->d[[1]][dCT["V_N2"]]+d[[1]][dCT["bias"]]+23,
				       "swp"->d[[1]][dCT["V_sw"]],
				       "V_2A"->d[[1]][dCT["V_2A"]],
				       "V_1D"->d[[1]][dCT["V_1D"]]|>;
potentialLines={
Line[{
{potentials[dCT["bias"]],Log[logPlotZero]},
{potentials[dCT["bias"]],3000}
}],
Line[{
{potentials["n2"],Log[logPlotZero]},
{potentials["n2"],3000}
}],
Line[{
{potentials["HeT"],Log[logPlotZero]},
{potentials["HeT"],3000}
}],
Line[{
{potentials["HeT"]+potentials["V_2A"],Log[logPlotZero]},
{potentials["HeT"]+potentials["V_2A"],3000}
}],
Line[{
{potentials["HeT"]-potentials["V_1D"],Log[logPlotZero]},
{potentials["HeT"]-potentials["V_1D"],3000}
}]
};

	(* It would be really nice to plot the bottom axis with the total He offset, and the top axis with the secondary electron energy. What follows is attempting to accomplish that *)

	vPR=Automatic;(*Vertical Plot Range *)
	dataPoints=d[All,{primaryAbcissa,#["CountRate"]&}];
	exfnPlotRaw=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{Automatic,vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
		FrameLabel->{{"Counts",Automatic},{"Total He Offset (V)",Automatic}},
			Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentialLines[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentialLines[[2]],
		colorBlindPalette[[3]],
		potentialLines[[3]],
		Thickness[sweepLineThickness],
		potentialLines[[4]],
		potentialLines[[5]]
		}
	];
	exfnPlotRawNoLabel=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{Automatic,vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
			Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentialLines[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentialLines[[2]],
		colorBlindPalette[[3]],
		potentialLines[[3]],
		Thickness[sweepLineThickness],
		potentialLines[[4]],
		potentialLines[[5]]
		}
	];
	<|"exFnPlotRaw"->exfnPlotRaw,"exFnPlotRawNoLabel"->exfnPlotRawNoLabel,"exFnRawData"->Normal[dataPoints]|>
];

PlotRFAFromDataset[dataset_,
		primaryAbcissa_:dCT["V_he"],
		secondaryAbcissa_:dCT["see"]]:=
Module[{rawImport,d,h,divisions,mmSA,vPR,maxCurrent,
mmPA,primaryGoodDiv,secondaryGoodDiv,dataPoints,
topAxis,plot,plotNoLabel,exp,potentials},
d=dataset;

potentials={
Line[{
{d[[1]][dCT["bias"]],Log[logPlotZero]},
{d[[1]][dCT["bias"]],3000}
}],
Line[{
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]],Log[logPlotZero]},
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]],3000}
}],
Line[{
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]+d[[1]][dCT["V_2A"]],Log[logPlotZero]},
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]+d[[1]][dCT["V_2A"]],3000}
}],
Line[{
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]-d[[1]][dCT["V_1D"]],Log[logPlotZero]},
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]-d[[1]][dCT["V_1D"]],3000}
}],
Line[{
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]+23,Log[logPlotZero]},
{d[[1]][dCT["bias"]]+d[[1]][dCT["V_N2"]]+23,3000}
}]
};
	(* It would be really nice to plot the bottom axis with the total He offset, and the top axis with the secondary electron energy. What follows is attempting to accomplish that *)

vPR={.1,30}; (*Vertical Plot Range *)
dataPoints=d[All,{primaryAbcissa,(#[dCT["I_f"]])/nA &}];
maxCurrent=Max[Transpose[dataPoints][[2]]];
	plot=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{Automatic,vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
		FrameLabel->{{"Current (nA)",Automatic},{"Total He Offset (V)",Automatic}},
		Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentials[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentials[[2]],
		Thickness[sweepLineThickness],
		potentials[[3]],
		potentials[[4]]
		}
	];
	plotNoLabel=ListLogPlot[dataPoints,
		Frame->True,
		GridLines->Automatic,
		PlotRange->{Automatic,vPR},
		FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}},
		Epilog->{Thickness[potentialsLineThickness],
		Dashed,potentials[[1]],
		colorBlindPalette[[2]],
		DotDashed,potentials[[2]],
		Thickness[sweepLineThickness],
		potentials[[3]],
		potentials[[4]]
		}
	];
	<|"rfaPlot"->plot,"rfaPlotNoLabel"->plotNoLabel,"rfaData"->Normal[dataPoints],"maxCurrent"->maxCurrent|>
];




(* ::Chapter:: *)
(*Related Functions*)


ProcessMonitorCurrent[fileName_]:=Module[
{f,h,d,exp,iavg,dwell,
istd,istdLong,noiseLevel,
currentPlot,countPlot,currentAndCountPlot,
results},

f=ImportFile[fileName];
h=f[[1]];
d=f[[2]];


exp=h["MagnitudeOfCurrent(*10^-X)"];
dwell=h["DWELL"];
iavg=d[Mean,"Current"]*10^(-exp);
istd=d[Mean,"CurrentStDev"]*10^(-exp);
istdLong=d[StandardDeviation,"Current"]*10^(-exp);
noiseLevel=istd/iavg;
currentPlot=ListLogPlot[Normal[d[All,"Current"]]*10^(-exp)/nA,PlotRange->{1*^-12/nA,1*^-6/nA},FrameLabel->{"Measurement #","Current (nA)"}];
countPlot=ListLogPlot[Normal[d[All,"Count"]]/dwell,PlotRange->{1,10000},FrameLabel->{"Measurement #","CountRate (Hz)"}];
currentAndCountPlot=TwoAxisListPlot[{Normal[d[All,"Current"]]*10^(-exp),Normal[d[All,"Count"]]/dwell}];
results=<|"AverageCurrent"->iavg,"StdCurrentLong"->istdLong,"StdCurrent"->istd,"NoiseLevel"->noiseLevel,"currentPlot"->currentPlot,"countPlot"->countPlot,"combinedPlot"->currentAndCountPlot|>;
AppendTo[results,h]
];
(*SetAttributes[ProcessMonitorCurrent,Listable];*)


(* ::Chapter:: *)
(*Eventually go through and remove these legacy functions*)


(* Take input of folder name, output stitched data and graphs of excitation functions. *)
CombineExcitationFunctions[files_]:=Module[{constantsFromDataToInclude,constantsFromHeaderToInclude,i,header, magnitude, dat, newDat, allDat,currentData,countData,currentPlot,countPlot,output,sweep,filBias},
output=<||>;
For[i=1,i<=Length[files],i++,
header=GetFileHeaderInfo[files[[i]]];
magnitude=header["MagnitudeOfCurrent(*10^-X)"];
dat=GetFileDataset[files[[i]]];

If[i==1,allDat=newDat,allDat=Join[allDat,newDat]];
];
(*constantsFromHeaderToInclude={"CVGauge(N2)(Torr)","CurrTemp(Res)","CurrTemp(Targ)"};
For[i=1,i<=Length[constantsFromHeaderToInclude],i++,
AppendTo[output,constantsFromHeaderToInclude[[i]]->header[constantsFromHeaderToInclude[[i]]]];
];*)
AppendTo[output,header];
constantsFromDataToInclude={dCT["V_N2"],dCT["V_sw"],dCT["bias"]};
For[i=1,i<=Length[constantsFromDataToInclude],i++,
AppendTo[output,constantsFromDataToInclude[[i]]->dat[All,constantsFromDataToInclude[[i]]][[1]]];
];
sweep=dat[All,dCT["V_sw"]][[1]];
filBias=dat[All,dCT["bias"]][[1]];
currentData=allDat[All,{dCT["V_he"],"Current(nA)"}];
AppendTo[output,"currentData"->currentData];
countData=allDat[All,{dCT["V_he"],"CountRate"}];
AppendTo[output,"countData"->countData];
currentPlot=ListLogPlot[allDat[All,{dCT["V_he"],"Current(nA)"}],PlotRange->Full,GridLines->Automatic,Frame->True,Ticks->Automatic];
AppendTo[output,"currentPlot"->currentPlot];
countPlot=ListLogPlot[allDat[All,{dCT["V_he"],"CountRate"}]];
AppendTo[output,"countPlot"->countPlot];
output
];



