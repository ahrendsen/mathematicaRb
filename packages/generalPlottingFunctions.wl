(* ::Package:: *)

(* ::Title:: *)
(*Plotting Functions*)


(* ::Chapter:: *)
(*Plot Formatting Settings*)


imageSize=Small;
wongColorNames={"Black","Orange","Light Blue","Teal","Yellow","Royal Blue","Vermillion","Pink"};
wongColors={RGBColor[0/255,0/255,0/255], (* Black *)
					RGBColor[230/255,159/255,0], (*Orange*)
					RGBColor[86/255,180/255,233/255], (*Light Blue*)
						RGBColor[0/255,158/255,115/255],(* Teal *)
					 RGBColor[240/255,228/255,66/255], (* Yellow *)
					RGBColor[0/255,114/255,178/255],(* Royal Blue *)
					 RGBColor[213/255,94/255,0/255], (* Vermillion *)
					 RGBColor[204/255,121/255,167/255] (* Pink *)
};

colorBlindPalette={
RGBColor[0/255,0/255,0/255], (* Black *)
RGBColor[213/255,94/255,0], (*bamboo*)
RGBColor[34/255,113/255,178/255], (*Honolulu blue*)
RGBColor[53/255,155/255,115/255], (* Ocean Green *)
RGBColor[61/255,183/255,233/255],(* Summer sky *)
RGBColor[230/255,159/255,0/255],(* gamboge *)
RGBColor[247/255,72/255,165/255], (* barbie Pink *)
RGBColor[240/255,228/255,66/255] (* Paris Daisy *)
}

plotMarkAssoc=PlotMarkers->{Graphics[{colorBlindPalette[[2]],Polygon[{{1,0},{0,Sqrt[3]},{-1,0}}]}],.05}

fs=12;   
SetOptions[Plot,BaseStyle->{FontSize->fs},Frame->True,ImageSize->imageSize,LabelStyle->fs,GridLines->Automatic];
SetOptions[ListPlot,BaseStyle->{FontSize->fs},Frame-> True,ImageSize->imageSize,LabelStyle->fs,ImageSize->{584,360},GridLines->Automatic];
SetOptions[ListLogPlot,BaseStyle->{FontSize->fs},ImageSize->imageSize,LabelStyle->fs,ImageSize->{584,360},GridLines->Automatic];
SetOptions[BarChart,BaseStyle->{FontSize->fs,PointSize[Large]},ImageSize->imageSize,GridLines->Automatic];

$WORK=FileNameJoin[{$HomeDirectory,"Box Sync","Gay Group","Project - Rb Spin Filter","karl"}];
$DATARUNS=FileNameJoin[{$HomeDirectory,"Box Sync","Gay Group","Project - Rb Spin Filter","karl","DataRuns"}];
$PACKAGES=FileNameJoin[{$HomeDirectory,"Box Sync","Gay Group","Project - Rb Spin Filter","karl","mathematicaRb","packages"}];


(* ::Chapter::Closed:: *)
(*Labels*)


lbl=<|  "ste"      -> {"Rubidium Spin Polarization, \!\(\*SubscriptBox[\(P\), \(Rb\)]\)","Electron Polarization, \!\(\*SubscriptBox[\(P\), \(e\)]\)"},
		"rfa"      ->{"Helium Potential (V)","Current (nA)"},
		"exfn"      ->{"Helium Potential (V)","Counts (Hz/nA)"},
		"exfnRaw"      ->{"Helium Potential (V)","Counts (Hz)"},
		"twoPumps" ->{"S+ Pump","S- Pump"},
		"n_Rb"     ->"\!\(\*SubscriptBox[\(n\), \(Rb\)]\) (x \!\(\*SuperscriptBox[\(10\), \(12\)]\) \!\(\*SuperscriptBox[\(cm\), \(-3\)]\))",
		"rbPol"    -> {"1/\[Delta] (\!\(\*SuperscriptBox[\(GHz\), \(-1\)]\))","\[CapitalDelta]\[Theta]=\!\(\*SubscriptBox[\(\[Theta]\), \(Pumped\)]\)-\!\(\*SubscriptBox[\(\[Theta]\), \(Magnets\\\ On\)]\)(rad)"}|>


(* ::Chapter::Closed:: *)
(*Ticks*)


oneRevolutionTicks={{0,0},{\[Pi]/2,\[Pi]/2}}


(* ::Chapter::Closed:: *)
(*GridLines*)


oneRevolutionGridlines={{0,Pi/2,Pi,3 Pi/2,2 Pi},Automatic}


(* ::Chapter::Closed:: *)
(*Legending*)


(* Often, there is enough empty space in a plot 
for the legend, but Mathematica doesn't have  a
good setup for placing the legend inside the 
boundaries. You must input the following code:

PlotLegends\[Rule]Placed[PointLegend[{"Label 1","Label2"},
								LegendFunction\[Rule]legendFrame],
					{.2,.75}]
					*)
legendFrame[legend_]:=Framed[legend,
						FrameStyle->Black,
						RoundingRadius->0,
						FrameMargins->0,
						Background->White];


(* ::Chapter::Closed:: *)
(*Specific Plot Settings*)


dataPoints=16;
xPlotRange={-\[Pi]/dataPoints,33\[Pi]/dataPoints};
yPlotRange=Automatic;
anglePos=Range[0,31\[Pi]/dataPoints,2\[Pi]/dataPoints];
(* Generates a list of markers commonly used for identifying rotation. *)
tickM=Range[0,2\[Pi],\[Pi]/2];
xLabel="Polarimeter Position (steps)";
yLabel="Intensity (Hz)";
polarimeterPlotOptions={
FrameTicks->{Partition[Riffle[tickM,tickM],2],Automatic,Automatic,Automatic},
FrameLabel->{xLabel,yLabel},
PlotStyle->colorBlindPalette,
PlotRange->{xPlotRange,yPlotRange},
GridLines->oneRevolutionGridlines
};


(*
polarimeterPlotOptions={FrameTicks\[Rule]{Automatic,{{0,\[Pi]/2,\[Pi],3\[Pi]/2,2\[Pi]},None}},
FrameLabel\[Rule]{xLabel,yLabel},PlotStyle\[Rule]colorBlindPalette,
GridLines\[Rule]oneRevolutionGridlines};
*)


(* ::Chapter:: *)
(*Error Plotting Helpers*)


ListEPlot[xValues_,yValues_,errorValues_,plotOptions_:{}]:=Module[{x,y,yWithErr,eb,op,ope},
x=xValues;
y=yValues;
yWithErr=Apply[Around[#1,#2]&,Transpose[{y,errorValues}],{1}];
op=Partition[Riffle[x,yWithErr],2];(*Ordered Pairs*)
ListPlot[op,PlotRange->{{Min[x],Max[x]},{Min[y-errorValues]*.995,Max[y+errorValues]*1.005}},plotOptions]
];


(* ::Chapter::Closed:: *)
(*Linear Plotting Helpers*)


GetFitPlotFromLinearData[ordPairs_]:=Module[{mm,listPlot,fitPlot,lmf},
mm=MinMax[Transpose[ordPairs][[1]]];
listPlot=ListPlot[ordPairs];
lmf=LinearModelFit[ordPairs,b,b];
Print[lmf];
fitPlot=Plot[lmf[x],{x,mm[[1]],mm[[2]]}];
Show[{listPlot,fitPlot}]
];
GetSlopeFromLinearFit[ordPairs_]:=LinearModelFit[ordPairs,x,x]["BestFitParameters"][[2]];
Assoc2OrdPairs[assoc_]:=Transpose[{Keys[assoc],Values[assoc]}];



(* Needs ordered pairs that are Around[] objects. *)
GayLinearFit[orderedPairs_]:=Module[
{xMeasurements,yMeasurements,
xVal,xErr,yVal,yErr,
a,sa,
syf,af,saf},
{xMeasurements,yMeasurements}=Transpose[orderedPairs];
xVal=Map[#[[1]]&,xMeasurements];
xErr=Map[#[[2]]&,xMeasurements];
yVal=Map[#[[1]]&,yMeasurements];
yErr=Map[#[[2]]&,yMeasurements];
a=Total[(xVal*yVal)/(yErr)^2]/Total[xVal^2/yErr^2];
sa=1/Sqrt[Total[xVal^2/yErr^2]];
syf=Sqrt[(a*xErr)^2+yErr^2];
af=Total[(xVal*yVal)/(syf)^2]/Total[xVal^2/syf^2];
saf=1/Sqrt[Total[xVal^2/syf^2]];
Around[af,saf]
];

Clear[GayLinearFitNoError];
GayLinearFitNoError[orderedPairs_]:=Module[
{xMeasurements,yMeasurements,
xVal,xErr,yVal,yErr,
a,sa,
syf,af,saf},
{xMeasurements,yMeasurements}=Transpose[orderedPairs];
xVal=xMeasurements;
yVal=yMeasurements;
yErr=ConstantArray[.001,Length[yVal]];
xErr=ConstantArray[.001,Length[xVal]];
a=Total[(xVal*yVal)/(yErr)^2]/Total[xVal^2/yErr^2];
sa=Sqrt[(1/Sqrt[Total[xVal^2/yErr^2]])^2+Sqrt[Total[(a*xVal-yVal)^2/(a*xVal)^0]]];
syf=Sqrt[(a*xErr)^2+yErr^2];
af=Total[(xVal*yVal)/(syf)^2]/Total[xVal^2/syf^2];
saf=Sqrt[(1/Sqrt[Total[xVal^2/syf^2]])^2+Sqrt[Total[(af*xVal-yVal)^2/(af*xVal)^0]]];
Around[af,saf]
];


(* ::Chapter::Closed:: *)
(*Two Axis Plot*)


TwoAxisPlot[{f_, g_}, {x_, x1_, x2_}] := 
 Module[{fgraph, ggraph, frange, grange, fticks, 
   gticks}, {fgraph, ggraph} = 
   MapIndexed[
    Plot[#, {x, x1, x2}, Axes -> True, 
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

(* The fourier coefficients that are input should be the output of the DFT I created *)
FourierCoefficientBarChart2[fc_]:=Module[
{bcLabels,blanks,bc,
sinValues,sinValuesNoZero,
cosValues,cosValuesNoZero,
bcLabelsNoZero,bcLabelsSparse,
results,sparser},

results=<||>;
sparser=16;
bcLabels=Take[Keys[fc[[1]]],{2,-1,sparser}];
sinValues=Values[fc[[1]]]/Values[fc[[2]]][[1]];
sinValuesNoZero=Drop[sinValues,{1}];
blanks=ConstantArray[ConstantArray["",sparser-1],Length[bcLabels]*16];
bcLabelsNoZero=Drop[bcLabels,{1}];
bcLabels=Riffle[blanks,bcLabelsNoZero];
bc=BarChart[sinValuesNoZero,ChartStyle->colorBlindPalette[[2]],ChartLabels->Flatten[bcLabels],ImageSize->{720,450},Frame->True];
AppendTo[results,"SinCoefficients"->bc];

bcLabels=Take[Keys[fc[[2]]],{2,-1,sparser}];
cosValues=Values[fc[[2]]]/Values[fc[[2]]][[1]];
cosValuesNoZero=Drop[cosValues,{1}];
bcLabelsNoZero=Drop[bcLabels,{1}];
blanks=ConstantArray[ConstantArray["",sparser],Length[bcLabelsNoZero]];
bcLabelsSparse=Riffle[blanks,bcLabelsNoZero];
bc=BarChart[cosValuesNoZero,ChartStyle->colorBlindPalette[[3]],ChartLabels->Flatten[bcLabelsSparse],ImageSize->{720,450},Frame->True];
AppendTo[results,"CosCoefficients"->bc];
results
];
