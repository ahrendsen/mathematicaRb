(* ::Package:: *)

BeginPackage["asymmetry`"];
ASYMAsymmetry::usage="ASYMAsymmetry[filename,detector] Returns a list of the calculated asymmetries between successive blocks."
ASYMAsymmetry2::usage="ASYMAsymmetry2[filename,detector] Returns a list of the calculated asymmetries between successive blocks, skipping the first half of the first block"
ASYMAsymmetry3::usage="ASYMAsymmetry3[filename,detector] Returns a list of the calculated asymmetries between successive blocks, averaging methods 1 and 2"
ASYMRawOverview::usage="ASYMRawOverview[filename,detector,plotLabel] Returns Plot of all datapoints in order they were collected."
ASYMBlockAverage::usage="ASYMBlockAverage[filename,detector,plotLabel] Returns Plot of average of blocks of data"
ASYMBlockAverageShift::usage="ASYMBlockAverageShift[filename,detector,plotLabel] Returns Plot of data and data shifted so opposite pump collected first."


Begin["`Private`"];
<<dataManipulation`; (* For SDM *)

ASYMRawOverview[fileName_,detector_:"cb",plotLabel_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72
},
Switch[detector,
"fd",detectorColumn=3,
"ct",detectorColumn=5,
"cb",detectorColumn=7,
"he",detectorColumn=9
];
justData=Import[fileName][[firstDataLine;;]];
logNumber=Range[Length[justData]];
justData=Map[PrependTo[justData[[#]],#]&,logNumber];
spinUp=Flatten[Partition[Cases[justData,{_,sPlusPos,___}][[All,{1,detectorColumn}]],10],1];
spinDown=Flatten[Partition[Cases[justData,{_,sMinusPos,___}][[All,{1,detectorColumn}]],10],1];
ListPlot[{{{Null,Null}},spinUp,spinDown},
AspectRatio->1/4,
PlotRange->Full,
ImageSize->Large,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]<>plotLabel]
]


ASYMBlockAverage[fileName_,detector_:"cb",plotLabel_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72
},
Switch[detector,
"fd",detectorColumn=3,
"ct",detectorColumn=5,
"cb",detectorColumn=7,
"he",detectorColumn=9
];
justData=Import[fileName][[firstDataLine;;]];
logNumber=Range[Length[justData]];
justData=Map[PrependTo[justData[[#]],#]&,logNumber];
spinUp=Cases[justData,{_,sPlusPos,___}][[All,{detectorColumn}]]//Flatten;
spinDown=Cases[justData,{_,sMinusPos,___}][[All,{detectorColumn}]]//Flatten;
spinUp=Partition[spinUp,10];
spinDown=Partition[spinDown,10];
spinUp=Map[Around,spinUp];
spinDown=Map[Around,spinDown];
ListPlot[{{{Null,Null}},spinUp,spinDown},
ImageSize->Medium,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]<>plotLabel]
]


ASYMBlockAverageShift[fileName_,detector_:"cb",plotLabel_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72,
shiftUp,shiftDown,
spinDownShift,spinUpShift,
setSpacing=3
},
Switch[detector,
"fd",detectorColumn=3,
"ct",detectorColumn=5,
"cb",detectorColumn=7,
"he",detectorColumn=9
];
justData=Import[fileName][[firstDataLine;;]];
logNumber=Range[Length[justData]];
justData=Map[PrependTo[justData[[#]],#]&,logNumber];
spinUp=Cases[justData,{_,sPlusPos,___}][[All,{detectorColumn}]]//Flatten;
spinDown=Cases[justData,{_,sMinusPos,___}][[All,{detectorColumn}]]//Flatten;
spinUp=Partition[spinUp,10];
spinDown=Partition[spinDown,10];
spinUp=Map[Around,spinUp];
spinDown=Map[Around,spinDown];
shiftUp=Range[Length[spinUp]]+Length[spinUp]+setSpacing-1;
shiftDown=Range[Length[spinUp]]+Length[spinUp]+setSpacing;
spinUpShift=Transpose[{shiftUp,spinUp}];
spinDownShift=Transpose[{shiftDown,spinDown}];
spinUp=Transpose[{Range[Length[spinUp]],spinUp}];
spinDown=Transpose[{Range[Length[spinUp]],spinDown}];
ListPlot[{{{Null,Null}},Join[spinUp,spinUpShift],Join[spinDown,spinDownShift]},
ImageSize->Medium,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]<>plotLabel,
PlotLegends->Placed[{None,"QWP="<>ToString[sPlusPos],"QWP="<>ToString[sMinusPos]},{.5,.5}]]
]


Clear[ASYMAsymmetry];
ASYMAsymmetry[fileName_,detector_:"cb"]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72
},
Switch[detector,
"fd",detectorColumn=3,
"ct",detectorColumn=5,
"cb",detectorColumn=7,
"he",detectorColumn=9
];
justData=Import[fileName][[firstDataLine;;]];
logNumber=Range[Length[justData]];
justData=Map[PrependTo[justData[[#]],#]&,logNumber];
spinUp=Cases[justData,{_,sPlusPos,___}][[All,{detectorColumn}]]//Flatten;
spinDown=Cases[justData,{_,sMinusPos,___}][[All,{detectorColumn}]]//Flatten;
spinUp=Partition[spinUp,10];
spinDown=Partition[spinDown,10];
spinUp=Map[Around[Mean[#],SDM[#]]&,spinUp];
spinDown=Map[Around[Mean[#],SDM[#]]&,spinDown];
asymmetry=(spinUp-spinDown)/(spinUp+spinDown)
]


Clear[ASYMAsymmetry2];
ASYMAsymmetry2[fileName_,detector_:"cb"]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72
},
Switch[detector,
"fd",detectorColumn=3,
"ct",detectorColumn=5,
"cb",detectorColumn=7,
"he",detectorColumn=9
];
justData=Import[fileName][[firstDataLine;;]];
logNumber=Range[Length[justData]];
justData=Map[PrependTo[justData[[#]],#]&,logNumber];
spinUp=Cases[justData,{_,sPlusPos,___}][[All,{detectorColumn}]]//Flatten;
spinDown=Cases[justData,{_,sMinusPos,___}][[All,{detectorColumn}]]//Flatten;
spinUp=Partition[spinUp,10];
spinDown=Partition[spinDown,10];
spinUp=Drop[spinUp,1];
spinDown=Drop[spinDown,-1];
spinUp=Map[Around[Mean[#],SDM[#]]&,spinUp];
spinDown=Map[Around[Mean[#],SDM[#]]&,spinDown];
asymmetry=(spinUp-spinDown)/(spinUp+spinDown)
]


Clear[ASYMAsymmetry3];
ASYMAsymmetry3[fileName_,detector_:"cb"]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72
},
Switch[detector,
"fd",detectorColumn=3,
"ct",detectorColumn=5,
"cb",detectorColumn=7,
"he",detectorColumn=9
];
(Drop[ASYMAsymmetry[fileName,detector],1]+ASYMAsymmetry2[fileName,detector])/2
]


End[];
EndPackage[];
