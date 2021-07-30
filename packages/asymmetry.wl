(* ::Package:: *)

BeginPackage["asymmetry`"];
ASYMImport::usage="Imports the file and header into one dataset object"
ASYMAsymmetry::usage="ASYMAsymmetry[filename,detector] Returns a list of the calculated asymmetries between successive blocks."
ASYMAsymmetry2::usage="ASYMAsymmetry2[filename,detector] Returns a list of the calculated asymmetries between successive blocks, skipping the first half of the first block"
ASYMAsymmetry3::usage="ASYMAsymmetry3[filename,detector] Returns a list of the calculated asymmetries between successive blocks, averaging methods 1 and 2"
ASYMRawOverview::usage="ASYMRawOverview[filename,detector,plotLabel] Returns Plot of all datapoints in order they were collected."
ASYMBlockAverage::usage="ASYMBlockAverage[filename,detector,plotLabel] Returns Plot of average of blocks of data"
ASYMBlockAverageShift::usage="ASYMBlockAverageShift[filename,detector,plotLabel] Returns Plot of data and data shifted so opposite pump collected first."
ASYMAsymmetrySingleValue::usage="ASYMAymmetrySingleValue[filename,detector]"


Begin["`Private`"];
<<dataManipulation`; (* For SDM *)
<<fileManipulation`; (* For importing files and headers *)

detColumns={"K617_1","PUMP_1","PROBE_1","REF_1"};
detColumns={"fd","ct","cb","he"};
qwpColumn="qwppos";

ASYMImport[fileName_]:=Module[{header,dataset,meas,append,t,d},
{header,dataset}=ImportFile[fileName];
meas=Range[Length[dataset]];
append=Map[<|"measurementNumber"->#|>&,meas];
t=Transpose[{append,dataset//Normal}];
d=Apply[Join[#1,#2]&,t, {1}];
d=Dataset[d][All,<|#,"fd"-> #["fd"]/(1*^-9),"ct"-> #["ct"]/(1*^-9),"cb"-> #["cb"]/(1*^-9),"he"-> #["he"]/(1*^-9)|>&];
Append[header,"data"->d//Normal]
]

ASYMRawOverview[fileName_,detector_:"cb",plotLabel_:""]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72,
nA=1*^-9,
d
},
Switch[detector,
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
ListPlot[{{{Null,Null}},spinUp,spinDown},
AspectRatio->1/4,
PlotRange->Full,
ImageSize->Large,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]<>plotLabel]
]


ASYMBlockAverage[fileName_,detector_:"cb",plotLabel_:""]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72,
t,d,meas,append,header,dataset,
nA=1*^-9
},
Switch[detector,
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
{header,dataset}=ImportFile[fileName];
meas=Range[Length[dataset]];
append=Map[<|"measurementNumber"->#|>&,meas];
t=Transpose[{append,dataset//Normal}];
d=Apply[Join[#1,#2]&,t, {1}];
spinUp=Dataset[d][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;

spinUp=Partition[spinUp,10];
spinDown=Partition[spinDown,10];
spinUp=Map[Around,spinUp];
spinDown=Map[Around,spinDown];
ListPlot[{{{Null,Null}},spinUp,spinDown},
ImageSize->Medium,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]<>plotLabel]
]


ASYMBlockAverageShift[fileName_,detector_:"cb",plotLabel_:""]:=Module[
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
setSpacing=3,
meas,append,t,d,
nA=1*^-9
},
Switch[detector,
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
{header,dataset}=ImportFile[fileName];
meas=Range[Length[dataset]];
append=Map[<|"measurementNumber"->#|>&,meas];
t=Transpose[{append,dataset//Normal}];
d=Apply[Join[#1,#2]&,t, {1}];
spinUp=Dataset[d][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
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
sMinusPos=72,
meas,append,t,d,
header,dataset,
nA=1*^-9
},
Switch[detector,
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
{header,dataset}=ImportFile[fileName];
meas=Range[Length[dataset]];
append=Map[<|"measurementNumber"->#|>&,meas];
t=Transpose[{append,dataset//Normal}];
d=Apply[Join[#1,#2]&,t, {1}];
spinUp=Dataset[d][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
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
sMinusPos=72,
header,dataset,meas,append,t,d,
nA=1*^-9
},
Switch[detector,
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
{header,dataset}=ImportFile[fileName];
meas=Range[Length[dataset]];
append=Map[<|"measurementNumber"->#|>&,meas];
t=Transpose[{append,dataset//Normal}];
d=Apply[Join[#1,#2]&,t, {1}];
spinUp=Dataset[d][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
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
(Drop[ASYMAsymmetry[fileName,detector],1]+ASYMAsymmetry2[fileName,detector])/2
]


Clear[ASYMAsymmetrySingleValue];
ASYMAsymmetrySingleValue[fileName_,detector_:"cb"]:=Module[
{
asymmetry,
firstDataLine=11,
detectorColumn,
sPlusPos=160,
sMinusPos=72
},
asymmetry=ASYMAsymmetry3[fileName,detector];
Mean[asymmetry]
]


End[];
EndPackage[];
