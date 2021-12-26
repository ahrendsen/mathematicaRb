(* ::Package:: *)

BeginPackage["asymmetry`"];
ASYMImport::usage="Imports the file and header into one dataset object"
ASYMAverageCurrent::usage="Reports the average current from all measurements"
ASYMAsymmetry::usage="ASYMAsymmetry[filename,detector] Returns a list of the calculated asymmetries between successive blocks."
ASYMAsymmetry2::usage="ASYMAsymmetry2[filename,detector] Returns a list of the calculated asymmetries between successive blocks, skipping the first half of the first block"
ASYMAsymmetry3::usage="ASYMAsymmetry3[filename,detector] Returns a list of the calculated asymmetries between successive blocks, averaging methods 1 and 2"
ASYMAsymmetryFiltered::usage="ASYMAsymmetry[filename,detector] Returns a list of the calculated asymmetries between successive blocks."
ASYMAsymmetry2Filtered::usage="ASYMAsymmetry2[filename,detector] Returns a list of the calculated asymmetries between successive blocks, skipping the first half of the first block"
ASYMAsymmetry3Filtered::usage="ASYMAsymmetry3[filename,detector] Returns a list of the calculated asymmetries between successive blocks, averaging methods 1 and 2"
ASYMRawOverview::usage="ASYMRawOverview[filename,detector,plotLabel] Returns Plot of all datapoints in order they were collected."
ASYMRawOverviewMethod2::usage="ASYMRawOverviewMethod2[filename,detector,plotLabel] Returns Plot of all datapoints in order they were collected."
ASYMRawOverviewFiltered::usage="ASYMRawOverviewFiltered[filename,detector,plotLabel] Returns Plot of all datapoints in order they were collected."
ASYMRawOverviewOtherPos::usage="ASYMRawOverviewOtherPos[fileName_,detector_] Returns plot at two positions data was collected at."
ASYMBlockAverage::usage="ASYMBlockAverage[filename,detector,plotLabel] Returns Plot of average of blocks of data"
ASYMBlockAverageShift::usage="ASYMBlockAverageShift[filename,detector,plotLabel] Returns Plot of data and data shifted so opposite pump collected first."
ASYMAsymmetrySingleValue::usage="ASYMAymmetrySingleValue[filename,detector]"
ASYMAsymmetrySingleValueFiltered::usage="ASYMAymmetrySingleValue[filename,detector]"


Begin["`Private`"];
<<dataManipulation`; (* For SDM *)
<<fileManipulation`; (* For importing files and headers *)

detColumns={"K617_1","PUMP_1","PROBE_1","REF_1"};
detColumns={"fd","ct","cb","he"};
qwpColumn="qwppos";
sPlusPos=164;
sMinusPos=72;
(* These values must be uncommented to run some historical sets of data 
sPlusPos=160;
sMinusPos=72;
*)
nA=1*^-9;

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
d
},
Switch[ToLowerCase[detector],
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

ASYMRawOverviewMethod2[fileName_,detector_:"cb",plotLabel_:""]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
firstPump,
d
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
firstPump=sPlusPos(*d["data"][[1]][qwpColumn]*);
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
(*
spinUp=Partition[spinUp,d["MeasurementsPerCycle"]];
spinDown=Partition[spinDown,d["MeasurementsPerCycle"]];
*)
spinUp=Partition[spinUp,10];
spinDown=Partition[spinDown,10];
If[firstPump==sPlusPos,
spinUp=Drop[spinUp,1];
spinDown=Drop[spinDown,-1];,
spinUp=Drop[spinUp,-1];
spinDown=Drop[spinDown,1];
];
spinUp=Flatten[spinUp,1];
spinDown=Flatten[spinDown,1];
ListPlot[{{{Null,Null}},spinUp,spinDown},
AspectRatio->1/4,
PlotRange->Full,
ImageSize->Large,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]<>plotLabel]
]

ASYMRawOverviewOtherPos[fileName_,detector_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
qwpPos,
d
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
qwpPos=DeleteDuplicates[Normal[d["data"][[All,qwpColumn]]]];
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==qwpPos[[1]]&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==qwpPos[[2]]&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
ListPlot[{spinUp,spinDown},
AspectRatio->1/4,
PlotRange->Full,
ImageSize->Large,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotLabel->FileBaseName[fileName]]
]


ASYMBlockAverage[fileName_,detector_:"cb",plotLabel_:""]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
t,d,meas,append,header,dataset,
nA=1*^-9
},
Switch[ToLowerCase[detector],
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


ASYMBlockAverageShift[fileName_,detector_,plotLabel_:""]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
shiftUp,shiftDown,
spinDownShift,spinUpShift,
setSpacing=3,
meas,append,t,d,
nA=1*^-9
},
Switch[ToLowerCase[detector],
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
PlotLabel->FileBaseName[fileName]<>plotLabel
,
PlotLegends->Placed[{None,"QWP="<>ToString[sPlusPos],"QWP="<>ToString[sMinusPos]},{.5,.5}]
]
]


Clear[ASYMAsymmetry];
ASYMAsymmetry[fileName_,detector_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn,
meas,append,t,d,
header,dataset,
nA=1*^-9
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinUp=Partition[spinUp,d["MeasurementsPerCycle"]];
spinDown=Partition[spinDown,d["MeasurementsPerCycle"]];
spinUp=Map[Around[Mean[#],SDM[#]]&,spinUp];
spinDown=Map[Around[Mean[#],SDM[#]]&,spinDown];
asymmetry=(spinUp-spinDown)/(spinUp+spinDown)
]


Clear[ASYMAsymmetry2];
ASYMAsymmetry2[fileName_,detector_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn,
header,dataset,meas,append,t,d,
nA=1*^-9,
firstPump
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
firstPump=d["data"][[1]][qwpColumn];
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinUp=Partition[spinUp,d["MeasurementsPerCycle"]];
spinDown=Partition[spinDown,d["MeasurementsPerCycle"]];
If[firstPump==sPlusPos,
spinUp=Drop[spinUp,1];
spinDown=Drop[spinDown,-1];,
spinUp=Drop[spinUp,-1];
spinDown=Drop[spinDown,1];
];
spinUp=Map[Around[Mean[#],SDM[#]]&,spinUp];
spinDown=Map[Around[Mean[#],SDM[#]]&,spinDown];
asymmetry=(spinUp-spinDown)/(spinUp+spinDown)
]


Clear[ASYMAsymmetry3];
ASYMAsymmetry3[fileName_,detector_]:=Module[
{
},
(Drop[ASYMAsymmetry[fileName,detector],1]+ASYMAsymmetry2[fileName,detector])/2
]


Clear[ASYMAsymmetrySingleValue];
ASYMAsymmetrySingleValue[fileName_,detector_]:=Module[
{
asymmetry
},
asymmetry=ASYMAsymmetry3[fileName,detector];
Around[Mean[asymmetry]["Value"],SDM[asymmetry]]
]


Clear[ASYMAverageCurrent];
ASYMAverageCurrent[fileName_, detector_] := Module[{justData, logNumber, spinUp, spinDown, asymmetry, detectorColumn, meas, append, t, d, header, dataset, nA = 1*^-9},
    Switch[ToLowerCase[detector],
        "fd",
            detectorColumn = detColumns[[1]]
        ,
        "ct",
            detectorColumn = detColumns[[2]]
        ,
        "cb",
            detectorColumn = detColumns[[3]]
        ,
        "he",
            detectorColumn = detColumns[[4]]
    ]; 
    d = ASYMImport[fileName]; 
    d = ASYMBlockLevelFilter[d, False]; 
    currents = Dataset[d["data"]][All, #[detectorColumn]&] // Normal; 
    currents = Around[Mean[currents], StandardDeviation[currents]]
]


(* ::Chapter:: *)
(*Filtering out Jumps in Data*)


CooksMask[list_]:=Module[{meansOneDropped,averageDeviation,smallUncertainty},
meansOneDropped=Map[Around[Drop[list,{#}]]&,Range[Length[list]]];
averageDeviation=Mean[Map[#["Uncertainty"]&,meansOneDropped]];
(* Gives a list where each element is True or False, representing if the
error bar for that datapoint is less than half the size of the average 
error bar. *)
smallUncertainty =Map[averageDeviation > 2*#["Uncertainty"]&,meansOneDropped];
(* Collapses the list into a single Boolean value of whether or not
we should reject the block of data.
*)
AnyTrue[smallUncertainty,#==True&]
];


Clear[ASYMBlockLevelFilter]
(* Takes a dataset from the asymmetry program and rejects "blocks" of data
based on whether there are "jumps" in the current *)
ASYMBlockLevelFilter[d_,rejects_]:=Module[
{blockedUp,blockedDown,
upRejectMask,downRejectMask,overallMask,
pickedUp,pickedDown,picked,d2,
vals,as,normAs,rejectMask,rejectSTD},
d["Cycles"];
d["MeasurementsPerCycle"];
(* Partition the dataset into the "blocks" they were collected in*)
blockedUp=Partition[Dataset[d["data"]][Select[#["qwppos"]==sPlusPos&]]//Normal,d["MeasurementsPerCycle"]];
blockedDown=Partition[Dataset[d["data"]][Select[#["qwppos"]==sMinusPos&]]//Normal,d["MeasurementsPerCycle"]];
(* Get just the current values for all measurements *)
vals=Dataset[d["data"]][All,ToLowerCase[d["detector"]]]//Normal;
(* Convert the currents into a "moving" asymmetry *)
as=Differences[vals]/(MovingAverage[vals,2]*2);
(* Normalize the asymmetry to the standard deviation *)
normAs=(as-Mean[as])/StandardDeviation[as];
(* Taking Differences or a moving average reduces the number of data points by 1,
add a zero to the end of the dataset so that we can partition it into equal groups *)
normAs=Prepend[normAs,0];
(* The number of standard deviations a datapoint must be from the 
average asymmetry in order to reject it. *)
rejectSTD=5;
rejectMask=Map[Abs[#]>rejectSTD&,normAs];
overallMask=Partition[rejectMask,2*d["MeasurementsPerCycle"]];
overallMask=Map[AnyTrue[#,TrueQ]&,overallMask];

(* Take the logical Or between the rejection Masks so that if the Up block is rejected, the down is as well. *)
pickedUp=Pick[blockedUp,overallMask,rejects];
pickedDown=Pick[blockedDown,overallMask,rejects];
picked=Flatten[Join[pickedUp,pickedDown]];
(* we picked up the dataset from the last element of d, so we put back the filtered data in the same spot *)
Append[Drop[d,-1],"data"->picked]
];


Clear[ASYMAsymmetryFiltered];
ASYMAsymmetryFiltered[fileName_,detector_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
detectorColumn,
meas,append,t,d,
header,dataset,
nA=1*^-9
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
d=ASYMBlockLevelFilter[d,False];
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinUp=Partition[spinUp,d["MeasurementsPerCycle"]];
spinDown=Partition[spinDown,d["MeasurementsPerCycle"]];
spinUp=Map[Around[Mean[#],SDM[#]]&,spinUp];
spinDown=Map[Around[Mean[#],SDM[#]]&,spinDown];
asymmetry=(spinUp-spinDown)/(spinUp+spinDown)
]


Clear[ASYMAsymmetry2Filtered];
ASYMAsymmetry2Filtered[fileName_,detector_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn,
header,dataset,meas,append,t,d,
nA=1*^-9,
firstPump
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
firstPump=d["data"][[1]][qwpColumn];
d=ASYMBlockLevelFilter[d,False];
spinUp=Dataset[d["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinDown=Dataset[d["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,#[detectorColumn]/nA&]//Normal;
spinUp=Partition[spinUp,d["MeasurementsPerCycle"]];
spinDown=Partition[spinDown,d["MeasurementsPerCycle"]];
If[firstPump==sPlusPos,
spinUp=Drop[spinUp,1];
spinDown=Drop[spinDown,-1];,
spinUp=Drop[spinUp,-1];
spinDown=Drop[spinDown,1];
];
spinUp=Map[Around[Mean[#],SDM[#]]&,spinUp];
spinDown=Map[Around[Mean[#],SDM[#]]&,spinDown];
asymmetry=(spinUp-spinDown)/(spinUp+spinDown)
]


Clear[ASYMAsymmetry3Filtered];
ASYMAsymmetry3Filtered[fileName_,detector_]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
asymmetry,
firstDataLine=11,
detectorColumn
},
(Drop[ASYMAsymmetryFiltered[fileName,detector],1]+ASYMAsymmetry2Filtered[fileName,detector])/2
]


ASYMRawOverviewFiltered[fileName_,detector_,plotLabel_:""]:=Module[
{justData,
logNumber,
spinUp,
spinDown,
firstDataLine=11,
detectorColumn,
nA=1*^-9,
d,d1,d2,spinUpRej,spinDownRej
},
Switch[ToLowerCase[detector],
"fd",detectorColumn=detColumns[[1]],
"ct",detectorColumn=detColumns[[2]],
"cb",detectorColumn=detColumns[[3]],
"he",detectorColumn=detColumns[[4]]
];
d=ASYMImport[fileName];
d1=ASYMBlockLevelFilter[d,False]; (* Pick Those that were accepted *)
d2=ASYMBlockLevelFilter[d,True]; (* Pick those that were rejected *)
spinUp=Dataset[d1["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinDown=Dataset[d1["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinUpRej=Dataset[d2["data"]][Select[#[qwpColumn] ==sPlusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinDownRej=Dataset[d2["data"]][Select[#[qwpColumn] ==sMinusPos&]][All,{#["measurementNumber"],#[detectorColumn]}&]//Normal;
spinUpRej=If[Length[spinUpRej]>0,spinUpRej,{Null,Null}];
spinDownRej=If[Length[spinDownRej]>0,spinDownRej,{Null,Null}];

ListPlot[{{Null,Null},spinUp,spinDown,spinUpRej,spinDownRej},
AspectRatio->1/4,
PlotRange->Full,
ImageSize->Large,
FrameLabel->{"Measurement Number","Current (nA)"},
PlotStyle->
{Global`colorBlindPalette[[1]],
Global`colorBlindPalette[[2]],
Global`colorBlindPalette[[3]],
Blend[{Global`colorBlindPalette[[2]],White},.75],
Blend[{Global`colorBlindPalette[[3]],White},.75]
},
PlotLabel->FileBaseName[fileName]<>plotLabel]
]


Clear[ASYMAsymmetrySingleValueFiltered];
ASYMAsymmetrySingleValueFiltered[fileName_,detector_]:=Module[
{
asymmetry
},
asymmetry=ASYMAsymmetry3Filtered[fileName,detector];
Around[Mean[asymmetry]["Value"],SDM[asymmetry]]
]


End[];
EndPackage[];
