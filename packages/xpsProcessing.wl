(* ::Package:: *)

ExtractXPSData[fn_String(*fileName *), sn_String(* Sheet name *)] :=
    Module[{ri, jd, op, lds, dc},
        (* Need to specify that we want the "Data" from the excel sheet, rather than the sheetnames or some other metadata.*)
        ri = Import[fn, {"Data", sn}];(* ri = (R)aw (I)mport *)
        lds = 17;(* The (l)ine the (d)ata (s)tarts on *)
        jd = ri[[lds ;; ]];(* jd = (J)ust (D)ata *)
        dc = {1, 3};(* (d)ata (c)olumns  *)
        op = jd[[All, dc]];(* op= (O)rdered (P)airs *)
       (* <|sn->op|>*)
       op
    ]
SetAttributes[ExtractXPSData, Listable]

Clear[ExtractXPSDataMultiLine];
ExtractXPSDataMultiLine[fn_String(*fileName *), sn_String(* Sheet name *),lines_,normalization_:""] :=
    Module[{ri, jd, op, lds, dc,dcNum,normChars,i},
        (* Need to specify that we want the "Data" from the excel sheet, rather than the sheetnames or some other metadata.*)
        ri = Import[fn, {"Data", sn}];(* ri = (R)aw (I)mport *)
        lds = 20;(* The (l)ine the (d)ata (s)tarts on *)
        jd = ri[[lds ;; ]];(* jd = (J)ust (D)ata *)
dcNum=Range[3,3+lines-1];
        dc = {1, 3, 4};(* (d)ata (c)olumns  *)
        op = jd[[All, {1,#}]]&/@ dcNum;(* op= (O)rdered (P)airs *)
       (* <|sn->op|>*)
       
       normChars=StringSplit[normalization,""];
       For[i=1,i<=Length[normChars],i++,
       Switch[normChars[[i]],
       "l",op=Map[NormalizeLeftSubtract,op];,
       "r",op=Map[NormalizeRightSubtract,op];,
       "s",op=Map[NormalizeSlopeSubtract,op];,
       "p",op=Map[NormalizeToPeak,op];]
       ];
       
       op
    ]


NormalizeRightSubtract[data_]:=Apply[{#1,#2-data[[-1,2]]}&,data,{1}];


NormalizeLeftSubtract[data_]:=Apply[{#1,#2-data[[1,2]]}&,data,{1}];


NormalizeSlopeSubtract[data_] :=
    Module[{rise, run, slope, temp},
        {run, rise} = data[[-1]] - data[[1]];
        slope = Echo[rise / run];
        temp = Apply[{#1, (#2 - slope * #1) / data[[-1, 2]]}&, data, {1}];
        Apply[{#1, #2 - temp[[-1, 2]]}&, temp, {1}]
    ];
ClearAttributes[NormalizeSlopeSubtract, Listable]


NormalizeToPeak[data_]:=Module[
{peak},
peak=Max[Transpose[data][[2]]];
Apply[{#1,#2/Max[peak]}&,data,{1}]
];
ClearAttributes[NormalizeToPeak,Listable]


OffsetDataVertically[
        data_, offset_
    ] :=
    Apply[
        {
            #1, #2 + offset
        }&, data, {
            1
        }
    ];


GetAllYPositionsFromFile[fileName_,numYPositions_,normalizationString_:"r"]:=Module[{sn},
sn=Import[fileName,"Sheets"];
AssociationMap[ExtractXPSDataMultiLine[fileName,#,numYPositions,normalizationString]&, sn[[;;-5]]]
]


GetAllPositionsFromArea[folderName_,normalizationString_:"r"]:=Module[{fn,allData,xpsAll},
SetDirectory[folderName];
fn=FileNames[];
allData=Map[GetAllYPositionsFromFile[#,Length[fn],normalizationString]&,fn];
ResetDirectory[];
xpsAll=Merge[{allData},Catenate];
xpsAll
];


XPSGetListOfMaxes[areaBindingEnergyData_]:=Module[{},
(* Maps the Max function to each of the datasets representing the points in the area scan. 
The part notation at the end ( [[All,All,2]] ) picks out just the intensity data so we don't
get a maximum that is a binding energy. *)
maxes=Map[Max,areaBindingEnergyData[[All,All,2]]]
];


XPSPlotListOfMaxes[listOfMaxes_]:=Module[{},
ListPlot[{listOfMaxes,{{Length[listOfMaxes]+1,Around[Mean[listOfMaxes],StandardDeviation[listOfMaxes]/(Length[listOfMaxes])]}}},PlotRange->Full,FrameLabel->{"Point Number","Intensity"}]
];
