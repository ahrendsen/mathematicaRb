(* ::Package:: *)

ExtractXPSData[fn_String(*fileName *), sn_String(* Sheet name *)] :=
    Module[{ri, jd, op, lds, dc},
        (* Need to specify that we want the "Data" from the excel sheet, rather than the sheetnames or some other metadata.*)
        ri = Import[fn, {"Data", sn}];(* ri = (R)aw (I)mport *)
        lds = 17;(* The (l)ine the (d)ata (s)tarts on *)
        jd = ri[[lds ;; ]];(* jd = (J)ust (D)ata *)
        dc = {1, 3};(* (d)ata (c)olumns  *)
        op = jd[[All, dc]] (* op= (O)rdered (P)airs *)
    ]
SetAttributes[ExtractXPSData, Listable]


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
