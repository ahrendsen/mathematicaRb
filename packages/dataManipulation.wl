(* ::Package:: *)

(* ::Title:: *)
(*Data Manipulation*)


Clear[DFT];
DFT[list_List,rev_:1]:=Module[{numItems,numItemsSin,numItemsCos,i,k,j,reconstructionCos,reconstructionSin,returnCos,errorsCos,returnSin,errorsSin,intensity},
returnCos=<||>;
errorsCos=<||>;
returnSin=<||>;
errorsSin=<||>;
reconstructionCos={};
reconstructionSin={};
numItems=Length[list];

For[j=0,j<=numItems/2,j++,
AppendTo[returnCos,j->0];
AppendTo[returnSin,j->0];

];

For[i=1,i<=numItems,i++,
intensity=list[[i]];
For [k=0,k<=numItems/2,k++,
(* This IF statement plays the role of the delta function in Berry's Paper *)
If[k==0 ||k==(numItems/2),
(* If K is 0, or L, one of the delta functions is 1, so we divide by 2 and the factor of 2 is not present *)
returnCos[k]=N[returnCos[k]+1/numItems*intensity*Cos[2*\[Pi] *k *rev* (i-1)/numItems]];
returnSin[k]=N[returnSin[k]+1/numItems*intensity*Sin[2*\[Pi] *k *rev* (i-1)/numItems]];
,
(* If K is not 0 or L, there is a factor of 2 in the numerator. *)
returnCos[k]=N[returnCos[k]+2/numItems*intensity*Cos[2*\[Pi] *k* rev *(i-1)/numItems]];
returnSin[k]=N[returnSin[k]+2/numItems*intensity*Sin[2*\[Pi] *k* rev *(i-1)/numItems]];
];
];
];
<|"Sin Coefficients"->returnSin,"Cos Coefficients"->returnCos|>
];

SDM[list_]:=N[StandardDeviation[list]/Sqrt[Length[list]]];


(* ::Chapter:: *)
(*Data Organization*)


CSVtoDataset[filename_]:=Module[{im,dat,h,d},
im=Import[filename];
h=im[[1]];
d=im[[2;;]];
dat=Map[AssociationThread[h,#]&,d];
Dataset[dat]
];
SetAttributes[CSVtoDataset,Listable]

XLSXtoDataset[filename_]:=Module[{im,dat,h,d},
im=Import[filename][[1]];
h=im[[1]];
d=im[[2;;]];
dat=Map[AssociationThread[h,#]&,d];
Dataset[dat]
];
SetAttributes[CSVtoDataset,Listable]

ListFileNamesAndComments[fileNames_]:=Module[{fnAndComments,i,f},
	fnAndComments={};
	For[i=1,i<= Length[fileNames],i++,
		f=ImportFile[fileNames[[i]]];
		AppendTo[fnAndComments,<|f[[1]]["File"]->"FileNumber:"<>ToString[i]<>"--"<>f[[1]]["Comments"]|>];
	];
	Dataset[fnAndComments]
]

TakeColumn[data_,columnNumber_]:=Flatten[Drop[Drop[data,None,{1,columnNumber-1}],None,{2,-1}]]

TakeAllDatasetColumnsBut[dataset_,columnNames_]:=Module[{length,names,cn,i},
	If[Length[columnNames]==0,length=1;cn={columnNames};,length=Length[columnNames]];
	names=Normal[Keys[dataset[[1]]]];
	For[i=1,i<=length,i++,
		names=DeleteCases[names,columnNames[[i]]]
	];
	dataset[All,names]
]

(* I used this function to assist in generating tables of data where we have a list of 
multiple values in some columns but not in others. It "pads" the list so that you have 
to lists of equal length that can be transposed. *)
AddBlanksToList[list_List,blanksBetweenEntries_Integer]:=Module[
{sparseList},
sparseList=Riffle[list,{ConstantArray[Null,blanksBetweenEntries]}];
Flatten[AppendTo[sparseList,ConstantArray[Null,blanksBetweenEntries]]]
]

RemoveErrorVariableColumnsFromDataset[dataset_Dataset]:=dataset[All,DeleteCases[Normal[Keys[dataset[[1]]]],t_String/;StringMatchQ[t,"*err"]]]


(* ::Subtitle:: *)
(*Data Manipulation*)


MultiplyColumn[data_,columnNumber_,value_]:=If[columnNumber==1,Transpose[{Transpose[data][[1]]*value,Transpose[data][[2]]}],Transpose[{Transpose[data][[1]],Transpose[data][[2]]*value}]];
MultiplyColumn[data_,columnNumber_,value_]:=If[columnNumber==1,Transpose[{Transpose[data][[1]]*value,Transpose[data][[2]]}],Transpose[{Transpose[data][[1]],Transpose[data][[2]]*value}]];
