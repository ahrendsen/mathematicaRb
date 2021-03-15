(* ::Package:: *)

(* ::Title:: *)
(*Data Manipulation*)


Clear[DFT];
DFT[list_,rev_:1]:=Module[{numItems,i,k,j,reconstructionCos,reconstructionSin,returnCos,errorsCos,returnSin,errorsSin,intensity},
returnCos=<||>;
errorsCos=<||>;
returnSin=<||>;
errorsSin=<||>;
reconstructionCos={};
reconstructionSin={};
numItems=Length[list];

For[j=0,j<numItems/2,j++,
AppendTo[returnCos,j->0];
AppendTo[returnSin,j->0];

];

For[i=1,i<=numItems,i++,
intensity=list[[i]];
For [k=0,k<numItems/2,k++,
If[k==0 ||k==(numItems/2-1),
returnCos[k]=N[returnCos[k]+1/numItems*intensity*Cos[2*\[Pi] *k *rev* (i-1)/numItems]];
returnSin[k]=N[returnSin[k]+1/numItems*intensity*Sin[2*\[Pi] *k *rev* (i-1)/numItems]];,
returnCos[k]=N[returnCos[k]+2/numItems*intensity*Cos[2*\[Pi] *k* rev *(i-1)/numItems]];
returnSin[k]=N[returnSin[k]+2/numItems*intensity*Sin[2*\[Pi] *k* rev *(i-1)/numItems]];
];
];
];
(*

For[j=0,j<numItems/2,j++,
AppendTo[errorsSin,0];
AppendTo[errorsCos,0];
For[i=1,i<=numItems,i++,
intensity=list[[i]];
AppendTo[reconstructionCos,0];
AppendTo[reconstructionSin,0];
reconstructionCos[[i]]=intensity/Cos[2*\[Pi]*j*rev*(i)/numItems] ;
If[j\[Equal]0,reconstructionSin[[i]]=0,reconstructionSin[[i]]=intensity/Sin[2*\[Pi]*1.0001*j*rev*(i)/numItems];];
For [k=0,k<numItems/2,k++,
If[k\[NotEqual]j,
reconstructionCos[[i]]-=(returnCos[[k+1]]*Cos[2*\[Pi]*k*rev*(i)/numItems] +returnSin[[k+1]]*Sin[2*\[Pi]*k*rev*(i)/numItems] )/Cos[2*\[Pi]*j*rev*(i)/numItems];
If[j\[Equal]0,reconstructionSin[[i]]=0,reconstructionSin[[i]]-=(returnCos[[k+1]]*Cos[2*\[Pi]*k*rev*(i)/numItems] +returnSin[[k+1]]*Sin[2*\[Pi]*k*rev*(1)/numItems] )/Sin[2*\[Pi]*1.0001*j*rev*(i)/numItems]];
];
];
errorsSin[[j+1]]+=Power[reconstructionSin[[i]]-returnSin[[j+1]],2];
errorsCos[[j+1]]+=Power[reconstructionCos[[i]]-returnCos[[j+1]],2];
];
errorsSin[[j+1]]=Sqrt[errorsSin[[j+1]]/numItems];
errorsCos[[j+1]]=Sqrt[errorsCos[[j+1]]/numItems];
];
*)
<|"Sin Coefficients"->returnSin,"Cos Coefficients"->returnCos,"Sin Coefficients Error"->errorsSin,"Cos Coefficients Error"->errorsCos|>
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
