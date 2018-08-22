(* ::Package:: *)

GenerateListOfFileNames[prefix_,date_,runtimes_]:=Module[{j,fileNames},
fileNames={};
For[j=1,j<=Length[runtimes],j++,
AppendTo[fileNames,StringJoin[prefix,date,"_",ToString[runtimes[[j]]],".dat"]]
];
fileNames
];
(* If given a list of files*)
GetIndexOfFilenameFromTimestamp[files_,timestamp_]:=Module[{i,index},
	For[i=1,i<=Length[files],i++,
		If[StringContainsQ[ files[[i]],timestamp],index=i]
	];
	index
];

GetFileHeaderInfo[dataFileName_]:=Module[{rawFile,data,association,k,rawImport,removedNulls,joinedWithSpaces},
rawFile=Import[dataFileName];
data={};
association=Association[];
k=1;
While[StringPart[rawFile[[k]][[1]],1]=="#",
If[Length[rawFile[[k]]]> 2,
rawImport=Take[rawFile[[k]],{2,-1}];
removedNulls=DeleteCases[rawImport,Null];
joinedWithSpaces=StringRiffle[removedNulls];
association=Append[association,StringReplace[rawFile[[k]][[1]],{"#"-> "",":"->""}]->joinedWithSpaces];
,
association=Append[association,StringReplace[rawFile[[k]][[1]],{"#"->"",":"->""}]->rawFile[[k]][[2]]];
];
k++;
];
association
];
Clear[GetFileDataset]
GetFileDataset[dataFileName_]:=Module[{alldata,dataSections,voltage,dataSection,wavelength,rawFile,headerStrippedData,columnHeaderLineNumber,lastColumnToImport,j,k,ass,tabularData,ds,rbAbsorptionRef,rbAbsorptionPro},
rawFile=Import[dataFileName];
k=1;
wavelength=0;
dataSections=False;
While[StringPart[rawFile[[k]][[1]],1]=="#",
k++;
];
headerStrippedData=Take[rawFile,{k,-1},All];(* Removed the header data *)
(*Echo["Removed Header Data"];
Echo[headerStrippedData];*)
alldata=<||>;
tabularData={};
dataSection=<||>;
For[j=2,j<=Length[headerStrippedData],j++,
	ass=Association[];
	If[Length[headerStrippedData[[j]]]>1, (*If there are at least two columns*)
		(*Read the data in like normal*)
		For[k=1,k<=Length[headerStrippedData[[j]]],k++,
			If[dataSections==True,
				ass=Append[ass,"WAV"->ToExpression[wavelength]];
				ass=Append[ass,"VOLT"->ToExpression[voltage]];
			];
			ass=Append[ass,headerStrippedData[[1]][[k]]->headerStrippedData[[j]][[k]]];
		];
		tabularData=Append[tabularData,ass];
		,
		
		(* IF there is one column, it's either a blank line or a data header *)
		If[Length[headerStrippedData[[j]]]>0 ,(* If it has an octothorpe at the beginning of the line *)
			(* It's a header, add it to the dataset *)
			(* But first, record the previous dataset *)
			(* Unless there hasn't been a data section yet *)
			If[dataSections,AppendTo[alldata,c/100/ToExpression[wavelength]-\[Nu]0 ->tabularData];];
			dataSections=True;
			wavelength=StringDelete[ StringDelete[headerStrippedData[[j]][[1]],"#VOLT:"~~DigitCharacter..~~"."~~DigitCharacter..~~"("],")"];
			voltage=StringDelete[ StringDelete[headerStrippedData[[j]][[1]],"#VOLT:"],"("~~DigitCharacter..~~"."~~DigitCharacter..~~")"];
			(* Then continue adding lines like usual *)
			tabularData={};
		];
	];
];
Echo[dataSections];
If[dataSections,AppendTo[alldata,c/100/ToExpression[wavelength]-\[Nu]0 ->tabularData];,alldata=tabularData];
ds=Dataset[alldata]
];

ImportFile[dataFileName_]:={GetFileHeaderInfo[dataFileName],GetFileDataset[dataFileName]}
SetAttributes[GetFileDataset,Listable];
SetAttributes[GetFileHeaderInfo,Listable];
SetAttributes[ImportFile,Listable];


GetListOfPropertiesFromHeadersOfListOfFileNames[fileName_,property_]:=Module[{f},
f=ImportFile[fileName];
Normal[f[[1]][property]]
];

ExtractTimeInfoFromFileNameString[fileNameString_]:=Module[{lastChars,date,year,month,day,hour,min,sec,do},
lastChars=StringTake[fileNameString,-21];
date=StringTake[lastChars,17];
year=StringTake[date,{1,4}];
month=StringTake[date,{6,7}];
day=StringTake[date,{9,10}];
hour=StringTake[date,{12,13}];
min=StringTake[date,{14,15}];
sec=StringTake[date,{16,17}];
do=DateObject[{Interpreter["Number"][year],Interpreter["Number"][month],Interpreter["Number"][day],Interpreter["Number"][hour],Interpreter["Number"][min],Interpreter["Number"][sec]}]
];

SetAttributes[GetListOfPropertiesFromHeadersOfListOfFileNames,Listable];

