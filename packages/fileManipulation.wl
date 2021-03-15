(* ::Package:: *)

(* ::Chapter::Closed:: *)
(*Just like general stuff, or whatever*)


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

Clear[GetFileHeaderInfo];
GetFileHeaderInfo[dataFileName_]:=
Module[{rawFile,data,association,
k,rawImport,removedNulls,
joinedWithSpaces,variableAssociations,as},
rawFile=Import[dataFileName,"tsv"];
data={};
association=Association[];

k=1;
(* 
If the file we're importing is an excel file, then the structure is 

File={Sheet,Sheet,Sheet,...}
Sheet={Row, Row, Row,...}
Row={Col,Col,Col,...}

So I need to select the sheet the data is in. For now, I'll just always
assume the data is in the first sheet 
*)
If[FileExtension[dataFileName]=="xlsx",rawFile=rawFile[[1]]];
While[StringPart[rawFile[[k]][[1]],1]=="#",
If[Length[rawFile[[k]]]> 2,
rawImport=Take[rawFile[[k]],{2,-1}];
removedNulls=DeleteCases[rawImport,Null];
joinedWithSpaces=StringRiffle[removedNulls];
association=Append[association,StringReplace[rawFile[[k]][[1]],{"#"-> "",":"->""}]->joinedWithSpaces];
,
association=Append[association,StringReplace[rawFile[[k]][[1]],{"#"->"",":"->""}]->rawFile[[k]][[2]]];
];
(* This block of code converts the comments to associations. 
It recognizes patterns of <someText>\[Rule]<otherText>, comma
separated, and makes them into associations that are then 
included in the header of the dataset. It's really handy
when you're changing a variable that isn't usually recorded
in the file. 
*)
If[StringMatchQ[Keys[association][[-1]],"Comments"],
(* This makes a list of each association statement*)
variableAssociations=StringCases[association[["Comments"]],RegularExpression["[^ =,]*->[^,]*"]];
(* If there are any associations to add to the dataset *)
If[Length[variableAssociations]>0,
(* Then each association statement is made into a list where the first item is the variable, and the second is the value. *)
variableAssociations=StringSplit[variableAssociations,"->"];
(* Each item is converted to an expression, so numbers will be numbers, not strings of text.*)
(* 2019-05-10: Removing this line and putting it instead only on the value not the key (see next line) *)
(*variableAssociations=ToExpression[variableAssociations];*)
(* Then we need to convert the variable names into strings. We simulataneously make the list of lists into an associations.*)
variableAssociations=AssociationThread[Map[ToString,Transpose[variableAssociations][[1]]]->ToString[Transpose[variableAssociations][[2]]]];
(* We append the associations to the existing collections of associations for the header.*)
AppendTo[association,variableAssociations];
];
];

(* End comments to associations.
*)
k++;
];
association=Append[association,"Time"->GetTimeInfoFromFileNameString[dataFileName]];
association
];

Clear[GetFileDataset]
GetFileDataset[dataFileName_]:=Module[{alldata,dataSections,voltage,dataSection,wavelength,rawFile,
										headerStrippedData,columnHeaderLineNumber,lastColumnToImport,
										j,k,ass,tabularData,ds,rbAbsorptionRef,rbAbsorptionPro,
										frequency},
rawFile=Import[dataFileName,"TSV"];
k=1;
wavelength=0;
dataSections=False;
(* 
If the file we're importing is an excel file, then the structure is 

File={Sheet,Sheet,Sheet,...}
Sheet={Row, Row, Row,...}
Row={Col,Col,Col,...}

So I need to select the sheet the data is in. For now, I'll just always
assume the data is in the first sheet 
*)
If[FileExtension[dataFileName]=="xlsx",rawFile=rawFile[[1]]];
While[StringPart[rawFile[[k]][[1]],1]=="#",
k++;
];
headerStrippedData=Take[rawFile,{k,-1},All];(* Removed the header data *)
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
			(* In case I need to process data with wavelength in the future
			If[dataSections,AppendTo[alldata,c/100/wavelength-\[Nu]0->tabularData];]; *)
			If[dataSections,AppendTo[alldata,frequency-\[Nu]0->tabularData];]; 
			dataSections=True;
			(* In case I need to process data with wavelength in the future 
			wavelength=ToExpression[StringDelete[ StringDelete[headerStrippedData[[j]][[1]],"#VOLT:"~~DigitCharacter..~~"."~~DigitCharacter..~~"("],")"]];*)
			frequency=ToExpression[StringDelete[ StringDelete[headerStrippedData[[j]][[1]],"#VOLT:"~~DigitCharacter..~~"."~~DigitCharacter..~~"("],")"]];
			voltage=ToExpression[StringDelete[ StringDelete[headerStrippedData[[j]][[1]],"#VOLT:"],"("~~DigitCharacter..~~"."~~DigitCharacter..~~")"]];
			(* Then continue adding lines like usual *)
			tabularData={};
		];
	];
];
(* If for some reason I need to process data where I record wavelength again...
If[dataSections,AppendTo[alldata,c/100/wavelength-\[Nu]0->tabularData];,alldata=tabularData];*)
If[dataSections,AppendTo[alldata,frequency-\[Nu]0->tabularData];,alldata=tabularData];
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

GetTimeInfoFromFileNameString[fileNameString_]:=Module[{lastChars,date,year,month,day,hour,min,sec,do},
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

GetTimeStringFromFileNameString[fileNameString_]:=Module[{lastChars,date,time},
lastChars=StringTake[fileNameString,-21];
date=StringTake[lastChars,17]
]

SetAttributes[GetListOfPropertiesFromHeadersOfListOfFileNames,Listable];


(* ::Chapter::Closed:: *)
(*Exporting Data*)


ExportDataToExcel[fileName_,dataset_]:=Module[
{rowNames,colNames,
columnsNamed,colAndRowNamed},
If[Head[Normal[dataset]]==List,
	(* Then we just have column names *)
	rowNames=Transpose[{Join[Normal[Keys[dataset]]]}];
	colNames=Normal[Keys[dataset[[1]]]];
	columnsNamed=Join[{colNames},dataset[All,Values],1];
	Export[fileName,columnsNamed,"XLSX"];
	,
	(* Then we have column and row names *)
		rowNames=Transpose[{Join[{""},Normal[Keys[dataset]]]}];
	colNames=Normal[Keys[dataset[[1]]]];
	columnsNamed=Join[{colNames},dataset[Values,Values],1];
	colAndRowNamed=Join[rowNames,columnsNamed,2];
	Export[fileName,colAndRowNamed,"XLSX"];
];
];


(* ::Chapter::Closed:: *)
(*Moving Files*)


$WORK


raspberryPiIP="129.93.33.86";
$RBRAW=FileNameJoin[{$HomeDirectory,"Box Sync","Gay Group","Project - Rb Spin Filter","RbDataRaw"}];
(* Provide a timestring and the program will transfer over the background files into an appropriate folder *)
MovePolarimeterBeamOffBackground[timeString_String,dateString_String:DateString[DateObject["Today"],{"Year","-","Month","-","Day"}],
								iterations_Integer:3,pumps_List:{"noPump","s+Pump","s-Pump"}]:=
  Module[{rbData,allFileNames,selector,pos,backgroundFiles,pumpOff,pumpOn,beamOffDirectory,laserOffDirectory,laserOnDirectory},
rbData=FileNameJoin[{"\\\\"<>raspberryPiIP<>"\\PiShare","RbData",dateString}];
If[DirectoryQ[rbData],rbData=rbData,rbData=Echo[FileNameJoin[{$RBRAW,dateString}]]];
allFileNames=FileNames["POL*.dat",rbData,2];
selector=Map[StringMatchQ[#,"*"~~timeString~~"*"]&,allFileNames];
pos=Position[allFileNames,Pick[allFileNames,selector][[1]]][[1]][[1]];
backgroundFiles=Take[allFileNames,{pos,pos+iterations*Length[pumps]-1}];
pumpOff=Take[backgroundFiles,{1,-1,Length[pumps]}];
pumpOn=Drop[backgroundFiles,{1,-1,Length[pumps]}];
beamOffDirectory=FileNameJoin[{dataRunFolder,"polarimeterBackground",dateString,"beamOff"}];
laserOffDirectory=FileNameJoin[{beamOffDirectory,"0000"}];
laserOnDirectory=FileNameJoin[{beamOffDirectory,"1800"}];
Map[CreateDirectory[#]&,{laserOffDirectory,laserOnDirectory}];

Map[CopyFile[#,FileNameJoin[{laserOffDirectory,FileNameTake[#]}]]&,pumpOff];
Map[CopyFile[#,FileNameJoin[{laserOnDirectory,FileNameTake[#]}]]&,pumpOn];
  ]
  MovePolarimeterBeamOffBackground::usage=
  "MovePolarimeterBeamOffBackground[] transfers the beam off background data from a run by providing just the timestamp of the first file.
	The timestamp is the only required argument. The other options are: [timestamp, date (if other than today), iterations (if other than standard 3), pumpsList (if other than no, S+,S-)]"
	
	
	(* Provide a timestring and the program will transfer over the background files into an appropriate folder *)
MovePolarimeterGasOffBackground[timeString_String,dateString_String:DateString[DateObject["Today"],{"Year","-","Month","-","Day"}],
								iterations_Integer:5,pumps_List:{"noPump","s+Pump","s-Pump"}]:=
  Module[{rbData,allFileNames,selector,pos,backgroundFiles,pumpOff,pumpOn,beamOffDirectory,laserOffDirectory,laserOnDirectory},
rbData=FileNameJoin[{"\\\\"<>raspberryPiIP<>"\\PiShare","RbData",dateString}];
If[DirectoryQ[rbData],rbData=rbData,rbData=Echo[FileNameJoin[{$RBRAW,dateString}]]];
allFileNames=FileNames["POL*.dat",rbData,2];
selector=Map[StringMatchQ[#,"*"~~timeString~~"*"]&,allFileNames];
pos=Position[allFileNames,Pick[allFileNames,selector][[1]]][[1]][[1]];
backgroundFiles=Take[allFileNames,{pos,pos+iterations*Length[pumps]-1}];
pumpOff=Take[backgroundFiles,{1,-1,Length[pumps]}];
pumpOn=Drop[backgroundFiles,{1,-1,Length[pumps]}];
beamOffDirectory=FileNameJoin[{dataRunFolder,"polarimeterBackground",dateString,"gasOff"}];
laserOffDirectory=FileNameJoin[{beamOffDirectory,"0000"}];
laserOnDirectory=FileNameJoin[{beamOffDirectory,"1800"}];
Map[CreateDirectory[#]&,{laserOffDirectory,laserOnDirectory}];

Map[CopyFile[#,FileNameJoin[{laserOffDirectory,FileNameTake[#]}]]&,pumpOff];
Map[CopyFile[#,FileNameJoin[{laserOnDirectory,FileNameTake[#]}]]&,pumpOn];
  ]
  MovePolarimeterGasOffBackground::usage=
  "MovePolarimeterGasOffBackground[] transfers the beam off background data from a run by providing just the timestamp of the first file.
	The timestamp is the only required argument. The other options are: [timestamp, date (if other than today), iterations (if other than standard 5), pumpsList (if other than no, S+,S-)]"
