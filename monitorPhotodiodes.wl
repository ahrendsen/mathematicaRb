(* ::Package:: *)

MPDOrganizeMonitorPhotodiodeData[fileName_ ,column_:3]:=Module[{ri,rd,ct,timeList,to,times,dateOp},
(* Raw Import *)
ri=Import[fileName] ;
(* Raw Data *)
rd = ri[[11;;-2]];
(* chiral Target *)
ct = rd[[All,column]]/1*^-9;
timeList=StringTake[fileName,{ {-10,-9},{-8,-7},{-6,-5}}];
timeList=Map[ToExpression,timeList]; (* Need to convert strings to numbers, so TimeObject can interpret them *)
to =TimeObject[timeList];
times=to+Map[Quantity[#,"Second"]&,Range[Length[ct]]*1.5];
dateOp=Partition[Riffle[times,ct],2](* Date ordered paris *)
];

MPDGetFilesBetweenTimestamps[date_,timeString1_,timeString2_]:=Module[{rbData,username,fileNames,fileString1,fileString2,pos1,pos2,day},
username=SystemInformation["Kernel","Username"];
(* If the username is kahrendsen2, we're on the lab computer, so we can access the raspberry pi on the network. Otherwise, use the backup on Box. *)
If[username=="kahrendsen2",
rbData=FileNameJoin[{"\\\\129.93.33.132\\PiShare","RbData"}];,
rbData=FileNameJoin[{$HomeDirectory,"OneDrive - University of Nebraska-Lincoln","Gay Group","Project - Rb Spin Filter","RbDataRaw"}];
];
fileString1=Flatten[FileNames[FileNameJoin[{rbData,date,"*"<>timeString1<>"*.dat"}]]];
fileNames=FileNames[FileNameJoin[{rbData,date,"*.dat"}]];
pos1=Flatten[Position[fileNames,fileString1[[1]]]][[1]];
fileString2=Flatten[FileNames[FileNameJoin[{rbData,date,"*"<>timeString2<>"*.dat"}]]];
pos2=Flatten[Position[fileNames,fileString2[[1]]]][[1]];
fileNames=fileNames[[pos1;;pos2]];
day=StringTake[date,-2];

If[DirectoryQ[day],,CreateDirectory[day]];
Map[CopyFile[#,FileNameJoin[{day,FileNameTake[#]}]]&,fileNames];
Map[FileNameJoin[{day,FileNameTake[#]}]&,fileNames]
];
