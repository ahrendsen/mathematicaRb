(* ::Package:: *)

BeginPackage["colors`"]
Clear[w,x,y]
ConvertMKBrewerLineToColors[name_]:=
Apply[RGBColor,Map[ToExpression,StringCases[raw,StartOfLine~~name~~"-"~~Except["\n"]..~~"-"~~Except["\n"]..~~" = "~~values:Except["\n"]..->"{"~~values~~"}"]]/255,{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","brewer.txt"}]];
StringSplit[raw,"\n"];
brewerNames=DeleteDuplicates[StringCases[raw,w:Except["\n"]..~~"-"~~x:Except["\n"]..~~"-"~~y:Except["\n"]..~~"-"~~Except["\n"]~~" ="->{w,x,y}]];
brewerPaletteNames=DeleteDuplicates[Map[StringReplace[#,x__~~"-"~~___->x]&,brewerNames[[All,1]]]];
stringIdentifiers=Map[StringJoin[#[[1]],"-",#[[2]]]&,brewerNames];
colors=Map[ConvertMKBrewerLineToColors[#]&,stringIdentifiers];
keys=Map[StringRiffle[#[[;;2]],"-"]&,brewerNames];
overview=Flatten[DeleteCases[StringCases[keys,x__~~"-"~~"5"],{}]];
colorMap=AssociationThread[keys,colors];
mapOverview=AssociationMap[colorMap[#]&,overview];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","cividis.txt"}]];
cividisRGB=Partition[ToExpression[StringSplit[raw]],3];
cividisColors=Apply[RGBColor,cividisRGB,{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","plasma-table-float-0256.csv"}]];
plasmaColors=Apply[RGBColor,raw[[2;;,{2,3,4}]],{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","black-body-table-float-0256.csv"}]];
blackBodyColors=Apply[RGBColor,raw[[2;;,{2,3,4}]],{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","viridis-table-float-0256.csv"}]];
viridisColors=Apply[RGBColor,raw[[2;;,{2,3,4}]],{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","inferno-table-float-0256.csv"}]];
infernoColors=Apply[RGBColor,raw[[2;;,{2,3,4}]],{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","smooth-cool-warm-table-float-0256.csv"}]];
smoothCoolWarmColors=Apply[RGBColor,raw[[2;;,{2,3,4}]],{1}];

raw=Import[FileNameJoin[{$PACKAGES,"palettes","colorValuesSpreadsheet.csv"}]];
allPalettes=SequenceSplit[Take[raw,{3,-1}],{{"","","",""},_}];
colorNumbers=Map[Take[#,{1,-1},{1,3}]&,allPalettes];
NumbersToColors[rgbList_List]:=Apply[RGBColor,rgbList/255,{1}];
colors=Map[NumbersToColors,colorNumbers];
strings=Map[ToString,Flatten[Take[raw,{2,-1}]]];
paletteNames=Flatten[StringCases[strings,"#"~~x___->x]];
colorPalettes=AssociationThread[paletteNames,colors];


EndPackage[];
