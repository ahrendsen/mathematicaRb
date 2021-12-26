(* ::Package:: *)

(** User Mathematica initialization file **)
Needs["PolygonPlotMarkers`"]
size=Offset[20];(* Scaled[7] Would change with the plot size *)
markers=Map[Graphics[{EdgeForm[],PolygonMarker[#,size]},AlignmentPoint->{0,0}]&,{"Circle","UpTriangle","DownTriangle","Square","Diamond","LeftTriangle","RightTriangle","Cross","DiagonalCross","TripleCross","Y"}];
plotMarkers=Partition[Riffle[markers,sizes],2];
fm[name_String, size_ : 6] := Graphics[{EdgeForm[], PolygonMarker[name, Offset[size]]}];
plotMarkers=fm/@{"Circle","UpTriangle","DownTriangle","Square","Diamond","LeftTriangle","RightTriangle","Cross","DiagonalCross","TripleCross","Y"};

imageSize={UpTo[72*3],UpTo[72*2]};

colorBlindPaletteNames={"Black","Bamboo","Honolulu Blue","Ocean Green","Summer Sky","Gamboge","Barbie Pink","Paris Daisy"};
colorBlindPalette={
RGBColor[000/255,000/255,000/255], (* Black         *)
RGBColor[213/255,094/255,000/255], (* Bamboo        *)
RGBColor[034/255,113/255,178/255], (* Honolulu blue *)
RGBColor[053/255,155/255,115/255], (* Ocean Green   *)
RGBColor[061/255,183/255,233/255], (* Summer sky    *)
RGBColor[230/255,159/255,000/255], (* Gamboge       *)
RGBColor[247/255,072/255,165/255], (* Barbie Pink   *)
RGBColor[240/255,228/255,066/255]  (* Paris Daisy   *)
};
qualitativePalette=colorBlindPalette;

tolIridescentStrings={"FEFBE9","FCF7D5","F5F3C1","EAF0B5","DDECBF","DDE7CA","C2E3D2","B5DDD8","A8D8DC","98D2E1","8DCBE4","81C4E7","7BBCE7","7EB2E4","88A5DD","9398D2","9B8AC4","9D7DB2","9A709E","906388","805770","684957","46353A","999999"};
tolIridescent=Map[StringJoin["#",#]&,tolIridescentStrings];
tolIridescentPalette=Map[RGBColor,tolIridescent];
sequentialPalette=tolIridescentPalette;

tolPurpleGreenValues={{118/255,14/85,131/255},{3/5,112/255,57/85},{194/255,11/17,69/85},{77/85,212/255,232/255},{247/255,247/255,247/255},{217/255,16/17,211/255},{172/255,211/255,158/255},{6/17,58/85,97/255},{9/85,8/17,11/51},{1,14/15,3/5}};
tolPurpleGreenPalette=Map[RGBColor,tolPurpleGreenValues];
divergingPalette=tolPurpleGreenPalette;



labelStyle=10;
fontSize=10;	   
SetOptions[Plot,BaseStyle->{FontSize->fontSize},Frame->True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic];
SetOptions[LogPlot,BaseStyle->{FontSize->fontSize},Frame->True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic];

SetOptions[ListPlot,BaseStyle->{FontSize->fontSize},Frame-> True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic,PlotMarkers->plotMarkers,PlotStyle->colorBlindPalette];
SetOptions[ListLinePlot,BaseStyle->{FontSize->fontSize},Frame-> True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic,PlotMarkers->plotMarkers,PlotStyle->colorBlindPalette];
SetOptions[ListLogPlot,BaseStyle->{FontSize->fontSize},Frame-> True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic,PlotMarkers->plotMarkers,PlotStyle->colorBlindPalette];
SetOptions[ListLogLogPlot,BaseStyle->{FontSize->fontSize},Frame-> True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic,PlotMarkers->plotMarkers,PlotStyle->colorBlindPalette];
SetOptions[ListLogLinearPlot,BaseStyle->{FontSize->fontSize},Frame-> True,ImageSize->imageSize,LabelStyle->labelStyle,GridLines->Automatic,PlotMarkers->plotMarkers,PlotStyle->colorBlindPalette];

SetOptions[BarChart,BaseStyle->{FontSize->fontSize,PointSize[Large]},Frame->True,ImageSize->imageSize,GridLines->Automatic];
SetOptions[Histogram,BaseStyle->{FontSize->fontSize,PointSize[Large]},Frame->True,ImageSize->imageSize,GridLines->Automatic];


$WORK=FileNameJoin[{$HomeDirectory,"Box Sync","Gay Group","Project - Rb Spin Filter","karl"}];
$DATARUNS=FileNameJoin[{$WORK,"DataRuns"}];
$PACKAGES=FileNameJoin[{$WORK,"mathematicaRb","packages"}];
$PROCESSEDDATA=FileNameJoin[{$WORK,"ProcessedDataArchive"}];

AppendTo[$Path,$PACKAGES];
