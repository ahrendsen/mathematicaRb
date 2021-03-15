(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



(* :Title: FigText *)
(* :Context: SciDraw` *)
(* :Author: Mark A. Caprio, Department of Physics, University of Notre Dame *)
(* :Summary: Miscellaneous text formatting utilities. *)
(* :Copyright: Copyright 2011, Mark A. Caprio *)
(* :Package Version: 0.0 *)
(* :Mathematica Version: 7.0 *)
(* :History:
Package MCText started December 2003, based upon functions written for LevelScheme.
Distributed with the LevelScheme package.
October 2010. Context changed to LevelScheme`.
May 2011.  Renamed from MCText to FigText.  Context changed to SciDraw`.
*)


BeginPackage["SciDraw`",SciDraw`Private`$ExternalContexts];


Unprotect[Evaluate[$Context<>"*"]];


Begin["`Private`"];


textup[x_]:=StyleForm[x,FontSlant->"Plain"];
textsl[x_]:=StyleForm[x,FontSlant->"Oblique"];
textit[x_]:=StyleForm[x,FontSlant->"Italic"];
textmd[x_]:=StyleForm[x,FontWeight->"Plain"];
textbf[x_]:=StyleForm[x,FontWeight->"Bold"];
textrm[x_]:=StyleForm[x,FontFamily->"Times"];
texttt[x_]:=StyleForm[x,FontFamily->"Courier"];
textsf[x_]:=StyleForm[x,FontFamily->"Helvetica"];


hspace[Lems_]:=AdjustmentBox["",BoxMargins->{{0,Lems},{0,0}}];


textsize[s_,x_]:=StyleForm[x,FontSize->s];
textcolor[c_,x_]:=StyleForm[x,FontColor->c];
texttracking[t_,x_]:=StyleForm[x,FontTracking->t];
textfamily[f_,x_]:=StyleForm[x,FontFamily->f];
texthidden[x_]:=StyleForm[x,ShowContents->False];


textsubscript[x_]:=SubscriptBox["",x];
textsuperscript[y_]:=SuperscriptBox["",y];
textsubsuperscript[x_,y_]:=SubsuperscriptBox["",x,y];


Options[textit]={hspace->0};
textit[x_,Opts___]:=Module[
{FullOpts=Flatten[{Opts,Options[textit]}]},
AdjustmentBox[
StyleForm[x,FontSlant->"Italic"],
BoxMargins->{{0,hspace/.FullOpts},{0,0}}
]
];


StackText[Alignment_,Spacing_,Lines_List,Opts___?OptionQ]:=GridBox[{#}&/@Lines,ColumnAlignments->Alignment,RowSpacings->Spacing,Opts];


SuperPrimeBox[x_,n:(_Integer?NonNegative):1]:=SuperscriptBox[x,StringJoin[Table["\[Prime]",{n}]]];
(*SuperPrime[x_]:=Superscript[x,"\[Prime]"];*)
SuperPrime[x_,n:(_Integer?NonNegative):1]:=Superscript[x,StringJoin[Table["\[Prime]",{n}]]];


Options[MultipletLabel]={EntrySeparator->",",Delimiter->{"(",")"}};
MultipletLabel[Values_List,Opts:OptionsPattern[]]:=Module[
{DeducedEntrySeparator,DeducedDelimiter},

DeducedDelimiter=Switch[
OptionValue[Delimiter],
None,{"",""},
_,OptionValue[Delimiter]
];
DeducedEntrySeparator=Switch[
OptionValue[EntrySeparator],
None,"",
_,OptionValue[EntrySeparator]
];

Row[Join[{DeducedDelimiter[[1]]},Riffle[Values,DeducedEntrySeparator],{DeducedDelimiter[[2]]}]]
];


SignString[x_?NumericQ]:=Switch[
Sign[x],
+1,"+",
0,"",
-1,"-"
];


SetAttributes[Sqrtize,Listable];
Sqrtize[x_?NumericQ]:=Which[
(*RationalQ[x],x,*)
RationalQ[x^2],Sign[x]*SqrtBox[x^2],
True,x
];
SetAttributes[Radicalize,Listable];
Radicalize[x_?NumericQ,n_Integer]:=Which[
(*RationalQ[x],x,*)
RationalQ[x^n],Sign[x]*RadicalBox[Abs[x^n],n],
True,x
];


RationalQ[x_]:=(IntegerQ[x]||(Head[x]===Rational));


Options[AlignmentBox]={
AlignmentMarker->"&",
Align->True,
ColumnAlignments->{Right,Left},ColumnSpacings->0
};


BreakString[Separator_,Str_]:=Module[
{PosnList},
PosnList=Join[{{Null,0}},StringPosition[Str,Separator],{{StringLength[Str]+1,Null}}];

Table[StringTake[Str,{PosnList[[i]][[2]]+1,PosnList[[i+1]][[1]]-1}],
{i,1,Length[PosnList]-1}
]
];

AlignmentBox[Str_,Opts___?OptionQ]:=Module[
{
FullOpts=Flatten[{Opts,Options[AlignmentBox]}]
},

CheckOption[Align,True|False,FullOpts];
CheckOption[AlignmentMarker,_String,FullOpts];
StyleBox[
If[
Align/.FullOpts,
GridBox[
{BreakString[(AlignmentMarker/.FullOpts),Str]},
Sequence@@FilterRules[FullOpts,Options@GridBox]
],
StringReplace[Str,(AlignmentMarker/.FullOpts)->""]
],
Sequence@@FilterRules[FullOpts,Options@StyleBox]
]
];


LaTeXTableEntryValue[Value_?NumericQ]:=Value;
LaTeXTableEntryValue[Str_String]:=Module[
{
FirstNumericPosn,
ErrorBarsPosn,
Value
},

FirstNumericPosn=StringPosition[Str,{"+","-","0","1","2","3","4","5","6","7","8","9","."},1][[1,1]];
ErrorBarsPosn=If[
Length[StringPosition[Str,"(",1]]>=1,
StringPosition[Str,{"(","\[PlusMinus]"},1][[1,1]],
StringLength[Str]+1
];
Value=ToExpression[
StringReplace[
StringTake[Str,{FirstNumericPosn,ErrorBarsPosn-1}],
{"&"->""}
]
];
If[!NumericQ[Value],
Message[LaTeXTableEntryValue::notnumeric,Str,Value]
];
Value
];


Options[SubarrayEllipsis]={Padding->"\[CenterEllipsis]"};
SubarrayEllipsis[m_?MatrixQ,{Rows:(_Integer|Infinity),Columns:(_Integer|Infinity)},OptionsPattern[]]:=Module[
{
RowTrim,ColumnTrim,RowMax,ColumnMax,MatrixTrimmed
},
RowTrim=(Dimensions[m][[1]]>Rows);
RowMax=Min[Dimensions[m][[1]],Rows];
ColumnTrim=(Dimensions[m][[2]]>Columns);
ColumnMax=Min[Dimensions[m][[2]],Columns];

ArrayPad[
m[[1;;RowMax,1;;ColumnMax]],
{{0,If[RowTrim,1,0]},{0,If[ColumnTrim,1,0]}},
OptionValue[Padding]
]
];
SubarrayEllipsis[m_List,{Rows:(_Integer|Infinity)},OptionsPattern[]]:=Module[
{
RowTrim,RowMax,MatrixTrimmed
},
RowTrim=(Dimensions[m][[1]]>Rows);
RowMax=Min[Dimensions[m][[1]],Rows];

ArrayPad[
m[[1;;RowMax]],
{{0,If[RowTrim,1,0]}},
OptionValue[Padding]
]
];


PageBreak[]:=CellPrint[Cell["",PageBreakBelow->True]];


UnitsLabel[FactorSequence___]:=Module[
{EntryList},

(* create exponents *)
EntryList=Replace[{FactorSequence},{u_,p_}:>Superscript[u,p],{1}];

(* intersperse thin spaces *)
EntryList=Riffle[EntryList,"\[ThinSpace]"(*ThinSpace*)];

(* but undo spaces around solidus *)
EntryList=Replace[EntryList,{pre___,"\[ThinSpace]"(*ThinSpace*),"/","\[ThinSpace]"(*ThinSpace*),post___}:>{pre,"/",post},{0}];

(* create row *)
Row[EntryList]
];


SpectroscopicLetter[0]="s";
SpectroscopicLetter[1]="p";
SpectroscopicLetter[2]="d";
SpectroscopicLetter[3]="f";
SpectroscopicLetter[4]="g";
SpectroscopicLetter[5]="h";
SpectroscopicLetter[6]="i";
SpectroscopicLetter[7]="j";
SpectroscopicLetter[8]="k";


Options[ShellLabel]={Style->Italic};
ShellLabel[{n_,l_,j_},OptionsPattern[]]:=Row[{n,Subscript[Style[SpectroscopicLetter[l],OptionValue[Style]],SolidusFractionize[j]]}];
ShellLabel[{l_,j_},OptionsPattern[]]:=Row[{Subscript[Style[SpectroscopicLetter[l],OptionValue[Style]],SolidusFractionize[j]]}];


Options[NucleusBox]={NuclearA->"",NuclearZ->"",NuclearN->""};
NucleusBox[Element_,Opts___?OptionQ]:=Module[
{FullOpts=Flatten[{Opts,Options[NucleusBox]}]},
Row[{
SubsuperscriptBox["",NuclearZ/.FullOpts,NuclearA/.FullOpts],SubsuperscriptBox[Element,NuclearN/.FullOpts,""]  (* to match subsuperscript on left for alignment *)
}]
];


(* LIMITATION: A and Z are left aligned *)


Isotope[A_:None,Z_:None,N_:None,Sup_:None,Element_String]:=Which[
(Z===None)&&(N===None),
(* use subsuperscript on both sides if use on either, to match alignment *)
Row[{
Superscript["",Replace[A,None->""]],Superscript[Element,Replace[Sup,None->""]]  
}],
True,
Row[{
Subsuperscript["",Replace[Z,None->""],Replace[A,None->""]],Subsuperscript[Element,Replace[N,None->""],Replace[Sup,None->""]]  
}]
];


Isotope[PriorArgs___,(Z_Integer)?NonNegative]:=Isotope[PriorArgs,ElementAbbreviation[Z]];


Isotope[Sup_:None,{Z_Integer,N_Integer}]:=Isotope[Z+N,None,None,Sup,ElementAbbreviation[Z]];


ElementAbbreviations=
	{"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P",
	 "S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
	 "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc",
	 "Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La",
	 "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	 "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
	 "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf",
	 "Es","Fm","Md","No","Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
	"Uun","Uuu","Uub"}; (* data circa Mathematica 5.0 *)
ElementAbbreviationsLength=Length[ElementAbbreviations];


(* from ChemicalElements package 1.4 Elements list, converted to strings *)
ElementNames={"Hydrogen","Helium","Lithium","Beryllium","Boron","Carbon","Nitrogen","Oxygen","Fluorine","Neon","Sodium","Magnesium","Aluminium","Silicon","Phosphorus","Sulfur","Chlorine","Argon","Potassium","Calcium","Scandium","Titanium","Vanadium","Chromium","Manganese","Iron","Cobalt","Nickel","Copper","Zinc","Gallium","Germanium","Arsenic","Selenium","Bromine","Krypton","Rubidium","Strontium","Yttrium","Zirconium","Niobium","Molybdenum","Technetium","Ruthenium","Rhodium","Palladium","Silver","Cadmium","Indium","Tin","Antimony","Tellurium","Iodine","Xenon","Caesium","Barium","Lanthanum","Cerium","Praseodymium","Neodymium","Promethium","Samarium","Europium","Gadolinium","Terbium","Dysprosium","Holmium","Erbium","Thulium","Ytterbium","Lutetium","Hafnium","Tantalum","Tungsten","Rhenium","Osmium","Iridium","Platinum","Gold","Mercury","Thallium","Lead","Bismuth","Polonium","Astatine","Radon","Francium","Radium","Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium","Americium","Curium","Berkelium","Californium","Einsteinium","Fermium","Mendelevium","Nobelium","Lawrencium","Rutherfordium","Dubnium","Seaborgium","Bohrium","Hassium","Meitnerium","Ununnilium","Unununium","Ununbium"};
ElementNamesLength=Length[ElementNames];


StableIsotopes["Hydrogen"]={1,2};
StableIsotopes["Helium"]={3,4};
StableIsotopes["Lithium"]={6,7};
StableIsotopes["Beryllium"]={9};
StableIsotopes["Boron"]={10,11};
StableIsotopes["Carbon"]={12,13};
StableIsotopes["Nitrogen"]={14,15};
StableIsotopes["Oxygen"]={16,17,18};
StableIsotopes["Fluorine"]={19};
StableIsotopes["Neon"]={20,21,22};
StableIsotopes["Sodium"]={23};
StableIsotopes["Magnesium"]={24,25,26};
StableIsotopes["Aluminium"]={27};
StableIsotopes["Silicon"]={28,29,30};
StableIsotopes["Phosphorus"]={31};
StableIsotopes["Sulfur"]={32,33,34,36};
StableIsotopes["Chlorine"]={35,37};
StableIsotopes["Argon"]={36,38,40};
StableIsotopes["Potassium"]={39,41};
StableIsotopes["Calcium"]={40,42,43,44,46,48};
StableIsotopes["Scandium"]={45};
StableIsotopes["Titanium"]={46,47,48,49,50};
StableIsotopes["Vanadium"]={51};
StableIsotopes["Chromium"]={50,52,53,54};
StableIsotopes["Manganese"]={55};
StableIsotopes["Iron"]={54,56,57,58};
StableIsotopes["Cobalt"]={59};
StableIsotopes["Nickel"]={58,60,61,62,64};
StableIsotopes["Copper"]={63,65};
StableIsotopes["Zinc"]={66,67,68,70};
StableIsotopes["Gallium"]={69,71};
StableIsotopes["Germanium"]={70,72,73,74,76};
StableIsotopes["Arsenic"]={75};
StableIsotopes["Selenium"]={74,76,77,78,80,82};
StableIsotopes["Bromine"]={79,81};
StableIsotopes["Krypton"]={78,80,82,83,84,86};
StableIsotopes["Rubidium"]={85};
StableIsotopes["Strontium"]={84,86,87,88};
StableIsotopes["Yttrium"]={89};
StableIsotopes["Zirconium"]={90,91,92,94};
StableIsotopes["Niobium"]={93};
StableIsotopes["Molybdenum"]={92,94,95,96,97,98,100};
StableIsotopes["Technetium"]={};
StableIsotopes["Ruthenium"]={96,98,99,100,101,102,104};
StableIsotopes["Rhodium"]={103};
StableIsotopes["Palladium"]={102,104,105,106,108,110};
StableIsotopes["Silver"]={107,109};
StableIsotopes["Cadmium"]={106,108,110,111,112,113,114,116};
StableIsotopes["Indium"]={113};
StableIsotopes["Tin"]={112,114,115,116,117,118,119,120,122,124};
StableIsotopes["Antimony"]={121,123};
StableIsotopes["Tellurium"]={120,122,124,125,126,128,130};
StableIsotopes["Iodine"]={127};
StableIsotopes["Xenon"]={129,130,131,132,134,136};
StableIsotopes["Caesium"]={133};
StableIsotopes["Barium"]={130,132,134,135,136,137,138};
StableIsotopes["Lanthanum"]={138};
StableIsotopes["Cerium"]={136,138,140,142};
StableIsotopes["Praseodymium"]={141};
StableIsotopes["Neodymium"]={142,143,145,146,148,150};
StableIsotopes["Promethium"]={};
StableIsotopes["Samarium"]={144,150,152,154};
StableIsotopes["Europium"]={151,153};
StableIsotopes["Gadolinium"]={154,155,156,157,158,160};
StableIsotopes["Terbium"]={159};
StableIsotopes["Dysprosium"]={156,158,160,161,162,163,164};
StableIsotopes["Holmium"]={165};
StableIsotopes["Erbium"]={162,164,166,167,168,170};
StableIsotopes["Thulium"]={169};
StableIsotopes["Ytterbium"]={168,170,171,172,173,174,176};
StableIsotopes["Lutetium"]={175};
StableIsotopes["Hafnium"]={176,177,178,179,180};
StableIsotopes["Tantalum"]={181};
StableIsotopes["Tungsten"]={180,182,183,184,186};
StableIsotopes["Rhenium"]={185};
StableIsotopes["Osmium"]={184,186,187,188,189,190,192};
StableIsotopes["Iridium"]={191,193};
StableIsotopes["Platinum"]={194,195,196,198};
StableIsotopes["Gold"]={197};
StableIsotopes["Mercury"]={196,198,199,200,201,202,204};
StableIsotopes["Thallium"]={203,205};
StableIsotopes["Lead"]={204,206,207,208};
StableIsotopes["Bismuth"]={209};
StableIsotopes["Polonium"]={};
StableIsotopes["Astatine"]={};
StableIsotopes["Radon"]={};
StableIsotopes["Francium"]={};
StableIsotopes["Radium"]={};
StableIsotopes["Actinium"]={};
StableIsotopes["Thorium"]={};
StableIsotopes["Protactinium"]={};
StableIsotopes["Uranium"]={};
StableIsotopes["Neptunium"]={};
StableIsotopes["Plutonium"]={};
StableIsotopes["Americium"]={};
StableIsotopes["Curium"]={};
StableIsotopes["Berkelium"]={};
StableIsotopes["Californium"]={};
StableIsotopes["Einsteinium"]={};
StableIsotopes["Fermium"]={};
StableIsotopes["Mendelevium"]={};
StableIsotopes["Nobelium"]={};
StableIsotopes["Lawrencium"]={};
StableIsotopes["Rutherfordium"]={};
StableIsotopes["Dubnium"]={};
StableIsotopes["Seaborgium"]={};
StableIsotopes["Bohrium"]={};
StableIsotopes["Hassium"]={};
StableIsotopes["Meitnerium"]={};
StableIsotopes["Ununnilium"]={};
StableIsotopes["Unununium"]={};
StableIsotopes["Ununbium"]={};



ElementAbbreviation[Z:0]="n";
ElementAbbreviation[(Z_Integer)?Positive]/;(Z<=ElementAbbreviationsLength):=ElementAbbreviations[[Z]];
ElementAbbreviation[(Z_Integer)?Positive]/;(Z>ElementAbbreviationsLength):=ToString[Z];


ElementName[Z:0]="Neutron";
ElementName[(Z_Integer)?Positive]/;(Z<=ElementNamesLength):=ElementNames[[Z]];
ElementName[(Z_Integer)?Positive]/;(Z>ElementNamesLength):=ToString[Z];


(*IsotopeIsStable[{Z:0,(A_Integer)?Positive}]:=False;
IsotopeIsStable[{(Z_Integer)?Positive,(A_Integer)?Positive}]/;(Z\[LessEqual]ElementNamesLength):=MemberQ[StableIsotopes[ElementName[Z]],A];
IsotopeIsStable[{(Z_Integer)?Positive,(A_Integer)?Positive}]/;(Z>ElementNamesLength):=False;*)


StableIsotopes[Z:0]:={};
StableIsotopes[(Z_Integer)?Positive]/;(Z<=ElementNamesLength):=StableIsotopes[ElementName[Z]];
StableIsotopes[(Z_Integer)?Positive]/;(Z>ElementNamesLength):={};


IsotopeIsStable[{Z:0,(N_Integer)?Positive}]:=False;
IsotopeIsStable[{(Z_Integer)?Positive,(N_Integer)?NonNegative}]/;(Z<=ElementNamesLength):=MemberQ[StableIsotopes[ElementName[Z]],Z+N];
IsotopeIsStable[{(Z_Integer)?Positive,(N_Integer)?NonNegative}]/;(Z>ElementNamesLength):=False;


Isotope[Args___,(Z_Integer)?NonNegative]:=Isotope[Args,ElementAbbreviation[Z]];


Options[LabelJiP]={Rational->SolidusFractionize};
LabelJiP[J_,i_,P:(+1|-1|None):+1,OptionsPattern[]]:=SubsuperscriptBox[
J/.{x_Rational:>OptionValue[Rational][x]},
i,
Switch[P,+1,"+",-1,"-",None,""]
];
Options[LabelJP]={Rational->SolidusFractionize};
LabelJP[J_,P:(+1|-1|None):+1,OptionsPattern[]]:=SuperscriptBox[
J/.{x_Rational:>OptionValue[Rational][x]},
Switch[P,+1,"+",-1,"-",None,""]
];


Options[LevelLabel]={Rational->SolidusFractionize,Parity->None};
LevelLabel[{J_,i_,P_},Opts:OptionsPattern[]]:=Module[
{JText,iValue,PValue,PText},

(* construct entry text *)
JText=J/.{x_Rational:>OptionValue[Rational][x]};
iValue=ReplaceSequential[i,{Automatic->None}];
PValue=ReplaceSequential[P,{Automatic->OptionValue[Parity]}];
PText=Switch[PValue,+1,"+",-1,"-",_,PValue];

(* generate scriptbox as appropriate *)
Which[
(iValue===None)&&(PValue===None),
JText,
(iValue===None),
Superscript[JText,PText],
(PValue===None),
Subscript[JText,iValue],
True,
Subsuperscript[JText,iValue,PText]
]
];
LevelLabel[{J_,i_},Opts:OptionsPattern[]]:=LevelLabel[{J,i,Automatic},Opts];
LevelLabel[{J_},Opts:OptionsPattern[]]:=LevelLabel[{J,Automatic,Automatic},Opts];


LevelLabel[J_,i:Except[_?OptionQ],P:Except[_?OptionQ],Opts:OptionsPattern[]]:=LevelLabel[{J,i,P},Opts];
LevelLabel[J_,i:Except[_?OptionQ],Opts:OptionsPattern[]]:=LevelLabel[{J,i},Opts];
LevelLabel[J_,Opts:OptionsPattern[]]:=LevelLabel[{J},Opts];


Options[EnergyLabel]={Symbol->textit["E"]};


EnergyLabel[Level1_,Opts:OptionsPattern[]]:=Module[
{},
Row[{OptionValue[Symbol],"(",Switch[Level1,_List,LevelLabel[Level1],_,Level1],")"}]
];


Options[RTPLabel]={MultipolaritySymbols->{"E","M"},MultipolarityStyle->Italic};


RTPLabel[Multipolarity:{sigma:(+1|-1),lambda_Integer},Level1_,Level2_,Opts:OptionsPattern[]]:=Module[
{Mode,StyleDirective},
Mode=Switch[
sigma*(-1)^lambda,
+1,OptionValue[MultipolaritySymbols][[1]],  (* "E" *)
-1,OptionValue[MultipolaritySymbols][[2]]  (* "M" *)
];
StyleDirective=OptionValue[MultipolarityStyle];
RTPLabel[Row[{Style[Mode,StyleDirective],lambda}],Level1,Level2,Opts]
];
RTPLabel[lambda_Integer,Level1_,Level2_,Opts:OptionsPattern[]]:=RTPLabel[{+1,lambda},Level1,Level2,Opts];


RTPLabel[Multipolarity:{sigma:(+1|-1),lambda_Integer},Opts:OptionsPattern[]]:=Module[
{Mode,StyleDirective},
Mode=Switch[
sigma*(-1)^lambda,
+1,OptionValue[MultipolaritySymbols][[1]],  (* "E" *)
-1,OptionValue[MultipolaritySymbols][[2]]  (* "M" *)
];
StyleDirective=OptionValue[MultipolarityStyle];
RTPLabel[Row[{Style[Mode,StyleDirective],lambda}],Opts]
];
RTPLabel[lambda_Integer,Opts:OptionsPattern[]]:=RTPLabel[{+1,lambda},Opts];


RTPLabel[Multipolarity:Except[(_Integer)|(_List)],Level1_,Level2_,OptionsPattern[]]:=Row[{
textit["B"],"(",Multipolarity,";",
Switch[Level1,_List,LevelLabel[Level1],_,Level1],
"\[Rule]",
Switch[Level2,_List,LevelLabel[Level2],_,Level2],
")"
}];
RTPLabel[Multipolarity:Except[(_Integer)|(_List)],OptionsPattern[]]:=Row[{
textit["B"],"(",Multipolarity,
")"
}];


Options[RMELabel]={MultipolaritySymbols->{"Q","M"},MultipolarityStyle->Bold,Dividers->"||"};


RMELabel[Multipolarity:{sigma:(+1|-1),lambda_Integer},Level2_,Level1_,Opts:OptionsPattern[]]:=Module[
{Mode,StyleDirective},
Mode=Switch[
sigma*(-1)^lambda,
+1,OptionValue[MultipolaritySymbols][[1]],  (* "Q" *)
-1,OptionValue[MultipolaritySymbols][[2]]  (* "M" *)
];
StyleDirective=OptionValue[MultipolarityStyle];
RMELabel[Subscript[Style[Mode,StyleDirective],lambda],Level2,Level1,Opts]
];
RMELabel[lambda_Integer,Level2_,Level1_,Opts:OptionsPattern[]]:=RMELabel[{+1,lambda},Level2,Level1,Opts];


RMELabel[Multipolarity:{sigma:(+1|-1),lambda_Integer},Opts:OptionsPattern[]]:=Module[
{Mode},
Mode=Switch[
sigma*(-1)^lambda,
+1,OptionValue[MultipolaritySymbols][[1]],  (* "E" *)
-1,OptionValue[MultipolaritySymbols][[2]]  (* "M" *)
];
RMELabel[Subscript[Mode,lambda],Opts]
];
RMELabel[lambda_Integer,Opts:OptionsPattern[]]:=RMELabel[{+1,lambda},Opts];


RMELabel[Multipolarity:Except[(_Integer)|(_List)],Level2_,Level1_,OptionsPattern[]]:=Row[{
"\[LeftAngleBracket]",Switch[Level2,_List,LevelLabel[Level2],_,Level2],OptionValue[Dividers],Multipolarity,OptionValue[Dividers],Switch[Level1,_List,LevelLabel[Level1],_,Level1],"\[RightAngleBracket]"
}]
RMELabel[Multipolarity:Except[(_Integer)|(_List)],OptionsPattern[]]:=Row[{
"\[LeftAngleBracket]",OptionValue[Dividers],Multipolarity,OptionValue[Dividers],"\[RightAngleBracket]"
}]


MomentLabel[lambda:1,Level1_,Opts:OptionsPattern[]]:=Row[{"\[Mu]","(",Switch[Level1,_List,LevelLabel[Level1],_,Level1],")"}];
MomentLabel[lambda:2,Level1_,Opts:OptionsPattern[]]:=Row[{textit["Q"],"(",Switch[Level1,_List,LevelLabel[Level1],_,Level1],")"}];


Fractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=DisplayForm[x];
Fractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=FractionBox[Numerator[x],Denominator[x],Opts];
f:Fractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];


TextFractionBox[a_,b_,Opts:OptionsPattern[]]:=Style[
GridBox[
{{a},{b}},
Opts,
RowLines->True,RowSpacings->0,RowMinHeight->{1,1.2},RowAlignments->{Bottom,Baseline},ColumnSpacings->0
],
Dashing[{}],AbsoluteThickness[0.5],Smaller
];


Style[TextFractionBox[Row[{"3","\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)"}],"2"],FontFamily->Times]//DisplayForm


TextFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=DisplayForm[x];
TextFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=
TextFractionBox[Numerator[x],Denominator[x],Opts];
f:TextFractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];


Options[SolidusFractionBox]={Spacings->{0,0}};
SolidusFractionBox[x_,y_,Opts___?OptionQ]:=Module[
{FullOpts=Flatten[{Opts,Options[SolidusFractionBox]}]},
Grid[
{{x,"/",y}},
Spacings->{Flatten[{0,(Spacings/.FullOpts),0}]},
BaselinePosition->{{1,3},Baseline} (* use baseline of denominator *)
]
];


SolidusFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=DisplayForm[x];
SolidusFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=SolidusFractionBox[Numerator[x],Denominator[x],Opts];
f:SolidusFractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];


Options[DiagonalFractionBox]={Spacings->{-0.1,-0.3},Baseline->{0.5,0.0,0.0},KernForSuperscript->-0.15};
DiagonalFractionBox[x_,y_,Opts___?OptionQ]:=Module[
{FullOpts=Flatten[{Opts,Options[DiagonalFractionBox]}]},
TagBox[
StyleBox[
Grid[
{{
AdjustmentBox[SubscriptBox["",DisplayForm[x]],BoxBaselineShift->-((Baseline/.FullOpts)[[1]])],
(*AdjustmentBox[SubscriptBox["",StyleForm["/",FontSlant->"Oblique",Larger]],BoxBaselineShift\[Rule]-((Baseline/.FullOpts)[[2]])],*)
AdjustmentBox[StyleForm["/",FontSlant->"Oblique"],BoxBaselineShift->-((Baseline/.FullOpts)[[2]])],
AdjustmentBox[SubscriptBox["",DisplayForm[y]],BoxBaselineShift->-((Baseline/.FullOpts)[[3]])]
}},
Spacings->{Flatten[{0,(Spacings/.FullOpts),0}]},
BaselinePosition->{{1,3},Bottom}
],
ScriptBaselineShifts->{0,0}
],
DiagonalFractionBox[(KernForSuperscript/.FullOpts)]
]
];


Unprotect[TagBox];
TagBox/:SuperscriptBox[x:TagBox[_,DiagonalFractionBox[Adjustment_]],n_]:=SuperscriptBox[AdjustmentBox[x,BoxMargins->{{0,Adjustment},{0,0}}],n];
TagBox/:Superscript[x:TagBox[_,DiagonalFractionBox[Adjustment_]],n_]:=Superscript[AdjustmentBox[x,BoxMargins->{{0,Adjustment},{0,0}}],n];
Protect[TagBox];


DiagonalFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=DisplayForm[x];
DiagonalFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=DiagonalFractionBox[Numerator[x],Denominator[x],Opts];
f:DiagonalFractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];


FractionString[x_?NumericQ]:=Module[
{f,NumeratorString,DenominatorString},
f=Rationalize[x];
NumeratorString=ToString[Numerator[f]];
DenominatorString=ToString[Denominator[f]];
Which[
f==0,"0",
Denominator[f]==1,StringJoin[NumeratorString],
Denominator[f]!=1,StringJoin[NumeratorString,"/",DenominatorString]
]
];


PiFractionString[x_?NumericQ]:=Module[
{f,NumeratorString,DenominatorString},
f=Rationalize[x/Pi];
NumeratorString=If[
Numerator[f]==1,
"",
ToString[Numerator[f]]
];
DenominatorString=ToString[Denominator[f]];
Which[
f==0,"0",
Denominator[f]==1,StringJoin[NumeratorString,"\[Pi]"],
Denominator[f]!=1,StringJoin[NumeratorString,"\[Pi]","/",DenominatorString]
]
];


End[];


Protect[Evaluate[$Context<>"*"]];
Unprotect[Evaluate[$Context<>"$*"]];
EndPackage[];
