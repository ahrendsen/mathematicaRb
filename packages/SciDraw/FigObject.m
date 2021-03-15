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



(* :Title: FigObject *)
(* :Context: SciDraw`FigObject` *)
(* :Summary: Figure object definition infrastructure *)
(* :Author: Mark A. Caprio, Department of Physics, University of Notre Dame *)
(* :Copyright: Copyright FIGYEAR, Mark A. Caprio *)
(* :Package Version: FIGVERSION *)
(* :Mathematica Version: MATHVERSION *)
(* :Discussion: FIGDISCUSSION *)
(* :History: See main package file. *)


BeginPackage["SciDraw`",SciDraw`Private`$ExternalContexts];


Unprotect[Evaluate[$Context<>"*"]];


Begin["`Private`"];


$FigClassRegistry={};


RegisterFigOptions[s_Symbol]:=($FigClassRegistry=Union[$FigClassRegistry,{s}]);


RegisterFigOptions[LinTicks];
RegisterFigOptions[LogTicks];


RegisterFigOptions[FigObject];


$FigClassAttachedLabels[FigObject]={};


General::nocreate="Invalid or unrecognized arguments given to figure object.  Please check the arguments and try again.";


DeclareFigClass[
Class_Symbol,
ParentClass:_Symbol:FigObject,
DataMemberNames:{___String},
MethodNames:{___String},
AttachedLabelNames:{((_Symbol)|(_String))...}
]:=Module[
{},

(* declare class as descendent of FigObject (or possibly another parent) *)
DeclareClass[Class,ParentClass,DataMemberNames,MethodNames,Replace->True];

(* register class, so options can be scoped *)
$FigClassRegistry=Union[$FigClassRegistry,{Class}];

(* set hook so that error in creation (duplicate name or syntax error) triggers Abort *)
(* The error message provided by MathObject suffices. *)
OnCreationFailure[Class,Self_Object][___]:=Abort[];

(* record attached label names *)
$FigClassAttachedLabels[Class]=AttachedLabelNames;

];


DefineFigClassOptions[
Class_Symbol,
DefaultOverrides:_List:{},
NewOptionRules_List
]:=Module[
{ParentClass,NewAttachedLabels},

ParentClass=ClassAncestry[Class][[-2]];
NewAttachedLabels=Complement[$FigClassAttachedLabels[Class],$FigClassAttachedLabels[ParentClass]];
DefineOptions[
Class,
{ParentClass,DefaultOverrides,All},
{

(* attached labels *)
(* use Complement to allow override via NewOptionRules *)
Complement[
FigDerivedLabelOptions[NewAttachedLabels,True,Default],
NewOptionRules,
SameTest->(SameQ[First[#1],First[#2]]&)
],

(* all other special options *)
NewOptionRules

}
]
];


General::figunrecopts="One or more unrecognized options `1` were encountered.";


FigRealizeOptions[Self_Object,Class_Symbol,OptionList_List]:=Module[
{
UnrecognizedOptions,RealizedOptions,OverrideOptions,BaseStyleOptions,StyleList
},

(* trap illegal options *)
UnrecognizedOptions=Complement[First/@Flatten[OptionList],First/@Options[Class]];
If[
UnrecognizedOptions=!={},
FigError[Self,"figunrecopts",UnrecognizedOptions]
];


(* extract option overrides *)
(* extract option *override* rules whose LHS pattern matches ObjectName[Self], and assemble a list of the resulting *option* rules *)
(* but do not use these if we are resolving options for another class *)
OverrideOptions=If[
ObjectClass[Self]===Class,
ResolveOptionOverrides[ObjectName[Self]],
{}
];

(* pass 1: realize options sans Style override, to find value of Style option *)
RealizedOptions=RealizeOptions[Class,{OptionList,OverrideOptions}];

(* validate style name *)
(*Print[RealizedOptions];*)
FigCheckOption[Self,Style,StyleSpecifierPattern,RealizedOptions];

(* pass 2: realize options after setting default options from style *)
(* Note: This may affect the realized options by resetting the options for Class directly or by resetting the options for an inheritance parent. *)
RealizedOptions=WithStyle[
(Style/.RealizedOptions),
RealizeOptions[Class,{OptionList,OverrideOptions}]
];

RealizedOptions

];


BaseOutlineOptionList={
{ShowLine,LogicalPattern|Default},
{LineColor,None|ColorDirectivePattern|Default}, 
{LineOpacity,None|UnitIntervalPattern|Default},
{LineThickness,FigThicknessPattern},
{LineDashing,FigDashingPattern},
{LineCapForm,None|"Butt"|"Round"|"Square"}, 
{LineJoinForm,None|"Miter"|"Bevel"|"Round"|{"Miter",NonNegativePattern}},
{LineDirectives,FlatListPattern}
};
BaseFillOptionList={
{ShowFill,LogicalPattern|Default},
{FillColor,None|ColorDirectivePattern|Default},
{FillOpacity,None|UnitIntervalPattern|Default},
{FillDirectives,FlatListPattern}
};
BasePointOptionList={
{ShowPoint,LogicalPattern|Default},
{PointColor,None|ColorDirectivePattern|Default}, 
{PointOpacity,None|UnitIntervalPattern|Default},
{PointSize,FigPointSizePattern},
{PointDirectives,FlatListPattern}
};
BaseTextOptionList={

(* text appearance *)
{ShowText,LogicalPattern|Default},
{TextColor,None|ColorDirectivePattern|Default},
{TextOpacity,None|UnitIntervalPattern|Default},
{FontFamily,FontFamilyPattern},
{FontSize,FontSizePattern},
{FontWeight,FontWeightPattern},
{FontSlant,FontSlantPattern},
{FontTracking,FontTrackingPattern},
{FontVariations,FlatListPattern},
{TextStyleOptions,FlatListPattern},

(* text positioning *)
{TextBaseBuffer,Automatic|ScalarParameterPattern},
{TextBuffer,ScalarParameterPattern},
{TextNudge,IntervalParametersPattern},
{TextOffset,FigTextOffsetPattern},
{TextOrientation,FigTextOrientationPattern},

(* text background, frame, and geometry *)
{TextBackground,Automatic|None|ColorDirectivePattern},
{TextFrame,LogicalPattern},
{TextFrameColor,None|ColorDirectivePattern|Default},
{TextRoundingRadius,IntervalParametersPattern},
{TextMargin,NonNegativeRangeParametersPattern},
{TextPadding,LogicalPattern},
{TextRectify,LogicalPattern},

(* text callout -- duplicate line options *)
{TextCallout,LogicalPattern},
{TextCalloutColor,None|ColorDirectivePattern|Default}, 
{TextCalloutOpacity,None|UnitIntervalPattern|Default},
{TextCalloutThickness,FigThicknessPattern},
{TextCalloutDashing,FigDashingPattern},
{TextCalloutCapForm,None|"Butt"|"Round"|"Square"}, 
{TextCalloutJoinForm,None|"Miter"|"Bevel"|"Round"|{"Miter",NonNegativePattern}},
{TextCalloutDirectives,FlatListPattern}

};


FigCheckBaseOptions[Self_Object,FullOptions_List]:=Module[
{},

(* overall appearance *)
FigCheckOption[Self,Show,LogicalPattern,FullOptions];
FigCheckOption[Self,Color,None|ColorDirectivePattern,FullOptions]; 
FigCheckOption[Self,Opacity,None|UnitIntervalPattern,FullOptions];
FigCheckOption[Self,Directives,_List,FullOptions];
FigCheckOption[Self,Layer,Automatic|(_?NumericQ),FullOptions];

(* outline/fill/point/text *)
MapThread[
FigCheckOption[Self,#1,#2,FullOptions]&,
Transpose@Join[BaseOutlineOptionList,BaseFillOptionList,BasePointOptionList,BaseTextOptionList]
];

(* style *)
(* option already resolved *)

(* Prolog/Epilog -- no validation required *)

(* diagnostic *)
FigCheckOption[Self,Debug,LogicalPattern,FullOptions];
(*FigCheckOption[Self,PrintTiming,LogicalPattern,FullOptions];*)

];


FigProcessStandardAttachedLabels[Self_Object]:=Module[
{Side},

Do[
FigCheckDerivedLabelOptions[Self,Side,True,FigOptions];
FigSpawnAttachedLabel[Self,Side,FigResolveDerivedLabelOptions[Side,True,FigOptions],FigOptions],
{Side,$FigClassAttachedLabels[ObjectClass[Self]]}
]
];


SetAttributes[FigObjectWrapper,HoldAll];
FigObjectWrapper[Class_Symbol,Self_Object,OptionList_List,Body_]:=PrintTiming[Block[
{
FigOptions=FigRealizeOptions[Self,Class,OptionList]
},

(* prolog processing *)
PrintTiming[
{
(* check in figure *)
FigCheckInFigure[Self];

(* validate object name *)
FigCheckObjectName[Self];

(* validate core options *)
FigCheckBaseOptions[Self,FigOptions];

(* save options *)
(*Self@SetOptionValues[FigOptions];*)
},
Label->Row[{"...Object setup ", Class," ",Self}],
Print->$PrintTiming
];

(* evaluate body *)
(* Note: Prolog and Epilog as bare (Prolog/.FigOptions) would require reimplementation so that they are not multiply evaluated during option inheritance processing. Current model is to presume option is wrapped in Hold.  *)
PrintTiming[
WithStyle[
(Style/.FigOptions),
{
If[(Debug/.FigOptions),Print["Entering ",ObjectClass[Self],"[",ObjectName[Self],"]..."]];
ReleaseHold[(Prolog/.FigOptions)];
Body;
ReleaseHold[(Epilog/.FigOptions)];
If[(Debug/.FigOptions),Print["Exiting ",ObjectClass[Self],"[",ObjectName[Self],"]."]]
}
],
Label->Row[{"...Object body ", Class," ",Self}],
Print->$PrintTiming
];

(*If[(Debug/.FigOptions),ShowObjectInformation[Self]];*)


(* validate and spawn any standard attached label *)
PrintTiming[
{
FigProcessStandardAttachedLabels[Self]
},
Label->Row[{"...Object labels ", Class," ",Self}],
Print->$PrintTiming
];

(* return reference to created object -- for possible future use in stream operator overloads *)
Self
],
Label->Row[{"***Object wrapper ", Class," ",Self}],
Print->$PrintTiming
];


General::noanchor="An anchor named `2` with argument `3` cannot be generated for this object.  Please check the anchor names and arguments allowed for objects of class `1`.\n\n(Debugging information: `4`)";
General::noanchorn="An anchor named `2` (with no argument) cannot be generated for this object.  Please check the anchor names and arguments allowed for objects of class `1`.\n\n(Debugging information: `4`)";


FigMakeAnchorWrapper[Class_Symbol,Self_Object,Name_,Arg_,Expr_]:=If[
MatchQ[Expr,_Object],
Expr,
If[
Arg===None,
FigError[MakeAnchor,Self,"noanchorn",Class,Name,Arg,Expr],
FigError[MakeAnchor,Self,"noanchor",Class,Name,Arg,Expr]
]
];


General::objectanchor="An anchor named `2`, with position `3`, orientation `4`, and offset `5`, could not be constructed for this object.  Please check the anchor names and parameters allowed for objects of class `1`.\n\n(Debugging information: `6`)";


SetAttributes[FigObjectAnchorWrapper,HoldAll];
FigObjectAnchorWrapper[
Class_Symbol,Self_Object,Name_,OptionsList_List,
Body_]:=Block[
{
FigAnchorOptions=Join[OptionsList,{Position->Automatic,Orientation->Automatic,Offset->Automatic,Buffer->None,Nudge->None}],
GeneratedAnchor,GeneratedPoint,GeneratedAngle,GeneratedOffset,BufferDirection,UsedPoint
},

Print[FigAnchorOptions];
(* validate buffer and nudge *)
FigCheckOption[MakeAnchor,Self,Buffer,ScalarParameterPattern,FigAnchorOptions];
FigCheckOption[MakeAnchor,Self,Nudge,IntervalParametersPattern,FigAnchorOptions];

(* evaluate body *)
GeneratedAnchor=Body;

(* validate body return *)
If[
!MatchQ[GeneratedAnchor,FigPointPattern],
PrintValue@FigError[MakeAnchor,Self,"objectanchor",Class,Name,(Position/.FigAnchorOptions),(Orientation/.FigAnchorOptions),(Offset/.FigAnchorOptions),Body]
];

(* extract components of anchor *)
{GeneratedPoint,GeneratedAngle,GeneratedOffset}={GeneratedAnchor@GetPoint[],GeneratedAnchor@GetAngle[],GeneratedAnchor@GetOffset[]};

(* apply buffer and nudge *)
(* determine buffered and nudged position *)
(* note that classic LevelScheme Nudge\[Rule]dy form is supported here, with UpgradePairVertical, but for derived labels it might be disallowed in the option validation of the calling object (verify?),
as idiosyncratic and incompatible with panel edge upgrade *)
BufferDirection={0,0};  (* TODO implement in generated anchor *)
UsedPoint=GeneratedPoint+UpgradeScalar[(Buffer/.FigAnchorOptions)]*BufferDirection+UpgradePairVertical[(Nudge/.FigAnchorOptions)];

(* return anchor *)
FigAnchor[Canvas[UsedPoint],GeneratedOffset,GeneratedAngle]
];


General::nobb="No bounding region has been implemented for this object.";


FigMakeBoundingBoxWrapper[Class_Symbol,Self_Object,Expr_]:=If[
MatchQ[Expr,NumericalRegionPattern],
Expr,
FigError[MakeBoundingBox,Class,"nobb"]
];


SidifyOptionName[Side:(_Symbol)|(_String)][Option_Symbol]:=Module[
{SideName},

SideName=Switch[
Side,
_Symbol,ToString[Side],
_String,Side
];

ToExpression["SciDraw`"<>SideName<>ToString[Option]]
];


FigDerivedLabelOptions[Side:(_Symbol)|(_String),IncludeLabelOptions:LogicalPattern,DefaultFontSize_]:=Flatten[{

(* label content *)
If[
IncludeLabelOptions,
{SidifyOptionName[Side]@Label->None,SidifyOptionName[Side]@LabelPosition->Automatic},
{}
],

(* override font size default *)
SidifyOptionName[Side]@FontSize->DefaultFontSize,

(* standard text options *)
Map[
(SidifyOptionName[Side][#]->Default)&,
Complement[First/@BaseTextOptionList,{FontSize}]
]
}];


FigDerivedLabelOptions[SideList:{((_Symbol)|(_String))...},IncludeLabelOptions:LogicalPattern,DefaultFontSize_]:=Flatten[(FigDerivedLabelOptions[#,IncludeLabelOptions,DefaultFontSize]&)/@SideList];


FigCheckDerivedLabelOptions[Self_Object,Side:(_Symbol)|(_String),IncludeLabelOptions:LogicalPattern,FullOptions_List]:=Module[
{},

(* label content *)
 (* anything goes! *)

(* label position *)
(* anything goes, but MakeAnchor method is free to balk later *)


(* standard text options *)
(* all initially acceptable patterns are accepted, plus Default *)
MapThread[
FigCheckOption[Self,SidifyOptionName[Side]@#1,#2|Default,FullOptions]&,
Transpose@BaseTextOptionList
];

];


FigResolveDerivedLabelOptions[Side:(_Symbol)|(_String),IncludeLabelOptions:LogicalPattern,FullOptions_List]:=Module[
{Content,UsedArgument,OptionsList},

Content=If[
IncludeLabelOptions,
ResolveOption[SidifyOptionName[Side][Label],{},FullOptions],
None
];

UsedArgument=If[
IncludeLabelOptions,
ResolveOption[SidifyOptionName[Side][LabelPosition],{Automatic->None},FullOptions],
None
];

OptionsList=(#->ResolveOption[SidifyOptionName[Side][#],{Default:>(#/.FullOptions)},FullOptions])&/@(First/@BaseTextOptionList);

{Content,UsedArgument,OptionsList}
];


General::figlabelpos="Cannot resolve `1` attached label position.  (Debugging information: Perhaps the option specification `2`->`3` is not of a supported form for this class of object.  The anchor `4` evaluated to `5`.)";


FigSpawnAttachedLabel[Self_Object,Side:(_Symbol)|(_String),{Content_,LabelPositionArgument_,OptionsList_},FullOptions_List]:=Module[
{Anchor},

(* extract label parameters *)
(*{Content,LabelPositionArgumentList,OptionsList}=FigResolveDerivedLabelOptions[Side,{Label,LabelPosition},FullOptions];*)

(* short circuit *)
If[
(Content===None)||((Show/.OptionsList)===False)||((ShowText/.OptionsList)===False),
Return[]
];

(* generate point *)
Anchor=Self@MakeAnchor[Side,LabelPositionArgument];

(* make text *)
FigTextElement[
Anchor,
Content,
Flatten[{OptionsList,FullOptions}]
];

];


SetAttributes[ScopeOptions,HoldAll];
ScopeOptions[Body_]:=BlockOptions[
$FigClassRegistry,
Body
];
DeclareFigFallThroughError[ScopeOptions];


SetOptionOverrides[rl:{((Rule|RuleDelayed)[_,_List])...}]:=($FigOptionRules=Join[$FigOptionRules,rl]);
SetOptionOverrides[r:((Rule|RuleDelayed)[_,_List])]:=SetOptionOverrides[{r}];
DeclareFigFallThroughError[SetOptionOverrides];


SetAttributes[ScopeOptionOverrides,HoldAll];
ScopeOptionOverrides[Body_]:=Block[
{$FigOptionRules=$FigOptionRules},
Body
];
DeclareFigFallThroughError[ScopeOptionOverrides];


ResolveOptionOverrides[Name_]:=Module[
{Overrides},
(*Print["resolving ",Name," with rules ",$FigOptionRules];*)
Overrides=OptionsUnion[Cases[$FigOptionRules,((patt_->opts_List)/;MatchQ[Name,patt]):>opts]];
If[
$DebugOptionOverrides,
If[Length[Overrides]>0,
Print["Option overrides for ",Name,": ",Overrides]
]
];
Overrides
]


End[];


Protect[Evaluate[$Context<>"*"]];
Unprotect[Evaluate[$Context<>"$*"]];
EndPackage[];
