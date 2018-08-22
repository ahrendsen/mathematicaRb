(* ::Package:: *)

(* ::Title:: *)
(*Excitation Function*)


(* ::Chapter:: *)
(*Initialization*)


dropBoxOn=FileNameJoin[{"Research","DataRuns"}];
If[SystemInformation["Kernel","OperatingSystem"]=="Unix",
researchFolder=FileNameJoin[{"/","home","karl","Research"}]; ,(* Home coputer is Linux, so needs "/" at beginning *)
researchFolder=FileNameJoin[{"C:","Users","kahrendsen2","Box Sync","Research"}];(* Lab computer is Windows, so needs "C:" at beginning"*)
];
SetDirectory[researchFolder];

SetDirectory[FileNameJoin[{"mathematicaRb","packages"}]];
Import["fileManipulation.wl"];
Import["physicalConstants.wl"];
ResetDirectory[];
ResetDirectory[];


(* ::Chapter:: *)
(*Equations To Work With Data*)





(* ::Chapter:: *)
(**)
