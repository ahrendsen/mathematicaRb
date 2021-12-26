(* ::Package:: *)

(* ::Title:: *)
(*System Functions and Constants*)


(* ::Chapter:: *)
(*Dimensions*)


collisionCellLength=3;(*cm*)
(*2020-01-28: Switched to 3 cm as I was previously neglecting the extra
length added by the entrance and exit electrodes and the kapton sheets. 
Before 2020-01-28: 2.794*)
pyrexCellLength=6.85; (*cm, Litaker's thesis, pg. 35*)
cellLength=collisionCellLength;


(* ::Chapter:: *)
(*Magnets*)


(* ::Subchapter:: *)
(*Regular System Magnets*)


(* The coefficients for the polynomial describing the magnetic field of magnet 1 [(c)oefficient of magent (1)]  based on the current.
Assumes x is measured in cm *)
c1x2=.085;
c1x1=-1.907;
c1x0=17.24;
c2x2=0.153;
c2x1=1.866;
c2x0=15.594;

(* The coefficients for the linear relationship between the current and the voltage of the magnets *)
c1v2iSlope=.0504;
c1v2iInt=-.0113;
c2v2iSlope=.0488;
c2v2iInt=-.0069;

(* Calculate the current in each of the magnets based on the voltage *)
i1[voltage_]:=voltage*c1v2iSlope+c1v2iInt;
i2[voltage_]:=voltage*c2v2iSlope+c2v2iInt;


Clear[CalculateIntegratedMagnetic];
Clear[CalculateIntegratedMagneticField];
CalculateIntegratedMagneticField[voltageList_]:=Module[{x,fields,i},
(* Integrates the value of the magnetic field over the length of the cell. *)
fields={};
For[i=1,i<=Length[voltageList],i++,

	AppendTo[fields,Integrate[getBEq[i1[voltageList[[i]][[1]]],i2[voltageList[[i]][[2]]],cellPos],{cellPos,0,cellLength}]];
];
fields
(*Integrate[getBEqAdHoc[voltageMagnet1,voltageMagnet2,cellPos],{cellPos,0,cellLength}]*)
];

CalculateIntegratedMagneticField[voltageMagnet1_,voltageMagnet2_]:=Module[{x},
(* Integrates the value of the magnetic field over the length of the cell. *)
Integrate[getBEq[i1[voltageMagnet1],i2[voltageMagnet2],cellPos],{cellPos,0,cellLength}]
(*Integrate[getBEqAdHoc[voltageMagnet1,voltageMagnet2,cellPos],{cellPos,0,cellLength}]*)
];


(* Provides the Magnetic field strength in Gauss.

Given: i1, the current in magnet 1 in amps.
       i2, the current in magnet 2 in amps.
       cellPos, the desired location to know the magnetic field in the cell, in cm.
       *)
getBEq[i1_,i2_,cellPos_]:=((c1x2 )i1 +(c2x2)i2) cellPos^2 +((c1x1 )i1 +(c2x1)i2) cellPos +((c1x0 )i1 +(c2x0)i2);  (* Use this eq. when using the normal spin filter. *)
(*getBEq[i1_,i2_,cellPos_]:=170;*)




(* ::Subchapter:: *)
(*Ad-Hoc Polarization Setup*)


getBEqAdHoc[i1_,i2_,cellPos_]:=(i1*(5.91220679503176` -0.46227239116010416` cellPos+0.011298812156657146` cellPos^2)+i2*(4.4953230074688975` +0.13763345084459075` cellPos+0.06404804608042143` cellPos^2+0.00011932029103778548` cellPos^3-0.0005576642279085437` cellPos^4+0.000018991369967445234` cellPos^5));


(* ::Chapter:: *)
(*Lasers*)


(* ::Section:: *)
(*Toptica Laser*)


pumpDriveCurrentsXX={2300,3000,3700};
pumpPowersXX={994,1560,2090};
pumpDriveToPowerLMF=LinearModelFit[Transpose[{pumpDriveCurrentsXX,pumpPowersXX}],x,x];
PumpPowerFromDriveCurrent[x_]:=pumpDriveToPowerLMF["Function"][x];


(* ::Section:: *)
(*Vortex Laser*)


ProbePower[driveCurrent_,cleanupLP_]:=(-5.855+.2288*driveCurrent)*.52*(1+Cos[2*(1.358+cleanupLP)]);


(* ::Chapter:: *)
(*Polarimeter*)


polConstxxp1Const=1.0614;
polConstxxp1Coeff=.9386;
polConstxxp3Coeff=2.6409;
polConstxxunPolP1=.04;


getBEq[i1[65],i2[65],0]
