(* ::Package:: *)

(* Apparatus Constants *)
alphaOld=3.2*\[Pi]/180;
betaOld=-6.95*\[Pi]/180;
rev=1;
alpha=20.4*\[Pi]/180;
beta=68.7*\[Pi]/180;
deltaAB=alpha-beta;
delta=94.54*\[Pi]/180 (*Munir's reported 1.66 plus or minus .01*);

(*Physical constants *)
\[CapitalGamma]=5750056; (* Hz *)
fge=0.34231; (* dimensionless *)
k=4/3; (* dimensionless *)
re=2.8179*^-13; (* cm *)
cellLength=2.794;(*cm*)
\[Nu]0=377107.463;  (* GHz *)
\[Lambda]0=c/100/\[Nu]0  (* nm *)
c=2.99792458*^10; (* cm/s *)
h=6.6261*^-27; (* cm^2*g/s *)
\[Mu]=9.2740*^-21; (* g*cm^2*s^-2*G^-1 *)
gHzToHz=1*^9; 
nA=1*^-9;
mTorr=1*^-3;

BdotL; (* G*cm *)

Q; (* Constant factors in Faraday Rotation which won't change from model to model, basically anything but numerical constants and \[Pi]s *)
nDensC;(*The single constant that is multiplied by the fitting parameter and divided by the integrated magnetic field to give the number density.*)
nDensCLitaker;(*Litaker's version of the single constant that is multiplied by the fitting parameter and divided by the integrated magnetic field to give the number density.*)
