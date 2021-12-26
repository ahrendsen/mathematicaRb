(* ::Package:: *)

(* Apparatus Constants *)
alphaOld=3.2*\[Pi]/180; (* I can't even remember when I was using these. *)
betaOld=-6.95*\[Pi]/180;
rev=1;
alpha2=20.4*\[Pi]/180; (* These constants are the ones I used up until the chamber redesign in 2020. *)
beta2=68.7*\[Pi]/180;
shift=-4.4552*\[Pi]/180
alpha=20.4*\[Pi]/180+shift; (* These constants are the ones I used after the chamber redesign in 2020. *)
beta=68.7*\[Pi]/180+shift;
deltaAB=alpha-beta;
delta=94.54*\[Pi]/180 (*Munir's reported 1.66 plus or minus .01*);

(* These are the constants I calculated after a more extensive polarimeter analysis (Rb9 pg. 48-69), I'm choosing to use the 
previously determined ones because they give me a P2 closer to zero. *)
alpha=-25.292*\[Pi]/180;
deltaAB=-13.06*\[Pi]/180;
beta=alpha-deltaAB;
delta=111.3*\[Pi]/180;

(*Physical constants *)
\[CapitalGamma]=5750056; (* Hz *)
fge=0.34231; (* dimensionless *)
k=4/3; (* dimensionless *)
re=2.8179*^-13; (* cm *)

\[Nu]0=377107.463;  (* GHz *)
\[Lambda]0=c/100/\[Nu]0  (* nm *)
c=2.99792458*^10; (* cm/s *)
h=6.6261*^-27; (* cm^2*g/s *)
\[Mu]=9.2740*^-21; (* g*cm^2*s^-2*G^-1 *)
gHzToHz=1*^9; 
nA=1*^-9;
mTorr=1*^-3;
kB=1.38064852*^-23

transGHz={-3.072,-1.372,1.84,4.576};

BdotL; (* G*cm *)

Q; (* Constant factors in Faraday Rotation which won't change from model to model, basically anything but numerical constants and \[Pi]s *)
nDensC;(*The single constant that is multiplied by the fitting parameter and divided by the integrated magnetic field to give the number density.*)
nDensCLitaker;(*Litaker's version of the single constant that is multiplied by the fitting parameter and divided by the integrated magnetic field to give the number density.*)









