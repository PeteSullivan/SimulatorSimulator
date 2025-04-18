(* ::Package:: *)

(* Extracting command-line arguments by skipping unnecessary items *)

args = $ScriptCommandLine;  

NIons = ToExpression[args[[2]]];  (* First argument *)
Vrf = ToExpression[args[[3]]];  (* Second argument *)
V1 = ToExpression[args[[4]]];  (* ... *)
V2 = ToExpression[args[[5]]];  
V3 = ToExpression[args[[6]]];  
ImportComsolNow =ToExpression[args[[7]]];  
V := {V1,V2,V3};
(* Print to verify *)
(*
Print["args: ", args]
Print["NIons: ", NIons];
Print["Vrf: ", Vrf];
Print["V1: ", V1];
Print["V2: ", V2];
Print["V3: ", V3];
Print["ImportComsolNow: ", ImportComsolNow]
*)


SetDirectory["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\MathematicaScripts\\TrapData"]


(* Import Comsol files with trap data, only needed first time running a trap *)
ImportComsol[filename_,InterpOrder_]:=Module[{A,Elec,ElecInterp,Interp},
A=AbsoluteTime[];
Elec=Import[filename,"Table"];
If[Elec==Null,Break[]];
Print["Import of " <>filename<>": "<>ToString[(AbsoluteTime[]-A)]<>" seconds"];
ElecInterp=Function[{x,y,z,f},{{x,y,z},f}]@@@Elec;
Interp=Interpolation[ElecInterp,InterpolationOrder->InterpOrder];
Print["Interpolation of " <>filename<>": "<>ToString[(AbsoluteTime[]-A)]<>" seconds"];
Export[StringDrop[filename,-4]<>".m",Interp];
Print["Export of " <>filename<>": "<>ToString[(AbsoluteTime[]-A)]<>" seconds"]];

If[ImportComsolNow==True,
ImportComsol["DC_XY.txt",{3,3,0}]
ImportComsol["DC_YZ.txt",{0,3,3}]
ImportComsol["DC_XZ.txt",{3,0,3}]

ImportComsol["RF_XY.txt",{3,3,0}];
ImportComsol["RF_YZ.txt",{0,3,3}];
ImportComsol["RF_XZ.txt",{3,0,3}];
]



DC1xy=Get["DC_XY.m"];
DC1zy=Get["DC_YZ.m"];
DC1zx=Get["DC_XZ.m"];

(* Fields were exported in V/m*)
RFxy=Get["RF_XY.m"];
RFzy=Get["RF_YZ.m"];
RFzx=Get["RF_XZ.m"];



(* Units constants *)
d = 1000;  (* mm/m *)
q = 1.602 10^(-19);
q
(* Atomic Constants *)
m = 171 1.67 10^(-27);
m
(* Experimental Constants *) 
\[CapitalOmega] = 2\[Pi] 15.0 10^6; (* RF drive frequency *)
Vrf:=1000; 
DCMax:=2500;


xCtr=0.;
xRange=.6;
yCtr=0.;
yRange=0.06;
zCtr=0.;
zRange=0.06;
xMin=xCtr-xRange/2;
xMax=xCtr+xRange/2;
yMin=yCtr-yRange/2;
yMax=yCtr+yRange/2;
zMin=zCtr-zRange/2;
zMax=zCtr+zRange/2;


appliedVoltages


DC1[t_]:=V[[1]];


(* Note: in v11, voltages applied to all endcap electrodes (DC1) and all center electrodes (RF1) in COMSOL *)
DCxy[x_,y_,z_,t_] :=DC1[t]*DC1xy[x,y,0]
DCzx[x_,y_,z_,t_] :=DC1[t]*DC1zx[x,0,z]
DCzy[x_,y_,z_,t_] :=DC1[t]*DC1zy[0,y,z]


RFVxy[x_,y_,z_,t_] :=(q (Vrf^2 RFxy[x,y,0]^2))/(4 m (\[CapitalOmega]^2) );
RFVzx[x_,y_,z_,t_] :=(q ( Vrf^2 RFzx[x,0,z]^2))/(4 m (\[CapitalOmega]^2) );
RFVzy[x_,y_,z_,t_] :=(q ( Vrf^2 RFzy[0,y,z]^2))/(4 m (\[CapitalOmega]^2) );


Uxy[x_,y_,z_,t_] :=DCxy[x,y,0,t]+RFVxy[x,y,0,t];
Uzx[x_,y_,z_,t_] :=DCzx[x,0,z,t]+RFVzx[x,0,z,t];
Uzy[x_,y_,z_,t_] :=DCzy[0,y,z,t]+RFVzy[0,y,z,t];


(*Uxy[x,y,z,t]
Uzx[x,yCtr,zCtr,0]*)


Quad[x_,y_,a1_,b1_,a2_,b2_,c_]:=a2 x^2+b2 y^2+a1 x+b1 y+c x y;
QuadRot[x_,y_,\[Alpha]_,\[Beta]_,\[Theta]_]:=\[Alpha] (x Cos[\[Theta]]+y Sin[\[Theta]])^2+\[Beta] (-x Sin[\[Theta]]+y Cos[\[Theta]])^2;


NumPts=100;
yzFits=LinearModelFit[Flatten[Table[{y,z,Uzy[0,y,z,0]},{y,yMin,yMax,(yMax-yMin)/NumPts},{z,zMin,zMax,(zMax-zMin)/NumPts}],1],{y^2,y z,z^2,y,z},{y,z}];

xyFit=LinearModelFit[Flatten[Table[{x,y,Uxy[x,y,0,0]},{x,-0.095,0.095,(0.095+0.095)/NumPts},{y,yMin,yMax,(yMax-yMin)/NumPts}],1],{x^2,x y,y^2,x,y},{x,y}];

yzFitParams=yzFits["BestFitParameters"];
xyFitParams=xyFit["BestFitParameters"];


Minimized=Minimize[yzFits[y,z],yMin<= y<= yMax&&zMin<= z<= zMax,{y,z}];
yzMins={y/.Minimized[[2,1]],z/.Minimized[[2,2]]};


\[Theta]x=ArcTan[(xyFitParams[[2]]-xyFitParams[[4]]),xyFitParams[[3]]]/2;
\[Theta]y=ArcTan[(yzFitParams[[2]]-yzFitParams[[4]]),yzFitParams[[3]]]/2;
\[Theta]z=\[Theta]y+\[Pi]/2;


PrincipleZaxis=Table[{yzMins[[1]]+ampl*Cos[\[Theta]z],yzMins[[2]]+ampl*Sin[\[Theta]z]},{ampl,zMin,zMax,0.001}];
PrincipleYaxis=Table[{yzMins[[1]]+ampl*Cos[\[Theta]y],yzMins[[2]]+ampl*Sin[\[Theta]y]},{ampl,yMin,yMax,0.001}];


\[Theta]x=\[Pi]/2;
FreqSecYaxis=1/(2\[Pi] )Sqrt[(d^2 q)2(yzFitParams[[2]]Cos[\[Theta]y]^2+yzFitParams[[3]]Sin[\[Theta]y]*Cos[\[Theta]y]+yzFitParams[[4]] Sin[\[Theta]y]^2)/m];
FreqSecZaxis=1/(2\[Pi] )Sqrt[(d^2 q)2(yzFitParams[[2]]Sin[\[Theta]y]^2-yzFitParams[[3]]Sin[\[Theta]y]*Cos[\[Theta]y]+yzFitParams[[4]] Cos[\[Theta]y]^2)/m];
FreqSecXaxis=1/(2\[Pi] )Sqrt[(d^2 q)2(xyFitParams[[2]]Sin[\[Theta]x]^2-xyFitParams[[3]]Sin[\[Theta]x]*Cos[\[Theta]x]+xyFitParams[[4]] Cos[\[Theta]x]^2)/m];
\[Kappa]=(2\[Pi]*FreqSecXaxis)^2*m*(200*10^-6)^2/(2*q*V[[1]]);

(*
Print["Radial secular frequency along principle axis Y (MHz): ", FreqSecYaxis*10^-6];
Print["Radial secular frequency along principle axis Z (MHz): ", FreqSecZaxis*10^-6];
Print["Radial secular frequency along X axis (MHz): ", FreqSecXaxis*10^-6];
Print["\[Kappa] = ", \[Kappa]];
*)


(*Export secular frequencies to a csv file for Python access *)
SecularFrequencies = {FreqSecXaxis*10^-6, FreqSecYaxis*10^-6, FreqSecZaxis*10^-6};
(*Print["Secular frequencies: ", SecularFrequencies];*)
Export["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\secular_frequencies_mid.csv", SecularFrequencies];
