(* ::Package:: *)

(* Extracting command-line arguments by skipping unnecessary items *)

args = $ScriptCommandLine;  


NIons = ToExpression[args[[2]]];  (* First argument *)
Vrf := ToExpression[args[[3]]];  (* Second argument *)
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
ImportComsol["Electrode1_XY.txt",{3,3,0}]
ImportComsol["Electrode1_YZ.txt",{0,3,3}]
ImportComsol["Electrode1_ZX.txt",{3,0,3}]

ImportComsol["Electrode2_XY.txt",{3,3,0}]
ImportComsol["Electrode2_YZ.txt",{0,3,3}]
ImportComsol["Electrode2_ZX.txt",{3,0,3}]

ImportComsol["RF1_XY.txt",{3,3,0}]
ImportComsol["RF1_YZ.txt",{0,3,3}]
ImportComsol["RF1_ZX.txt",{3,0,3}]

ImportComsol["RF_ENorm_XY.txt",{3,3,0}];
ImportComsol["RF_ENorm_YZ.txt",{0,3,3}];
ImportComsol["RF_ENorm_ZX.txt",{3,0,3}]; 
]


DC1xy=Get["Electrode1_XY.m"];
DC1zy=Get["Electrode1_YZ.m"];
DC1zx=Get["Electrode1_ZX.m"];

DC2xy=Get["Electrode2_XY.m"];
DC2zy=Get["Electrode2_YZ.m"];
DC2zx=Get["Electrode2_ZX.m"];

RF1xy=Get["RF1_XY.m"];
RF1zy=Get["RF1_YZ.m"];
RF1zx=Get["RF1_ZX.m"];

(* Fields were exported in V/m*)
RFxy=Get["RF_ENorm_XY.m"];
RFzy=Get["RF_ENorm_YZ.m"];
RFzx=Get["RF_ENorm_ZX.m"];



voltageToCoupled=1/2 ({
 {2, 0, 0},
 {0, 1, 1},
 {0, 1, -1}
});


coupledToVoltage=Inverse[voltageToCoupled];


V:={5,0,0};
U:=coupledToVoltage . V;
appliedVoltages:=U;


(* Units constants *)
d = 1000;  (* mm/m *)
q = 1.602 10^(-19);

(* Atomic Constants *)
m = 171 1.67 10^(-27);

(* Experimental Constants *) 
\[CapitalOmega] = 2\[Pi] 38.6 10^6; (* RF drive frequency *)
 
DCMax:=2500;


xCtr=0.;
xRange=.5;
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


V := {V1,V2,V3};
DC1[t_]:=V[[1]];
DC2[t_]:=V[[2]];
RF1[t_]:=V[[3]];


(* Note: in v11, voltages applied to all endcap electrodes (DC1) and all center electrodes (RF1) in COMSOL *)
(*DCxy[x_,y_,z_,t_] :=DC1[t]*DC1xy[x+xoffset1,y+yoffset1,0]+DC2[t]*DC2xy[x+xoffset2,y+yoffset2,0]++RF1[t]*RF1xy[x+RF1xoffset,y+RF1yoffset,0] ;
DCzx[x_,y_,z_,t_] :=DC1[t]*DC1zx[x+xoffset1,0,z+zoffset1]+DC2[t]*DC2zx[x+xoffset2,0,z+zoffset2]+RF1[t]*RF1zx[x+RF1xoffset,0,z+RF1zoffset] ;
DCzy[x_,y_,z_,t_] :=DC1[t]*DC1zy[0,y+yoffset1,z+zoffset1]+DC2[t]*DC2zy[0,y+yoffset2,z+zoffset2]+RF1[t]*RF1zy[0,y+RF1yoffset,z+RF1zoffset];*)
Quiet[
DCxy[x_,y_,z_,t_] :=DC1[t]*DC1xy[x,y,0]+DC2[t]*DC2xy[x,y,0]+RF1[t]*RF1xy[x,y,0] ;
DCzx[x_,y_,z_,t_] :=DC1[t]*DC1zx[x,0,z]+DC2[t]*DC2zx[x,0,z]+RF1[t]*RF1zx[x,0,z] ;
DCzy[x_,y_,z_,t_] :=DC1[t]*DC1zy[0,y,z]+DC2[t]*DC2zy[0,y,z]+RF1[t]*RF1zy[0,y,z];
]


Quiet[
RFVxy[x_,y_,z_,t_] :=(q (Vrf^2 RFxy[x,y,0]^2))/(4 m (\[CapitalOmega]^2) );
RFVzx[x_,y_,z_,t_] :=(q ( Vrf^2 RFzx[x,0,z]^2))/(4 m (\[CapitalOmega]^2) );
RFVzy[x_,y_,z_,t_] :=(q ( Vrf^2 RFzy[0,y,z]^2))/(4 m (\[CapitalOmega]^2) );
]
(*RFVxoffset =(xcoord-xCtr)/.NMinimize[{RFxy[xcoord,yCtr,zCtr],xMin\[LessEqual] xcoord\[LessEqual] xMax},xcoord][[2]];
RFVyoffset =(ycoord-yCtr)/.NMinimize[{RFxy[xCtr,ycoord,zCtr],yMin\[LessEqual] ycoord\[LessEqual] yMax},ycoord][[2]];
RFVzoffset =(zcoord-zCtr)/.NMinimize[{RFzy[xCtr,yCtr,zcoord],zMin\[LessEqual] zcoord\[LessEqual] zMax},zcoord][[2]];
RFVxy[x_,y_,z_,t_] :=(q (Vrf^2 RFxy[x+RFVxoffset,y+RFVyoffset,0]^2))/(4 m (\[CapitalOmega]^2) )
RFVzx[x_,y_,z_,t_] :=(q (Vrf^2 RFzx[x+RFVxoffset,0,z+RFVzoffset]^2))/(4 m (\[CapitalOmega]^2) )
RFVzy[x_,y_,z_,t_] :=(q (Vrf^2 RFzy[0,y+RFVyoffset,z+RFVzoffset]^2))/(4 m (\[CapitalOmega]^2) )*)



Quiet[
Uxy[x_,y_,z_,t_] :=DCxy[x,y,0,t]+RFVxy[x,y,0,t];
Uzx[x_,y_,z_,t_] :=DCzx[x,0,z,t]+RFVzx[x,0,z,t];
Uzy[x_,y_,z_,t_] :=DCzy[0,y,z,t]+RFVzy[0,y,z,t];
]


Quiet[
Quad[x_,y_,a1_,b1_,a2_,b2_,c_]:=a2 x^2+b2 y^2+a1 x+b1 y+c x y;
QuadRot[x_,y_,\[Alpha]_,\[Beta]_,\[Theta]_]:=\[Alpha] (x Cos[\[Theta]]+y Sin[\[Theta]])^2+\[Beta] (-x Sin[\[Theta]]+y Cos[\[Theta]])^2;
]


Quiet[
NumPts=100;
yzFits=LinearModelFit[Flatten[Table[{y,z,Uzy[0,y,z,0]},{y,yMin,yMax,(yMax-yMin)/NumPts},{z,zMin,zMax,(zMax-zMin)/NumPts}],1],{y^2,y z,z^2,y,z},{y,z}];

xyFit=LinearModelFit[Flatten[Table[{x,y,Uxy[x,y,0,0]},{x,xMin,xMax,(xMax-xMin)/NumPts},{y,yMin,yMax,(yMax-yMin)/NumPts}],1],{x^2,x y,y^2,x,y},{x,y}];

yzFitParams=yzFits["BestFitParameters"];
xyFitParams=xyFit["BestFitParameters"];
]


Quiet[
Minimized=Minimize[yzFits[y,z],yMin<= y<= yMax&&zMin<= z<= zMax,{y,z}];
yzMins={y/.Minimized[[2,1]],z/.Minimized[[2,2]]};
]


Quiet[
\[Theta]x=ArcTan[(xyFitParams[[2]]-xyFitParams[[4]]),xyFitParams[[3]]]/2;
\[Theta]y=ArcTan[(yzFitParams[[2]]-yzFitParams[[4]]),yzFitParams[[3]]]/2;
\[Theta]z=\[Theta]y+\[Pi]/2;
]


Quiet[
PrincipleZaxis=Table[{yzMins[[1]]+ampl*Cos[\[Theta]z],yzMins[[2]]+ampl*Sin[\[Theta]z]},{ampl,zMin,zMax,0.001}];
PrincipleYaxis=Table[{yzMins[[1]]+ampl*Cos[\[Theta]y],yzMins[[2]]+ampl*Sin[\[Theta]y]},{ampl,yMin,yMax,0.001}];
]


(*
principalAxesContourPlot=ContourPlot[Uzy[0,y,z,0],{y,yMin,yMax},{z,zMin,zMax},ColorFunction->Hue,Contours->50,PlotPoints->{10,10},ContourLines->False,Frame->True,
ContourLabels-> None,
FrameLabel->{StyleForm["z (mm)",FontFamily->"Arial",FontSize->14],StyleForm["y (mm)",FontFamily->"Arial",FontSize->14]}];
principalAxesContourPlotXY=ContourPlot[Uxy[x,y,zCtr,0],{x,xMin,xMax},{y,yMin,yMax},ColorFunction->Hue,Contours->50,PlotPoints->{10,10},ContourLines->False,Frame->True,
ContourLabels-> None,AspectRatio->1/5,
FrameLabel->{StyleForm["x (mm)",FontFamily->"Arial",FontSize->14],StyleForm["y (mm)",FontFamily->"Arial",FontSize->14]}];
*)


Quiet[
\[Theta]x=\[Pi]/2;
FreqSecYaxis=1/(2\[Pi] )Sqrt[(d^2 q)2(yzFitParams[[2]]Cos[\[Theta]y]^2+yzFitParams[[3]]Sin[\[Theta]y]*Cos[\[Theta]y]+yzFitParams[[4]] Sin[\[Theta]y]^2)/m];
FreqSecZaxis=1/(2\[Pi] )Sqrt[(d^2 q)2(yzFitParams[[2]]Sin[\[Theta]y]^2-yzFitParams[[3]]Sin[\[Theta]y]*Cos[\[Theta]y]+yzFitParams[[4]] Cos[\[Theta]y]^2)/m];
FreqSecXaxis=1/(2\[Pi] )Sqrt[(d^2 q)2(xyFitParams[[2]]Sin[\[Theta]x]^2-xyFitParams[[3]]Sin[\[Theta]x]*Cos[\[Theta]x]+xyFitParams[[4]] Cos[\[Theta]x]^2)/m];
\[Kappa]=(2\[Pi]*FreqSecXaxis)^2*m*(200*10^-6)^2/(2*q*V[[1]]);
]
(*
Print["Radial secular frequency along principle axis Y (MHz): ", FreqSecYaxis*10^-6];
Print["Radial secular frequency along principle axis Z (MHz): ", FreqSecZaxis*10^-6];
Print["Radial secular frequency along X axis (MHz): ", FreqSecXaxis*10^-6];
Print["\[Kappa] = ", \[Kappa]];
*)


(*Export secular frequencies to a csv file for Python access *)
SecularFrequencies = {FreqSecXaxis*10^-6, FreqSecYaxis*10^-6, FreqSecZaxis*10^-6};
(*Print["Secular frequencies: ", SecularFrequencies];*)
Export["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\secular_frequencies_ren.csv", SecularFrequencies];
