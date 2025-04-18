(* ::Package:: *)

hbar=6.626/(2\[Pi])*10^-34;
qe=1.602*10^-19;
mp=1.67*10^-27;
mYb=171*mp;
\[Epsilon]0=8.85*10^-12;
\[CapitalOmega]t=2\[Pi]*21.033*10^6; (*Drive Frequency*)
cf=qe/(mYb*\[CapitalOmega]t^2); (*Common factor in Mathieu parameters*)
Cou=qe^2/(4\[Pi]*\[Epsilon]0);
ct=4/(mYb*\[CapitalOmega]t^2);
kB=1.380649*\!\(\*SuperscriptBox[\(10\), \(\[Minus]23\)]\);
\[CapitalDelta]k=Sqrt[2]*2\[Pi]/(355*10^(-9)); (* The momentum kick difference between the 2 Raman beams. Assumes an angle of \[Pi]/2 between them. AK *)
(*R= (hbar*\[CapitalDelta]k^2)/(2*mYb);*)





IonPositions[NIons_,\[Omega]x_,\[Omega]y_,\[Omega]z_,ShowInfo_]:=Module[
{x0,y0,z0,maxT,Vinitial,friction,Fx,Fy,Fz,sol,xf,yf,zf},
x0=Table[RandomReal[{-20*10^(-6),20*10^(-6)}],{i,NIons}];
y0=Table[RandomReal[{-20*10^(-6),20*10^(-6)}],{i,NIons}];
z0=Table[RandomReal[{-20*10^(-6),20*10^(-6)}],{i,NIons}];
maxT=10*2*10^(-1); 
Vinitial=0;
friction= 300*10^-19;

Clear[Fx,Fy,Fz,t];
Table[Fx[i][t]=-mYb*\[Omega]x^2*XPos[i][t]+Cou*(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(i - 1\)]\(\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\)/
\*SuperscriptBox[\((
\*SuperscriptBox[\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\), \(2\)])\), \(3/2\)]\)\) +\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = i + 1\), \(NIons\)]\(\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\)/
\*SuperscriptBox[\((
\*SuperscriptBox[\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\), \(2\)])\), \(3/2\)]\)\) )-friction*XPos[i]'[t],{i,NIons}];
Table[Fy[i][t]=-mYb*\[Omega]y^2*YPos[i][t]+Cou*(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(i - 1\)]\(\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\)/
\*SuperscriptBox[\((
\*SuperscriptBox[\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\), \(2\)])\), \(3/2\)]\)\) +\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = i + 1\), \(NIons\)]\(\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\)/
\*SuperscriptBox[\((
\*SuperscriptBox[\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\), \(2\)])\), \(3/2\)]\)\) )-friction*YPos[i]'[t],{i,NIons}];
Table[Fz[i][t]=-mYb*\[Omega]z^2*ZPos[i][t]+Cou*(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(i - 1\)]\(\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\)/
\*SuperscriptBox[\((
\*SuperscriptBox[\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\), \(2\)])\), \(3/2\)]\)\) +\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = i + 1\), \(NIons\)]\(\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\)/
\*SuperscriptBox[\((
\*SuperscriptBox[\((\(XPos[i]\)[t] - \(XPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(YPos[i]\)[t] - \(YPos[j]\)[t])\), \(2\)] + 
\*SuperscriptBox[\((\(ZPos[i]\)[t] - \(ZPos[j]\)[t])\), \(2\)])\), \(3/2\)]\)\) )-friction*ZPos[i]'[t],{i,NIons}];

XVars=Table[XPos[i],{i,NIons}];
YVars=Table[YPos[i],{i,NIons}];
ZVars=Table[ZPos[i],{i,NIons}];
x0Eqns=Table[XPos[i][0]==x0[[i]],{i,NIons}];
y0Eqns=Table[YPos[i][0]==y0[[i]],{i,NIons}];
z0Eqns=Table[ZPos[i][0]==z0[[i]],{i,NIons}];
Vx0Eqns=Table[XPos[i]'[0]==Vinitial,{i,NIons}];
Vy0Eqns=Table[YPos[i]'[0]==Vinitial,{i,NIons}];
Vz0Eqns=Table[ZPos[i]'[0]==Vinitial,{i,NIons}];
FxEqns=Table[Fx[i][t] == mYb*XPos[i]''[t], {i, NIons}];
FyEqns=Table[Fy[i][t] == mYb*YPos[i]''[t], {i, NIons}];
FzEqns=Table[Fz[i][t] == mYb*ZPos[i]''[t], {i, NIons}];
(*PrintTemporary["Calculating..."];*)
MyTimer=AbsoluteTime[];
sol=NDSolve[Join[FxEqns,FyEqns,FzEqns,Vx0Eqns,Vy0Eqns,Vz0Eqns,x0Eqns,y0Eqns,z0Eqns],Join[XVars,YVars,ZVars],{t,0,maxT},MaxSteps->10000000,AccuracyGoal->8,PrecisionGoal->4];
If[ShowInfo==True,Print["Calculation took ", AbsoluteTime[]-MyTimer, " s"]];

xf=Flatten[Table[XPos[i][maxT]/.sol,{i,NIons}]];
yf=Flatten[Table[YPos[i][maxT]/.sol,{i,NIons}]];
zf=Flatten[Table[ZPos[i][maxT]/.sol,{i,NIons}]];

xyzfriffle=Partition[Flatten@Transpose[{xf,yf,zf}],3];
pm1 = Automatic;
(*If[
ShowInfo==True,
Print[ListPlot[Transpose[{xf,yf}]*10^6,Frame->True,Axes->False,PlotRange->{{-20,20},{-20,20}},AspectRatio->1,BaseStyle->{FontFamily->"Helvetica",FontSize->17},PlotMarkers->pm1,FrameLabel->{"x position (\[Mu]m)","y position (\[Mu]m)"},ImageSize->489.6]];
Print[ListPlot[Transpose[{zf,yf}]*10^6,Frame->True,Axes->False,PlotRange->{{-10,10},{-10,10}},AspectRatio->1,BaseStyle->{FontFamily->"Helvetica",FontSize->17},PlotMarkers->pm1,FrameLabel->{"z position (\[Mu]m)","y position (\[Mu]m)"},ImageSize->489.6/2]];
Print[Li "Z:\\Users\\Pete\\equil_pos3.csv"stPlot[Transpose[{zf,xf}]*10^6,Frame->True,Axes->False,PlotRange->{{-10,10},{-10,10}},AspectRatio->1,BaseStyle->{FontFamily->"Helvetica",FontSize->17},PlotMarkers->pm1,FrameLabel->{"z position (\[Mu]m)","x position (\[Mu]m)"},ImageSize->489.6/2]];
];*)
FinalPositions=Transpose[{xf,yf,zf}];
Export["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\equil_positions.csv", Chop@FinalPositions]
FinalPositions
]


(* Extracting command-line arguments by skipping unnecessary items *)

args = $ScriptCommandLine;  (* Skip the first three elements: -wlbanner, -script, and the script name *)

(* Now `args` should contain the actual parameters you need *)
NIons = ToExpression[args[[2]]];  (* First argument *)
wx = ToExpression[args[[3]]];  (* Second argument *)
wy = ToExpression[args[[4]]];  (* Third argument *)
wz = ToExpression[args[[5]]];  (* Fourth argument *)
ShowInfo = ToExpression[args[[5]]];   (* Convert ShowInfo as needed *)
(* Print to verify *)
(*
Print["args: ", args]
Print["NIons: ", NIons];
Print["omega_x: ", wx];
Print["omega_y: ", wy];
Print["omega_z: ", wz];
*)


IonPositions[NIons, wx, wy, wz, ShowInfo];

