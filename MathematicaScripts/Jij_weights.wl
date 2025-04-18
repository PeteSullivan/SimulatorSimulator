(* ::Package:: *)

(* Define helper functions *)
utmx[mat_] := Flatten[UpperTriangularize[mat, 1]]; (* Upper triangular, exclude diagonal *)
nodiag[mat_] := mat - DiagonalMatrix[Diagonal[mat]]; (* Zero out diagonal *)

(* Main function - now takes direct matrix inputs *)
CalculateFromMTXFiles[] := Module[
  {ions, Jall, bmx, Jijk, JijkMW1, JijkMW, possiblem, mmttable, mmijm, 
   Jijkjoin, Jijkjm, ml, ctable1, ctable, jcm1, jcm, JDes, NormalModeEigVecs, bstable, lpout},
  
  (* Import matrices directly from files *)
  JDes = Import["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\JDes.mtx", "MTX"];
  NormalModeEigVecs = Import["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\normal_modes.mtx", "MTX"];
  
  (* Debug prints - remove MatrixForm as it's only for display *)
  Print["JDes dimensions: ", Dimensions[JDes]];
  Print["NormalModeEigVecs dimensions: ", Dimensions[NormalModeEigVecs]];
  
  ions = Length[JDes];
  
  (* Calculate Jijk from normal modes *)
  bmx = NormalModeEigVecs // Chop;
  Jijk = Table[nodiag[Transpose[{bmx[[i]]}] . {bmx[[i]]}], {i, ions - 1}] // Chop;
  
  (* Build matrices for linear programming *)
  JijkMW1 = ions*Transpose[Table[utmx[Jijk[[i]]], {i, 1}]];
  JijkMW = ions*Transpose[Table[utmx[Jijk[[i]]], {i, ions - 1}]];
  
  (* Generate possible magnetization patterns *)
  possiblem = {ConstantArray[1, ions]};
  seed = possiblem[[1]]; seed[[1]] = -1;
  possiblem = Join[possiblem, Permutations[seed]];
  
  (* Prepare constraints *)
  mmttable = DeleteDuplicates[Table[Transpose[{possiblem[[i]]}] . {possiblem[[i]]}, {i, Length[possiblem]}]];
  mmijm = Transpose[utmx[#] & /@ mmttable];
  mmijm = Transpose[Join[Transpose[mmijm], -Transpose[mmijm]]];
  
  Jijkjoin = ions*Flatten[Table[(#*Jijk[[i]] & /@ mmttable), {i, ions - 1}], 1];
  Jijkjm = Transpose[utmx[#] & /@ Jijkjoin];
  Jijkjm = Transpose[Join[Transpose[Jijkjm], -Transpose[Jijkjm]]];
  
  (* Set up linear programming *)
  ml = Dimensions[mmijm][[2]];
  ctable1 = Table[Table[ToExpression["ci" <> ToString[k] <> "j" <> ToString[m]], {m, ml}], {k, 1}];
  ctable = Table[Table[ToExpression["ci" <> ToString[k] <> "j" <> ToString[m]], {m, ml}], {k, ions - 1}];
  
  jcm1 = (JijkMW1 . ctable1)*mmijm;
  jcm = (JijkMW . ctable)*mmijm;
  
  (* Solve *)
  bstable = Table[{utmx[JDes][[i]], 0}, {i, Length[utmx[JDes]]}];
  lpout = LinearProgramming[ConstantArray[1, Length[Flatten[ctable]]], Jijkjm, bstable] // Chop;
  
  (* Return results *)
  <|
    "Weights" -> lpout,
    "TotalWeights" -> Total[lpout],
    "JDesNorm" -> Total[Flatten[Abs[JDes]]]/2,
    "JDes" -> JDes,
    "NormalModeEigVecs" -> NormalModeEigVecs
  |>
]

(* Execute the function *)
result = CalculateFromMTXFiles[];
Print["Weights: ", result["Weights"]];
Print["Total weights: ", result["TotalWeights"]];

Export["Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\Jij_weights.csv",result["Weights"]];
