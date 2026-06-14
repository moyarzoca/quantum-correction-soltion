
RootFolder = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{RootFolder, "tools", "solver.wl"}]];
Get[FileNameJoin[{RootFolder, "config.wl"}]];

kInitial = 1/100;
kStep = 10^-4;
NpointsInk = 40;

frob = makeFrobeniusSolver[fineStage["orden"]];

Print["** starting Stage 1 to find seed lambda"];
(* ===== Stage 1: Find seed lambdas at ki for n=1..4 ===== *)

lambdaseeds = <| |>;
Do[
	result = RefineFindFirstLamb[n, kInitial, lambWindow, mediumStage, fineStage, frob];
	Print["n=", n, " → λ = ", result["lambda"]];
	lambdaseeds[n]=result["lambda"];
,{n, 1, 2}];

kGrid = makeKGrid[kInitial, kStep, NpointsInk];

narrowLambWindow = lambWindow;
narrowLambWindow["width"] = 1/80;
narrowLambWindow["Npoints"] = 30;

evolutions = Table[
	Print["** Computing running of lambda : n = ", n];
	trackMode[n, kGrid, narrowLambWindow, mediumStage, lambdaseeds, frob]
  ,{n, 1, 1}
];

Print["** Done"];

