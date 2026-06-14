
RootFolder = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{RootFolder, "tools", "equation.wl"}]];
Get[FileNameJoin[{RootFolder, "tools", "solver.wl"}]];
Get[FileNameJoin[{RootFolder, "tools", "modes.wl"}]];
Get[FileNameJoin[{RootFolder, "config.wl"}]];

Get[FileNameJoin[{RootFolder, "tools", "export.wl"}]];
runCfg = setupRunDirectory[RootFolder, outputConfig];
runDir = runCfg["dir"];
logStream = OpenWrite[FileNameJoin[{runDir, outputConfig["logFile"]}]];
AppendTo[$Output, logStream];
CopyFile[FileNameJoin[{RootFolder, "config.wl"}], FileNameJoin[{runDir, "config.wl"}]];


frob = makeFrobeniusSolver[fineStage["orden"]];

Print["Stage 1 to find seed lambda"];
(* ===== Stage 1: Find seed lambdas at ki for n=1..4 ===== *)

lambdaseeds = <| |>;
clock = AbsoluteTiming[
	Do[
		result = RefineFindFirstLamb[n, kInitial, lambWindow, mediumStage, fineStage, frob];
		Print["n=", n, " → λ = ", result["lambda"]];
		lambdaseeds[n]=result["lambda"];
	,{n, 1, 2}];
];

Print["Seeds computed  in ", clock[[1]], " sec."];

kGrid = makeKGrid[kInitial, kStep, NpointsInk];

evolutionsN = <||>;
fitsN = <||>;
clock = AbsoluteTiming[
	Do[
		Print["** Computing running of lambda : n = ", n];
		evolutionsN[n] = trackMode[n, kGrid, narrowLambWindow, mediumStage, lambdaseeds, frob];

		Print["fit n = ",n , ":  ",fitsN[n] = Fit[evolutionsN[n], {1, x}, x]];
	, {n, 1, 2}
	];
]

Print["Running of the modes computed in ", clock[[1]], " sec."];

Export[FileNameJoin[{runDir, outputConfig["evolutionsFile"]}], evolutionsN, "WL"];
Export[FileNameJoin[{runDir, outputConfig["fitsFile"]}], fitsN, "WL"];

Print["** Done"];

(*Close log stream*)
Close[logStream];

