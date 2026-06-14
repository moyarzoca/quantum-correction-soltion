
RootFolder = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{RootFolder, "tools", "equation.wl"}]];
Get[FileNameJoin[{RootFolder, "tools", "solver.wl"}]];
Get[FileNameJoin[{RootFolder, "tools", "modes.wl"}]];
Get[FileNameJoin[{RootFolder, "config.wl"}]];

kInitial = 1/100;
kStep = 10^-4;
NpointsInk = 40;

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

Print["Seeds computed  in ", clock[[1]], "sec."];

kGrid = makeKGrid[kInitial, kStep, NpointsInk];

narrowLambWindow = lambWindow;
narrowLambWindow["width"] = 1/80;
narrowLambWindow["Npoints"] = 30;

clock = AbsoluteTiming[
	evolutions = Table[
		Print["** Computing running of lambda : n = ", n];
		trackMode[n, kGrid, narrowLambWindow, mediumStage, lambdaseeds, frob]
	  ,{n, 1, 1}
	];
]

Print["Running of the modes computed in ", clock[[1]], "sec."];

Print["** Done"];

