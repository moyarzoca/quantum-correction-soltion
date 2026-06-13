
RootFolder = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{RootFolder, "tools", "solver.wl"}]];
Get[FileNameJoin[{RootFolder, "config.wl"}]];

ki = 1/1000;
kf = 10/1000;
paso = 1/3000;

(* One Frobenius solver — use the highest orden for both stages *)
frob = makeFrobeniusSolver[fineStage["orden"]];

(* ===== Stage 1: Find seed lambdas at ki for n=1..4 ===== *)
lambWindow = <|"center" -> 1/20, "width" -> 1/10, "Npoints" -> 20|>;

lambdaseeds = Table[
  result = RefineFindFirstLamb[n, ki, lambWindow, mediumStage, fineStage, frob];
  Print["n=", n, " → λ = ", result["lambda"]];
  result["lambda"],
  {n, 1, 4}
];
Abort[];

Print["** Done"];
(* ============================================================
	SEGUIMIENTO DEL MODO (tracking) Usa el resultado anterior como semilla para el siguiente k
============================================================*)
trackModo[nV_,kStart_,kEnd_,dk_,lambdaSeed_] := Module[
	{steps,result},
	steps=Range[kStart,kEnd,dk];
	result = FoldList[Function[{prev,kNow},
		Module[{lambdaPrev,lambdaNew},
			lambdaPrev = prev[[2]];(*NestWhileList reemplaza el For;FoldList acumula {k,lambda}*)
			lambdaNew = lambda/.FindRoot[SOL3[lambda, nV, kNow],{lambda, lambdaPrev},
									WorkingPrecision->prec, AccuracyGoal->prec/2, PrecisionGoal->prec/2
									];
			{kNow,lambdaNew}
		]
		],{kStart,lambdaSeed},(*valor inicial*)
			Rest[steps](*itera sobre los k siguientes*)
	];
	result  (*lista de pares {k,lambda}*)
];
(* ============================================================EJECUCIÓN============================================================*)
evolucionn1=trackModo[1,ki,kf,paso,lambdaseedn1];
evolucionn2=trackModo[2,ki,kf,paso,lambdaseedn2];
evolucionn3=trackModo[3,ki,kf,paso,lambdaseedn3];
evolucionn4=trackModo[4,ki,kf,paso,lambdaseedn4];

(*Extraer columnas*)
{kValsn1,lambdaValsn1}=Transpose[evolucionn1];
{kValsn2,lambdaValsn2}=Transpose[evolucionn2];
{kValsn3,lambdaValsn3}=Transpose[evolucionn3];
{kValsn4,lambdaValsn4}=Transpose[evolucionn4];

Save["evolucion_n1.m",evolucionn1]
Save["evolucion_n2.m",evolucionn2]
Save["evolucion_n3.m",evolucionn3]
Save["evolucion_n4.m",evolucionn4]
