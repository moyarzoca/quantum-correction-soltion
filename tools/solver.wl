
SOL[lambda_, n_, k_, frob_, stage_] := Module[
  {prec, z0, z1, acc, precGoal, maxStp, lambdap, np, kp, bc, sys},

  prec = stage["prec"];
  z0 = SetPrecision[stage["z0"], prec];
  z1 = 1 - SetPrecision[stage["z1Offset"], prec];
  acc = Replace[stage["accuracyGoal"], Automatic :> prec/2];
  precGoal = Replace[stage["precisionGoal"], Automatic :> prec/2];
  maxStp = Replace[stage["maxSteps"], Automatic :> 10000];
  {lambdap, np, kp} = SetPrecision[{lambda, n, k}, prec];

	sys = {
		Equ[lambdap, np, kp] == 0,
		F[z0]  == frob["solu"][z0, lambdap, np, kp],
		F'[z0] == frob["Dsolu"][z0, lambdap, np, kp]
	};

  NDSolve[sys, F, {z, z0, z1},
    WorkingPrecision -> prec - 4, AccuracyGoal -> acc,
    PrecisionGoal -> precGoal, MaxSteps -> maxStp,
    Method -> "StiffnessSwitching"]
];

SOL3[lambda_, n_, k_, frob_, stage_] := Module[{res},
  prec = stage["prec"];
  z1 = 1 - SetPrecision[stage["z1Offset"], prec];
  res = AbsoluteTiming[F[z1] /. SOL[lambda, n, k, frob, stage][[1]]];
  If[VERBO === True, Print[res[[1]]]];
  res[[2]]
];



