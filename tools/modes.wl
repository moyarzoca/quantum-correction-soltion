FindFirstLamb[n_, k_, lambWindow_, stage_, frob_] := Module[
  {prec, kp, lambMin, lambMax, lambInterval, Npoints, listlamb, solLamb, interp},

  prec = stage["prec"];
  kp = SetPrecision[k, prec];
  {lambMin, lambMax} = lambWindow["center"] + {-1, 1} * lambWindow["width"]/2;
  Npoints = lambWindow["Npoints"];
  lambInterval = Subdivide[lambMin, lambMax, Npoints];
  listlamb = Table[{lam, SOL3[lam, n, kp, frob, stage]}, {lam, lambInterval}];

  interp = Interpolation[listlamb];
  solLamb = FindRoot[interp[lamb], {lamb, lambWindow["center"]}];
  Print["  partial λ = ", (lamb /. solLamb), "  k = ", N[k], "  n = ", n];

  (*Print[Image[ListPlot[listlamb]]];*)

  Which[
	solLamb === {},
		Failure["NoRoot", <|"grid" -> listlamb|>],
	lambMin <= (lamb /. solLamb) <= lambMax,
		<|"lambda" -> lamb /. solLamb, "grid" -> listlamb|>,
	True,
		Print["OutOfRange!!!!"];
		Failure["OutOfRange", <|"grid" -> listlamb, "root" -> lamb /. solLamb|>]
  ]
];




RefineFindFirstLamb[n_, k_, lambWindow_, stage1_, stage2_, frob1_, frob2_: Automatic] := Module[
  {first, refinedWindow},
  first = FindFirstLamb[n, k, lambWindow, stage1, frob1];
  If[FailureQ[first], Return[first]];
  refinedWindow = <|
    "center" -> first["lambda"],
    "width" -> lambWindow["width"] / 6,
    "Npoints" -> lambWindow["Npoints"]
  |>;
  FindFirstLamb[n, k, refinedWindow, stage2, Replace[frob2, Automatic -> frob1]]
]; 

makeKGrid[kStart_, stepSize_, nPoints_] := Subdivide[kStart, kStart - stepSize*nPoints, nPoints];


trackMode[n_, kGrid_, NarrowLambWindow_, stage_, lambdaseeds_, frob_] := Module[
	{window = NarrowLambWindow, newlamb},
	lambdaSeed = lambdaseeds[n];
	listkLamb = {{kGrid[[1]], lambdaSeed}};

	window["center"] =  lambdaSeed;
	
	Do[
		newlamb = FindFirstLamb[n, newk, window, stage, frob]["lambda"];
		window["center"] =  newlamb;
		AppendTo[listkLamb, {newk, newlamb}]
	,{newk, Rest[kGrid]}];

	listkLamb
];
