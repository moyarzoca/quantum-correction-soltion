(*Global variables*)

kInitial = 1/100;


(*Solver stages*)
mediumStage = <|
  "prec" -> 35, "orden" -> 5,  "z0" -> 10^-5, "z1Offset" -> 10^-6,
  "accuracyGoal" -> Automatic, "precisionGoal" -> Automatic, "maxSteps" -> Automatic
|>;

fineStage = <|
  "prec" -> 38, "orden" -> 7,  "z0" -> 10^-5, "z1Offset" -> 10^-7,
  "accuracyGoal" -> Automatic, "precisionGoal" -> Automatic, "maxSteps" -> 20000
|>;

(*
===========================================
		windows of serch of lambda
===========================================
*)
lambWindow = <|"center" -> 1/20, "width" -> 1/10, "Npoints" -> 30|>;

(* 
================================
	Mode tracker
================================
	- center is override by the last value of k
	- used in the mode tracker
*)

kStep = 1.5*10^-4;
NpointsInk = 40;

narrowLambWindow = <|"width" -> 1/80, "Npoints" -> 30|>;

(*Output configuration*)
outputConfig = <|
  "baseDir" -> "output",
  "evolutionsFile" -> "evolutions.wl",
  "fitsFile" -> "fits.wl",
  "logFile" -> "run.log"
|>;


