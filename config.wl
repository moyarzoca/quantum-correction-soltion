mediumStage = <|
  "prec" -> 35, "orden" -> 5,  "z0" -> 10^-5, "z1Offset" -> 10^-6,
  "accuracyGoal" -> Automatic, "precisionGoal" -> Automatic, "maxSteps" -> Automatic
|>;

fineStage = <|
  "prec" -> 38, "orden" -> 7,  "z0" -> 10^-5, "z1Offset" -> 10^-7,
  "accuracyGoal" -> Automatic, "precisionGoal" -> Automatic, "maxSteps" -> 20000
|>;

lambWindow = <|"center" -> 1/20, "width" -> 1/10, "Npoints" -> 30|>;
