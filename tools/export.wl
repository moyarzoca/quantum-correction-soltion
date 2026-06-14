
setupRunDirectory[baseFolder_, config_] := Module[
  {outputBase, nums, runNum, runDir},

  outputBase = FileNameJoin[{baseFolder, config["baseDir"]}];

  If[DirectoryQ[outputBase], Null, CreateDirectory[outputBase]];

  nums = Cases[FileNames[DigitCharacter .., outputBase], s_String :> ToExpression[FileNameTake[s]]];
  runNum = If[nums === {}, 1, Max[nums] + 1];
  runDir = FileNameJoin[{outputBase, IntegerString[runNum, 10, 3]}];
  CreateDirectory[runDir];
  <|"dir" -> runDir|>
];
