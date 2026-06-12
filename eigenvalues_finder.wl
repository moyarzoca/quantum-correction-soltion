
{orden = 10, prec = 30, z0 = 10^-6, z1 = 1 - 10^-6};

R = 1 + Sum[a[i] z^(i), {i, 1, orden}] + O[z]^(orden + 1);

P2[z_] := (512*z*(z^2-3*z+3)^2*(z-1)^2*k^4+512*z*(z^2-3*z+3)^2*(z-1)^2*Sqrt[k^2+3]*k^3-768*(z-1)^2*(z^2-3)*(z^2-3*z+3)*z*k^2-768*z^2*(z^2-3*z+3)*(2*z-3)*(z-1)^2*Sqrt[k^2+3]*k+576*(z-1)^2*z^3*(2*z-3)^2);

P1[z_] := ((512*(z-1))*(z^2-3*z+3)*(n*z^3-4*n*z^2+z^3+6*n*z-2*z^2-3*n-3)*k^4+(512*(z-1))*(z^2-3*z+3)*(n*z^3-4*n*z^2+z^3+6*n*z-2*z^2-3*n-3)*Sqrt[k^2+3]*k^3-(768*(z-1))*(n*z^5-4*n*z^4+z^5+3*n*z^3-5*z^4+9*n*z^2+12*z^3-18*n*z-18*z^2+9*n+9*z+9)*k^2-(768*(z-1))*z*(2*n*z^4-11*n*z^3+2*z^4+24*n*z^2-10*z^3-24*n*z+21*z^2+9*n-27*z+18)*Sqrt[k^2+3]*k+(576*(z-1))*z^2*(2*z-3)*(2*n*z^2-5*n*z+2*z^2+3*n-7*z+9));

P0[z_] := ((128*n^2*z^5-1024*n^2*z^4+3584*n^2*z^3+512*n*z^4-2048*z^5-6912*n^2*z^2-3584*n*z^3+12032*z^4+512*lambda*z^2+6528*n^2*z+9216*n*z^2-28800*z^3-1536*lambda*z-2304*n^2-10752*n*z+36096*z^2+1536*lambda+4608*n-25728*z+11520)*k^4+(128*(n^2*z^5-8*n^2*z^4+28*n^2*z^3+4*n*z^4-16*z^5-54*n^2*z^2-28*n*z^3+94*z^4+4*lambda*z^2+51*n^2*z+72*n*z^2-225*z^3-12*lambda*z-18*n^2-84*n*z+282*z^2+12*lambda+36*n-201*z+90))*Sqrt[k^2+3]*k^3+(-192*n^2*z^5+960*n^2*z^4-1344*n^2*z^3+384*n*z^4+3072*z^5-1152*n^2*z^2-3840*n*z^3-11712*z^4+384*lambda*z^2+3456*n^2*z+13824*n*z^2+13824*z^3-2304*lambda*z-1728*n^2-20736*n*z-1152*z^2+3456*lambda+10368*n-13824*z+19008)*k^2-(192*(2*n^2*z^5-13*n^2*z^4+35*n^2*z^3+2*n*z^4-32*z^5-48*n^2*z^2-8*n*z^3+155*z^4+2*lambda*z^2+33*n^2*z-297*z^3-9*n^2+24*n*z+288*z^2-6*lambda-18*n-129*z-9))*Sqrt[k^2+3]*k-144*z*(-4*n^2*z^4+20*n^2*z^3-37*n^2*z^2+8*n*z^3+64*z^4+30*n^2*z-44*n*z^2-244*z^3+8*lambda*z-9*n^2+72*n*z+333*z^2-12*lambda-36*n-150*z-27));

Eq = P2[z] D[R, z, z] + P1[z] D[R, z] + P0[z] R;

Do[eqs[j] = SeriesCoefficient[Eq, j], {j, 0, orden}]

{minusb, m} = 
  CoefficientArrays[
   Thread[Array[eqs, orden - 1, 0] == ConstantArray[0, orden - 1]], 
   Array[a, orden - 1]];

b=-minusb;
sol=LinearSolve[m,b];
coefs=Thread[Array[a,orden-1]->sol];
solu=1+Sum[a[i]z^(i),{i,1,orden-1}] /. coefs;
Dsolu = D[solu,z];
bcsolu[lambda_?NumericQ,n_?NumericQ, k_?NumericQ]=SetPrecision[solu/.z-> z0,prec];
Dbcsolu[lambda_?NumericQ,n_?NumericQ, k_?NumericQ]=SetPrecision[Dsolu/.z-> z0,prec];
Equ[lambda_?NumericQ,n_?NumericQ, k_?NumericQ]=P1[z]D[F[z],z]+P2[z]D[F[z],z,z]+P0[z]F[z];
SOL[lambda_?NumericQ,n_?NumericQ, k_?NumericQ]:=Module[{lambdap, np, kp,bc,sys},{lambdap,np,kp}= SetPrecision[{lambda,n,k},prec];bc= {F[z0] ==bcsolu[lambdap,np,kp], F'[z0] ==Dbcsolu[lambdap,np,kp]};sys={Equ[lambdap,np,kp]==0, bc[[1]],bc[[2]]};NDSolve[ sys,F,{z,z0,z1},WorkingPrecision->prec, AccuracyGoal->prec/2, PrecisionGoal->prec/2, MaxStepSize->10^-3,Method->"StiffnessSwitching"]];
SOL3[lambda_?NumericQ,n_?NumericQ, k_?NumericQ]:=Evaluate[F[z]/.SOL[lambda,n,k]][[1]]/.z->z1
(*
===============
Seeds
===============
*)
BuscarModos[nVal_,kVal_,lambdamin_,lambdamax_,nPuntos_:200] := Module[
	{nP,kP,vals,signChanges,seeds,modos},
	{nP,kP}=SetPrecision[{nVal,kVal},prec];
	vals=Table[{lambda,SOL3[lambda,nVal,kVal]},{lambda,Subdivide[lambdamin,lambdamax,nPuntos]}];
	signChanges=Pick[vals[[;;-2,1]],Negative[vals[[;;-2,2]]*vals[[2;;,2]]]]; (*Detecta cambios de signo*)
	modos=Table[lambda/.FindRoot[SOL3[lambda,nVal,kVal],{lambda,seed},WorkingPrecision->prec,AccuracyGoal->10],{seed,signChanges}];(*Refina cada cero*)
	modos
];



ki = SetPrecision[1/1000,prec];
kf = SetPrecision[10/1000,prec];
paso = SetPrecision[1/3000,prec];

(*Semilla inicial - busca el modo en ki con barrido grueso*)

lambdaseedn1 = First@BuscarModos[1,ki,0,3/100];
lambdaseedn2 = First@BuscarModos[2,ki,0,3/100];
lambdaseedn3 = First@BuscarModos[3,ki,0,3/100];
lambdaseedn4 = First@BuscarModos[4,ki,0,3/100];
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
