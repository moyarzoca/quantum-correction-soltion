

P2[z_,lambda_, n_, k_] := (512*z*(z^2-3*z+3)^2*(z-1)^2*k^4+512*z*(z^2-3*z+3)^2*(z-1)^2*Sqrt[k^2+3]*k^3-768*(z-1)^2*(z^2-3)*(z^2-3*z+3)*z*k^2-768*z^2*(z^2-3*z+3)*(2*z-3)*(z-1)^2*Sqrt[k^2+3]*k+576*(z-1)^2*z^3*(2*z-3)^2);

P1[z_,lambda_, n_, k_] := ((512*(z-1))*(z^2-3*z+3)*(n*z^3-4*n*z^2+z^3+6*n*z-2*z^2-3*n-3)*k^4+(512*(z-1))*(z^2-3*z+3)*(n*z^3-4*n*z^2+z^3+6*n*z-2*z^2-3*n-3)*Sqrt[k^2+3]*k^3-(768*(z-1))*(n*z^5-4*n*z^4+z^5+3*n*z^3-5*z^4+9*n*z^2+12*z^3-18*n*z-18*z^2+9*n+9*z+9)*k^2-(768*(z-1))*z*(2*n*z^4-11*n*z^3+2*z^4+24*n*z^2-10*z^3-24*n*z+21*z^2+9*n-27*z+18)*Sqrt[k^2+3]*k+(576*(z-1))*z^2*(2*z-3)*(2*n*z^2-5*n*z+2*z^2+3*n-7*z+9));

P0[z_,lambda_, n_, k_] := ((128*n^2*z^5-1024*n^2*z^4+3584*n^2*z^3+512*n*z^4-2048*z^5-6912*n^2*z^2-3584*n*z^3+12032*z^4+512*lambda*z^2+6528*n^2*z+9216*n*z^2-28800*z^3-1536*lambda*z-2304*n^2-10752*n*z+36096*z^2+1536*lambda+4608*n-25728*z+11520)*k^4+(128*(n^2*z^5-8*n^2*z^4+28*n^2*z^3+4*n*z^4-16*z^5-54*n^2*z^2-28*n*z^3+94*z^4+4*lambda*z^2+51*n^2*z+72*n*z^2-225*z^3-12*lambda*z-18*n^2-84*n*z+282*z^2+12*lambda+36*n-201*z+90))*Sqrt[k^2+3]*k^3+(-192*n^2*z^5+960*n^2*z^4-1344*n^2*z^3+384*n*z^4+3072*z^5-1152*n^2*z^2-3840*n*z^3-11712*z^4+384*lambda*z^2+3456*n^2*z+13824*n*z^2+13824*z^3-2304*lambda*z-1728*n^2-20736*n*z-1152*z^2+3456*lambda+10368*n-13824*z+19008)*k^2-(192*(2*n^2*z^5-13*n^2*z^4+35*n^2*z^3+2*n*z^4-32*z^5-48*n^2*z^2-8*n*z^3+155*z^4+2*lambda*z^2+33*n^2*z-297*z^3-9*n^2+24*n*z+288*z^2-6*lambda-18*n-129*z-9))*Sqrt[k^2+3]*k-144*z*(-4*n^2*z^4+20*n^2*z^3-37*n^2*z^2+8*n*z^3+64*z^4+30*n^2*z-44*n*z^2-244*z^3+8*lambda*z-9*n^2+72*n*z+333*z^2-12*lambda-36*n-150*z-27));

Equ[lambda_, n_, k_] := P1[z, lambda, n, k]*D[F[z], z]+P2[z, lambda, n, k]*D[F[z],z,z]+P0[z, lambda, n, k]*F[z];

makeFrobeniusSolver[orden_] := Module[
	{R, Eq, RawEqs, Eqs, minusb, m, b, solCoeffs, evalCoefs, solu, Dsolu},
	R = 1 + Sum[a[i] z^(i), {i, 1, orden}] + O[z]^(orden + 1);
	Eq = P2[z, lambda, n, k]*D[R, z, z] + P1[z, lambda, n, k]*D[R, z] + P0[z, lambda, n, k]*R;
	Do[RawEqs[j] = SeriesCoefficient[Eq, j], {j, 0, orden}];
	Eqs = Thread[Array[RawEqs, orden - 1, 0] == ConstantArray[0, orden - 1]];
	{minusb, m} = CoefficientArrays[Eqs, Array[a, orden - 1]];
	b = -minusb;
	solCoeffs = LinearSolve[m, b];
	evalCoefs = Thread[Array[a, orden - 1] -> solCoeffs];
	solu = 1 + Sum[a[i] z^(i), {i, 1, orden - 1}] /. evalCoefs;
	Dsolu = D[solu, z];
	<|
	"solu"  -> Function[{z0, lam, nIn, kIn}, Evaluate[solu  /. {z -> z0, lambda -> lam, n -> nIn, k -> kIn}]],
	"Dsolu" -> Function[{z0, lam, nIn, kIn}, Evaluate[Dsolu /. {z -> z0, lambda -> lam, n -> nIn, k -> kIn}]]
	|>
];
