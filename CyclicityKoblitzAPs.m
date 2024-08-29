/////////////////////////////////////////////////
// All of our code was run on a machine with the
// following specifications.
// CPU: Apple M3 Pro
// Memory: 18 GB
// OS: macOS Sonoma Version 14.5
// Magma Version: V2.28-5
/////////////////////////////////////////////////

// We first load "OpenImage" by D. Zywina, which 
// is available online at 
// https://github.com/davidzywina/OpenImage
// * Required for KoblitzAP and CyclicityAP.
// * Requires Magma V2.27 or later.
ChangeDirectory("/Users/Jacob/OpenImage-master");
load "main/FindOpenImage.m";

// Given groups G and H, output the homomorphism
// pi: G -> H defined by coercing elements from
// G to H, provided this is possible.
Red := function(G,H)
    return hom<G -> H | [<g, H!g> : g in Generators(G)]>;
end function;

// Given a positive integer n, output the
// radical of n.
Rad := function(n)
    return &*PrimeDivisors(n);
end function;

// Computes L as defined in Eq. (6).
ComputeL := function(n,mE)
    return &*[ell^Max(1, Valuation(n, ell)) : ell in PrimeDivisors(mE)];
end function;

// Given a Serre curve E, output the adelic
// level m_E of E, computed using Prop. 2.4.
AdelicLevelSerreCurve := function(E)
    DeltaPrime := SquareFree(Integers() ! Discriminant(IntegralModel(E)));
    if DeltaPrime mod 4 eq 1 then
        return 2 * Abs(DeltaPrime);
    else
        return 4 * Abs(DeltaPrime);
    end if;
end function;

// Computes tau as in Definition 5.1.
tau := function(DeltaPrime,k)
    if DeltaPrime mod 4 eq 1 then
        return -1;
    end if;
    if DeltaPrime mod 4 eq 3 then
        if k mod 4 eq 1 then
            return -1;
        else
            return 1;
        end if;
    end if;
    if DeltaPrime mod 8 eq 2 then
        if k mod 8 in [1,7] then
            return -1;
        else
            return 1;
        end if;
    end if;
    if DeltaPrime mod 8 eq 6 then
        if k mod 8 in [1,3] then
            return -1;
        else
            return 1;
        end if;
    end if;
end function;

// Computes tau^cyc as in Definition 5.1.
tauCyc := function(L,DeltaPrime,n,k)
    Lodd := Integers() ! (L / 2^Valuation(L,2));
    N1 := &*Include([-1 : ell in PrimeDivisors(Lodd) | n mod ell ne 0],1);
    N2 := &*Include([KroneckerSymbol(k,ell) : ell in PrimeDivisors(GCD(n,Lodd)) | (k-1) mod ell ne 0],1);
    return tau(DeltaPrime,k) * N1 * N2;
end function;

// Computes tau^prime as in Definition 5.1.
tauPrime := function(L,DeltaPrime,n,k)
    Lodd := Integers() ! (L / 2^Valuation(L,2));
    N := &*Include([KroneckerSymbol(k,ell) : ell in PrimeDivisors(GCD(n,Lodd)) | (k-1) mod ell ne 0],1);
    return -tau(DeltaPrime,k) * N;
end function;

// Given n, output an associative array with the
// values of C^prime_{n,k}, as in Eq. 40, for all
// units k modulo n. The optional parameter B 
// specifies the largest prime used in the
// infinite product.
AvgKoblitzAP := function(n:B:=10^5)
	N1 := &*[1 - (ell^2-ell-1)/((ell-1)^3*(ell+1)) : ell in PrimesUpTo(B) | n mod ell ne 0];
	A := AssociativeArray();
	for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
		N2 := &*Include([1 - (ell^2+ell)/((ell^2-1)*(ell^2-ell)) : ell in PrimeDivisors(n) | (k-1) mod ell ne 0],1);
		N3 := &*Include([1 - ell/((ell^2-1)*(ell^2-ell)) : ell in PrimeDivisors(n) | (k-1) mod ell eq 0],1);
		A[k] := N1*N2*N3/EulerPhi(n);
	end for;
	return A;
end function;

// Given a non-CM elliptic curve E and positive 
// integer n, output an associative array with the 
// values of delta^prime_{E,n,k}(L),  as in Eq. 49, 
// for all units k modulo n. The optional parameter 
// B specifies the largest prime used in the 
// infinite product.
DeltaKoblitzAP := function(E,n:B:=10^5)
	G := FindOpenImage(E);
	mE := #BaseRing(G);
	R := Rad(mE);
	L := ComputeL(n,mE);
	pi1 := Red(GL(2,Integers(LCM(mE,L))), GL(2,Integers(mE)));
	Ga := G @@ pi1;
	pi2 := Red(GL(2,Integers(LCM(mE,L))),GL(2,Integers(L)));
	GLft := pi2(Ga);
	C := ConjugacyClasses(GLft);
	S := {* (Integers(GCD(L,n)) ! Determinant(c[3]))^^c[2] : c in C | IsUnit(Determinant(Matrix(One(GLft)) - Matrix(c[3]))) *};
	A := AssociativeArray();
	for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
		A[k] := Multiplicity(S,k)/#GLft;
	end for;
	return A, R;
end function;

// Computes delta^prime_{E,n,k} via Lemma 5.2 in 
// the simple case that ell divides n yet ell 
// does not divide m_E.
DeltaKoblitzAPSimple := function(n,k,ell)
    alpha := Valuation(n,ell);
    if k mod ell eq 1 then
          return 1/(ell^alpha - ell^(alpha-1)) * (1 - ell/((ell^2-1)*(ell^2-ell)));
        else
          return 1/(ell^alpha - ell^(alpha-1)) * (1 - (ell^2+ell)/((ell^2-1)*(ell^2-ell)));
    end if;
end function;

// Given a non-CM elliptic curve E and positive 
// integer n, use Prop. 4.4 to output an 
// associative array with the values of 
// C^prime_{E,n,k} for all units k modulo n.
// The optional parameter B specifies the largest 
// prime used in the infinite product.
KoblitzAP := function(E,n:B:=10^5)
	DeltaBiases, R := DeltaKoblitzAP(E,n:B:=B);
	D := &*[ 1 - 1/ell : ell in PrimeDivisors(R)];
	N1 := &*[1 - (ell^2-ell-1)/((ell-1)^3*(ell+1)) : ell in PrimesUpTo(B) | (n*R) mod ell ne 0];
	A := AssociativeArray();
	for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
		N2 := &*Include([DeltaKoblitzAPSimple(n,k,ell) : ell in PrimeDivisors(n) | R mod ell ne 0],1);
		A[k] := RealField() ! (DeltaBiases[k]/D * N1*N2);
	end for;
	return A;
end function;

// Given a Serre curve E and positive 
// integer n, use Thm. 1.7 to output an 
// associative array with the values of 
// C^prime_{E,n,k} for all units k modulo n.
// The optional parameter B specifies the largest 
// prime used in the infinite product.
SerreCurveKoblitzAP := function(E,n:B:=10^5)
    mE := AdelicLevelSerreCurve(E);
    L := ComputeL(n,mE);
    AvgKob := AvgKoblitzAP(n:B:=B);
    A := AssociativeArray();
    if L mod mE ne 0 then
        for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
            A[k] := RealField() ! AvgKob[k];
        end for;
    else
        P := &*Include([1/(ell^3-2*ell^2-ell+3) : ell in PrimeDivisors(L) | (2*n) mod ell ne 0],1);
        for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
            DeltaPrime := SquareFree(Integers() ! Discriminant(IntegralModel(E)));
            tauValue := tauPrime(L,DeltaPrime,n,k);
            A[k] := RealField() ! (AvgKob[k] * (1 + tauValue * P));
        end for;
    end if;
    return A;
end function;

// Given n, output an associative array with the
// values of C^cyc_{n,k}, as in Eq. 21, for all
// units k modulo n. The optional parameter B 
// specifies the largest prime used in the
// infinite product.
AvgCyclicityAP := function(n:B:=10^5)
	N1 := &*[1 - 1 / ((ell^2-1)*(ell^2-ell)) : ell in PrimesUpTo(B) | n mod ell ne 0];
	A := AssociativeArray();
	for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
		N2 := &*Include([1 - (ell-1)/((ell^2-1)*(ell^2-ell)) : ell in PrimeDivisors(n) | (k-1) mod ell eq 0],1);
		A[k] := N1*N2/EulerPhi(n);
	end for;
	return A;
end function;

// Given a non-CM elliptic curve E and positive 
// integer n, output an associative array with the 
// values of delta^cyc_{E,n,k}(L),  as in Eq. 48, 
// for all units k modulo n. The optional parameter 
// B specifies the largest prime used in the 
// infinite product.
DeltaCyclicityAP := function(E,n:B:=10^5)
	G := FindOpenImage(E);
	mE := #BaseRing(G);
	R := Rad(mE);
	L := ComputeL(n,mE);
	pi1 := Red(GL(2,Integers(LCM(mE,L))), GL(2,Integers(mE)));
	Ga := G @@ pi1;
	pi2 := Red(GL(2,Integers(LCM(mE,L))),GL(2,Integers(L)));
	GLft := pi2(Ga);
	C := ConjugacyClasses(GLft);
	S := {**};
	for c in C do
		g := c[3];
            flg := true;
            for ell in PrimeDivisors(R) do
                if Matrix(Integers(ell),g) eq Matrix(Integers(ell),[[1,0],[0,1]]) then
                    flg := false;
                    break;
                end if;
            end for;
		if flg then
			D2 := Integers(GCD(n,L)) ! Determinant(g);
			Include(~S,D2^^c[2]);
		end if;
	end for;
	A := AssociativeArray();
	for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
		A[k] := Multiplicity(S,k)/#GLft;
	end for;
	return A, R;
end function;

// Computes delta^cyc_{E,n,k} via Lemma 5.2 in 
// the simple case that ell divides n yet ell 
// does not divide m_E.
DeltaCyclicityAPSimple := function(n,k,ell)
    alpha := Valuation(n,ell);
    if k mod ell eq 1 then
    	return 1/(ell^alpha - ell^(alpha-1)) * (1 - (ell-1)/((ell^2-1)*(ell^2-ell)));
	else
		return 1/(ell^alpha - ell^(alpha-1));
    end if;
end function;

// Given a non-CM elliptic curve E and positive 
// integer n, use Prop. 4.12 to output an 
// associative array with the values of 
// C^cyc_{E,n,k} for all units k modulo n.
// The optional parameter B specifies the largest 
// prime used in the infinite product.
CyclicityAP := function(E,n:B:=10^5)
    DeltaBiases, R := DeltaCyclicityAP(E,n:B:=B); 
    N2 := &*Include([ 1 - 1/((ell^2-1)*(ell^2-ell)) : ell in PrimesUpTo(B) | (n*R) mod ell ne 0],1);
	A := AssociativeArray();
	for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
		N1 := &*Include([ DeltaCyclicityAPSimple(n,k,ell) : ell in PrimeDivisors(n) | R mod ell ne 0],1);
		A[k] := RealField() ! (DeltaBiases[k]*N1*N2);
        end for;
    return A;
end function;

// Given a Serre curve E and positive 
// integer n, use Thm. 1.7 to output an 
// associative array with the values of 
// C^prime_{E,n,k} for all units k modulo n.
// The optional parameter B specifies the largest 
// prime used in the infinite product.
SerreCurveCyclicityAP := function(E,n:B:=10^5)
    mE := AdelicLevelSerreCurve(E);
    L := ComputeL(n,mE);
    AvgCyc := AvgCyclicityAP(n:B:=B);
    A := AssociativeArray();
    if L mod mE ne 0 then
        for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
            A[k] := RealField() ! AvgCyc[k];
        end for;
    else
        P := &*Include([1/(ell^4 - ell^3 - ell^2 + ell - 1) : ell in PrimeDivisors(L) | (2*n) mod ell ne 0],1);
        for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
            DeltaPrime := SquareFree(Integers() ! Discriminant(IntegralModel(E)));
            tauValue := tauCyc(L,DeltaPrime,n,k);
            A[k] := RealField() ! (AvgCyc[k] * (1 + tauValue * 1/5 * P));
        end for;
    end if;
    return A;
end function;
