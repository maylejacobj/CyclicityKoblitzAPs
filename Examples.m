/////////////////////////////////////////////////
// Example 1 ////////////////////////////////////
/////////////////////////////////////////////////

// Let E be the curve y^2 = x^3 + 6x - 2
// Note that E does not have complex multiplication
E := EllipticCurve([6,-2]);
HasComplexMultiplication(E);

// Computing koblitz and cyclicity constants (n = 6)
KobAP := KoblitzAP(E,6);
CycAP := CyclicityAP(E,6);
for k in [1,5] do
    k, KobAP[k];
end for;

// Numerical data
NE := Conductor(E);
{* p mod 6 : p in PrimesUpTo(10^6) | NE mod p ne 0 and IsPrime(#ChangeRing(E,GF(p)))*};


/////////////////////////////////////////////////
// Example 2 ////////////////////////////////////
/////////////////////////////////////////////////

// Let E be the curve y^2 = x^3 + 5x - 10
// Note that E does not have complex multiplication
E := EllipticCurve([5,-10]);
HasComplexMultiplication(E);

// Computing koblitz and cyclicity constants (n = 8)
SerreKobAP := SerreCurveKoblitzAP(E,8);
SerreCycAP := SerreCurveCyclicityAP(E,8);
for k in [1,3,5,7] do
    k, SerreKobAP[k], SerreCycAP[k];
end for;

// Checking against general code
KobAP := KoblitzAP(E,8);
CycAP := CyclicityAP(E,8);
for k in [1,3,5,7] do
    k, KobAP[k], CycAP[k];
end for;

// Numerical data
NE := Conductor(E);
{* p mod 8 : p in PrimesUpTo(10^6) | NE mod p ne 0 and IsPrime(#ChangeRing(E,GF(p)))*};
{* p mod 8 : p in PrimesUpTo(10^6) | NE mod p ne 0 and IsCyclic(AbelianGroup(ChangeRing(E,GF(p))))*};

/////////////////////////////////////////////////
// Example 3 ////////////////////////////////////
/////////////////////////////////////////////////

// Let E be the curve y^2 = x^3 - 216x - 1296
// Note that E does not have complex multiplication
E := EllipticCurve([-216,-1296]);
HasComplexMultiplication(E);

// Computing koblitz and cyclicity constants (n = 12)
KobAP := KoblitzAP(E,12);
CycAP := CyclicityAP(E,12);
for k in [1,5,7,11] do
    k, KobAP[k], CycAP[k];
end for;

// Numerical data
NE := Conductor(E);
{* p mod 12 : p in PrimesUpTo(10^6) | NE mod p ne 0 and IsPrime(#ChangeRing(E,GF(p)))*};
{* p mod 12 : p in PrimesUpTo(10^6) | NE mod p ne 0 and IsCyclic(AbelianGroup(ChangeRing(E,GF(p))))*};


/////////////////////////////////////////////////
// Example 4 ////////////////////////////////////
/////////////////////////////////////////////////

// Let E be the curve y^2 = x^3 - 4
// Note that E has complex multiplication
E := EllipticCurve([0,-4]);
HasComplexMultiplication(E);

// Computing C_E^Prime
R<x> := PolynomialRing(Rationals());
K<a> := NumberField(x^2+3);
OK := RingOfIntegers(K);
DK := Integers()!Discriminant(OK);
Cprime := 1/2 * &* [1 - KroneckerSymbol(DK,ell) * (ell^2-ell-1) / ( (ell - KroneckerSymbol(DK,ell)) * (ell-1)^2 ) : ell in PrimesUpTo(10^6) | 6 mod ell ne 0 ];
RealField(10)!Cprime;

// Numerical data
NE := Conductor(E);
{* p mod 6 : p in PrimesUpTo(10^6) | NE mod p ne 0 and IsPrime(#ChangeRing(E,GF(p)))*};

// To compute the mod 6 image of E/K, where K = 
// Q(sqrt(-3)), we employ Sutherland's galrep.
// See #4 in the Readme for instructions.
E := EllipticCurve([0,-4]);
m := 6;
B:=TorsionBasis(WeierstrassModel(E),m);
G,S,phi:=AutomorphismGroup(Parent(B[1][1]));
K := Domain(phi(G!1));
R<x> := PolynomialRing(K);
alpha := Roots(x^2+3)[1][1];
G := sub<G | {g : g in G | phi(g)(alpha) eq alpha}>;
A:=TorsionPoints(B,m);
H := sub<GL(2,Integers(m))|[SigmaMatrix(phi(g),B,A,m): g in Generators(G)]>;
#H; // this gives that |G_E(6)| = 6
I := Matrix(One(H));
// the next line gives the size of Psi_{K,6,k}^{prime} \cap G_E(6)
{* Determinant(h) : h in H | IsUnit(Determinant(h - I))*};

/////////////////////////////////////////////////
// Table 1 //////////////////////////////////////
/////////////////////////////////////////////////

for n in [2..6] do
    A := AvgCyclicityAP(n : B := 10^6);
    for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
        k, RealField(5) ! A[k];
    end for;
end for;

/////////////////////////////////////////////////
// Table 2 //////////////////////////////////////
/////////////////////////////////////////////////

for n in [2..6] do
    A := AvgKoblitzAP(n);
    for k in [i : i in [0..n-1] | GCD(n,i) eq 1] do
        k, RealField(5) ! A[k];
    end for;
end for;
