# Cyclicity and Koblitz Constants for APs

This code is for the paper [Opposing Average Congruence Class Biases in the Cyclicity and Koblitz Conjectures for Elliptic Curves](https://arxiv.org/abs/2408.16641) by Sung Min Lee, Jacob Mayle, and Tian Wang.

The main functions are as follows:

* `AvgKoblitz` - Computes $C_{n,k}^{\text{prime}}$ as in Equation (40).
* `AvgCyclicity` - Computes $C_{n,k}^{\text{cyc}}$ as in Equation (21).
* `KoblitzAP` - For a non-CM elliptic curve $E$, computes $C_{E,n,k}^{\text{prime}}$ using Proposition 4.4.
* `CyclicityAP` - For a non-CM elliptic curve $E$, computes $C_{E,n,k}^{\text{cyc}}$ using Proposition 4.12.
* `SerreCurveKoblitzAP` - For a Serre curve $E$, computes $C_{E,n,k}^{\text{prime}}$ using Theorem 1.7.
* `SerreCurveCyclicityAP` - For a Serre curve $E$, computes $C_{E,n,k}^{\text{cyc}}$ using Theorem 1.7.

Note that `SerreCurveKoblitzAP` and `SerreCurveCyclicityAP` require $E$ to be a Serre curve, whereas `KoblitzAP` and `CyclicityAP` work for arbitrary non-CM $E$. The  Serre curve-specific functions use Theorem 1.7 instead of directly working with the adelic image, so are much more efficient for Serre curves with large adelic level compared to the general functions.

The file `Examples.m` includes code for the tables and examples that appear in the paper.

Installation instructions:
1. Install the latest version of Magma (at least V2.27) from [http://magma.maths.usyd.edu.au/magma/](http://magma.maths.usyd.edu.au/magma/).
2. Install D. Zywina's "OpenImage" from [https://github.com/davidzywina/OpenImage](https://github.com/davidzywina/OpenImage).
3. Download and run "CyclicityKoblitzAPs.m" from this repository.
4. For Example 4 only, download and run Sutherland's "galrep.m" from [https://math.mit.edu/~drew/galrep.html](https://math.mit.edu/~drew/galrep.html).

We welcome any questions, comments, or suggestions.
