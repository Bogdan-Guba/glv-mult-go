package main

import (
	"fmt"
	"math/big"

	"glv-multiplication/internal/ecc"
	"glv-multiplication/internal/glv"
)

func main() {

	// 1. Init params elliptic curve BN254
	params := ecc.NewBN254CurveParams()
	fmt.Printf("BN254 initialized:\n")
	fmt.Printf("  P : %s\n", params.P.Text(16))
	fmt.Printf("  Q : %s\n", params.Q.Text(16))
	fmt.Printf("  B : %s\n\n", params.B.String())

	// 2. Generate G1
	g1 := ecc.GetG1BN254()
	fmt.Printf(" G1 BN254:\n")
	fmt.Printf("  X: %s\n", g1.X.Text(16))
	fmt.Printf("  Y: %s\n\n", g1.Y.Text(16))

	// 3. Test is G1 o curve
	fmt.Println("Test is G1 on curve")
	ySq := ecc.FieldMul(g1.Y, g1.Y, params.P)               // Y^2
	xCubed := ecc.FieldMul(g1.X, g1.X, params.P)            // X^2
	xCubed = ecc.FieldMul(xCubed, g1.X, params.P)           // X^3
	xCubedPlusB := ecc.FieldAdd(xCubed, params.B, params.P) // X^3 + B

	if ySq.Cmp(xCubedPlusB) == 0 {
		fmt.Println("  [+] G1 is on the curve")
	} else {
		fmt.Println("  [-] G1 is NOT on the curve!")
		fmt.Printf("    Y^2 mod P: %s\n", ySq.Text(16))
		fmt.Printf("    X^3 + B mod P: %s\n", xCubedPlusB.Text(16))
	}
	fmt.Println()

	// 4. Test doubling point
	fmt.Println("Test doubling G1")
	g1Double := ecc.DoublePoint(g1, params)
	if g1Double.IsInfinity {
		fmt.Println("  [-] [2]G1 is on infinity point")
	} else {
		fmt.Printf("  [2]G1 X: %s\n", g1Double.X.Text(16))
		fmt.Printf("  [2]G1 Y: %s\n", g1Double.Y.Text(16))
	}
	fmt.Println()

	// 5. Test adding point (G1 + G1 = [2]G1)
	fmt.Println("Test adding point G1 + G1...")
	g1PlusG1 := ecc.AddPoints(g1, g1, params)
	if g1PlusG1.IsEqual(g1Double) {
		fmt.Println("  [+] G1 + G1 = [2]G1 ")
	} else {
		fmt.Println("  [-] G1 + G1 != [2]G1.")
		fmt.Printf("    (G1+G1) X: %s\n", g1PlusG1.X.Text(16))
		fmt.Printf("    (G1+G1) Y: %s\n", g1PlusG1.Y.Text(16))
	}
	fmt.Println()

	// 6. Тest ScalarMult
	fmt.Println("Test ScalarMult ([3]G1)...")
	alphaSmall := big.NewInt(3)
	g1Triple := ecc.ScalarMult(alphaSmall, g1, params)
	expectedG1Triple := ecc.AddPoints(g1Double, g1, params) // [2]G1 + G1 = [3]G1

	if g1Triple.IsEqual(expectedG1Triple) {
		fmt.Println("  [+] ScalarMult ([3]G1) is work correct")
		fmt.Printf("    [3]G1 X: %s\n", g1Triple.X.Text(16))
		fmt.Printf("    [3]G1 Y: %s\n", g1Triple.Y.Text(16))
	} else {
		fmt.Println("  [-] ScalarMult ([3]G1) isn`t correct")
		fmt.Printf("      exprected [3]G1 X: %s\n", expectedG1Triple.X.Text(16))
		fmt.Printf("      expected [3]G1 Y: %s\n", expectedG1Triple.Y.Text(16))
		fmt.Printf("      resultant [3]G1 X: %s\n", g1Triple.X.Text(16))
		fmt.Printf("      resultant [3]G1 Y: %s\n", g1Triple.Y.Text(16))
	}
	fmt.Println()

	// test mult by 0
	fmt.Println("test ScalarMult ([0]G1)...")
	zero := big.NewInt(0)
	g1Zero := ecc.ScalarMult(zero, g1, params)
	if g1Zero.IsInfinity {
		fmt.Println("  [+] ScalarMult ([0]G1) is point at infinity.")
	} else {
		fmt.Println("  [-] ScalarMult ([0]G1) is NOT point at infinity!")
	}
	fmt.Println()

	omega := glv.FindOmega(params.P)
	fmt.Printf("(nontrivial cube root of 1 mod p): %s\n", omega.Text(16))

	// Test phi(G1) = (omega * x, y)
	phiG1 := ecc.NewPoint(
		ecc.FieldMul(g1.X, omega, params.P),
		g1.Y,
	)

	// Test is phi(G1) on the curve
	phiYSq := ecc.FieldMul(phiG1.Y, phiG1.Y, params.P)
	phiXCubed := ecc.FieldMul(phiG1.X, phiG1.X, params.P)
	phiXCubed = ecc.FieldMul(phiXCubed, phiG1.X, params.P)
	phiXCubedPlusB := ecc.FieldAdd(phiXCubed, params.B, params.P)

	if phiYSq.Cmp(phiXCubedPlusB) == 0 {
		fmt.Println("  [+] phi(G1) is on the curve.")
	} else {
		fmt.Println("  [-] phi(G1) is NOT on the curve!")
	}

	//test find lamda
	lambda := glv.FindLambdaBN254(params.Q)
	fmt.Printf("  [+] lambda (BN254): %s\n", lambda.Text(16))

	l2 := new(big.Int).Mul(lambda, lambda)
	l2.Add(l2, lambda)
	l2.Add(l2, big.NewInt(1))
	l2.Mod(l2, params.Q)
	fmt.Printf("  [*] lambda^2 + lambda + 1 mod q = %s\n", l2.Text(16))

}
