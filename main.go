package main

import (
	"bufio"
	"fmt"
	"math/big"
	"os"
	"strings"
	"time"

	"glv-multiplication/internal/ecc"
	"glv-multiplication/internal/glv"
)

func main() {
	params := ecc.NewBN254CurveParams()
	reader := bufio.NewReader(os.Stdin)

	fmt.Println("BN254 curve (y² = x³ + 3)")

	// Get point coordinates
	fmt.Println("Enter point coordinates for BN254 curve:")
	fmt.Print("X coordinate (hex): ")
	xStr, _ := reader.ReadString('\n')
	fmt.Print("Y coordinate (hex): ")
	yStr, _ := reader.ReadString('\n')

	// Parse coordinates
	x, ok1 := new(big.Int).SetString(strings.TrimSpace(xStr), 16)
	y, ok2 := new(big.Int).SetString(strings.TrimSpace(yStr), 16)
	if !ok1 || !ok2 {
		fmt.Println("Error: Invalid coordinate format")
		return
	}

	P := ecc.NewPoint(x, y)

	// Debug info for curve equation
	ySq := ecc.FieldMul(y, y, params.P)
	xCube := ecc.FieldMul(x, x, params.P)
	xCube = ecc.FieldMul(xCube, x, params.P)
	xCubePlusB := ecc.FieldAdd(xCube, params.B, params.P)

	fmt.Printf("\nCurve equation check (y² = x³ + 3):")
	fmt.Printf("\ny² = %s", ySq.Text(16))
	fmt.Printf("\nx³ + 3 = %s\n", xCubePlusB.Text(16))

	// Verify point is on curve
	if !ecc.IsOnCurve(P, params) {
		fmt.Println("Error: Point is not on the curve!")
		return
	}

	fmt.Println("\nPoint verified to be on curve.")

	// Get scalar
	fmt.Print("\nEnter scalar α (hex): ")
	alphaStr, _ := reader.ReadString('\n')
	alpha, ok := new(big.Int).SetString(strings.TrimSpace(alphaStr), 16)
	if !ok {
		fmt.Println("Error: Invalid scalar format")
		return
	}

	// Classical double-and-add method
	start := time.Now()
	resultClassic := ecc.ScalarMult(alpha, P, params)
	classicTime := time.Since(start)

	// GLV method
	start = time.Now()
	resultGLV := glv.GLVMultiply(alpha, P, params)
	glvTime := time.Since(start)

	// Print results
	fmt.Println("\nResults:")
	fmt.Printf("Input point:  (%s,\n              %s)\n", x.Text(16), y.Text(16))
	fmt.Printf("Scalar α:     %s\n", alpha.Text(16))
	fmt.Printf("\nClassical:    (%s,\n              %s)\n", resultClassic.X.Text(16), resultClassic.Y.Text(16))
	fmt.Printf("GLV:         (%s,\n              %s)\n", resultGLV.X.Text(16), resultGLV.Y.Text(16))
	fmt.Printf("Match:        %v\n", resultClassic.IsEqual(resultGLV))

	fmt.Println("\nPerformance:")
	fmt.Printf("Classical:    %v\n", classicTime)
	fmt.Printf("GLV:         %v\n", glvTime)
	fmt.Printf("Speedup:      %.2fx\n", float64(classicTime)/float64(glvTime))
}
