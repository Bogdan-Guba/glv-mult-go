package ecc

import (
	"math/big"
)

type Point struct {
	X          *big.Int
	Y          *big.Int
	IsInfinity bool // true, if it neutrall element
}

func NewPoint(x, y *big.Int) *Point {
	if x == nil || y == nil {
		return &Point{IsInfinity: true}
	}
	return &Point{X: x, Y: y, IsInfinity: false}
}

func (p *Point) IsEqual(other *Point) bool {
	if p.IsInfinity && other.IsInfinity {
		return true
	}
	if p.IsInfinity != other.IsInfinity {
		return false
	}
	return p.X.Cmp(other.X) == 0 && p.Y.Cmp(other.Y) == 0
}

func (p *Point) Negation(params *CurveParams) *Point {
	if p.IsInfinity {
		return NewInfinity()
	}
	return NewPoint(p.X, new(big.Int).Neg(p.Y).Mod(new(big.Int).Neg(p.Y), params.P))
}

// AddPoints
func AddPoints(P1, P2 *Point, params *CurveParams) *Point {

	if P1.IsInfinity {
		return P2
	}

	if P2.IsInfinity {
		return P1
	}

	if P1.IsEqual(P2) {
		return DoublePoint(P1, params)
	}

	if P1.X.Cmp(P2.X) == 0 && FieldAdd(P1.Y, P2.Y, params.P).Cmp(big.NewInt(0)) == 0 {
		return NewPoint(nil, nil) // point at infinity
	}

	// Add two  not indentical points
	deltaY := FieldSub(P2.Y, P1.Y, params.P)
	deltaX := FieldSub(P2.X, P1.X, params.P)
	invDeltaX := FieldInv(deltaX, params.P)

	slope := FieldMul(deltaY, invDeltaX, params.P)

	// X3 = slope^2 - X1 - X2 mod P
	x3 := FieldMul(slope, slope, params.P)
	x3 = FieldSub(x3, P1.X, params.P)
	x3 = FieldSub(x3, P2.X, params.P)

	// Y3 = slope * (X1 - X3) - Y1 mod P
	y3 := FieldSub(P1.X, x3, params.P)
	y3 = FieldMul(y3, slope, params.P)
	y3 = FieldSub(y3, P1.Y, params.P)

	return NewPoint(x3, y3)
}

// dounle point
func DoublePoint(P *Point, params *CurveParams) *Point {
	// if p is infinity
	if P.IsInfinity {
		return NewPoint(nil, nil)
	}

	if P.Y.Cmp(big.NewInt(0)) == 0 {
		return NewPoint(nil, nil)
	}

	// Duble point
	xSq := FieldMul(P.X, P.X, params.P)                 // X^2
	numerator := FieldMul(big.NewInt(3), xSq, params.P) // 3 * X^2
	numerator = FieldAdd(numerator, params.B, params.P) // 3 * X^2 + B

	denominator := FieldMul(big.NewInt(2), P.Y, params.P) // 2 * Y
	invDenominator := FieldInv(denominator, params.P)     // (2 * Y)^(-1)

	slope := FieldMul(numerator, invDenominator, params.P)

	// X3 = slope^2 - 2 * X mod P
	x3 := FieldMul(slope, slope, params.P)
	x3 = FieldSub(x3, FieldMul(big.NewInt(2), P.X, params.P), params.P) // 2 * X

	// Y3 = slope * (X - X3) - Y mod P
	y3 := FieldSub(P.X, x3, params.P)
	y3 = FieldMul(y3, slope, params.P)
	y3 = FieldSub(y3, P.Y, params.P)

	return NewPoint(x3, y3)
}

// ScalarMult by metod of double-and-add.
func ScalarMult(k *big.Int, P *Point, params *CurveParams) *Point {
	if k.Sign() < 0 {
		kAbs := new(big.Int).Neg(k)
		res := ScalarMult(kAbs, P, params)
		return res.Negation(params)
	}
	result := NewPoint(nil, nil) // init inf point
	current := P

	for i := 0; i < k.BitLen(); i++ {
		if k.Bit(i) == 1 {
			result = AddPoints(result, current, params)
		}
		current = DoublePoint(current, params)
	}
	return result
}

// create inf point
func NewInfinity() *Point {
	return NewPoint(nil, nil)
}

func IsOnCurve(P *Point, params *CurveParams) bool {
	if P.IsInfinity {
		return true
	}
	// y² = x³ + b
	ySq := FieldMul(P.Y, P.Y, params.P)
	xCube := FieldMul(P.X, P.X, params.P)
	xCube = FieldMul(xCube, P.X, params.P)
	xCubePlusB := FieldAdd(xCube, params.B, params.P)
	return ySq.Cmp(xCubePlusB) == 0
}
