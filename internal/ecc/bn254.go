package ecc

import (
	"math/big"
)

type CurveParams struct {
	P *big.Int
	B *big.Int
	Q *big.Int
}

func NewBN254CurveParams() *CurveParams {

	p, ok := new(big.Int).SetString("30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47", 16)
	if !ok {
		panic("Failed to parse P for BN254 curve parameters.")
	}

	q, ok := new(big.Int).SetString("30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001", 16)
	if !ok {
		panic("Failed to parse Q for BN254 curve parameters.")
	}

	b := big.NewInt(3)

	return &CurveParams{
		P: p,
		B: b,
		Q: q,
	}
}

func GetG1BN254() *Point {

	x := big.NewInt(1)

	y := big.NewInt(2)

	return NewPoint(x, y)
}
