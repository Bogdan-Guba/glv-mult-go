package ecc

import (
	"fmt"
	"math/big"
)

// (a + b) mod p
func FieldAdd(a, b, p *big.Int) *big.Int {
	res := new(big.Int).Add(a, b)
	return res.Mod(res, p)
}

// (a - b) mod p
func FieldSub(a, b, p *big.Int) *big.Int {
	res := new(big.Int).Sub(a, b)
	res.Mod(res, p)
	if res.Sign() == -1 { // if data is negative
		res.Add(res, p)
	}
	return res
}

// (a * b) mod p
func FieldMul(a, b, p *big.Int) *big.Int {
	res := new(big.Int).Mul(a, b)
	return res.Mod(res, p)
}

// (-a) mod p
func FieldNeg(a, p *big.Int) *big.Int {
	res := new(big.Int).Neg(a)
	res.Mod(res, p)
	if res.Sign() == -1 { // if data is negative
		res.Add(res, p)
	}
	return res
}

// (a^exp) mod p.
func FieldExp(a, exp, p *big.Int) *big.Int {
	return new(big.Int).Exp(a, exp, p)
}

// (a^(-1)) mod p
func FieldInv(a, p *big.Int) *big.Int {
	res := new(big.Int).ModInverse(a, p)
	if res == nil {

		panic(fmt.Sprintf("Modular inverse does not exist for %s mod %s (a = 0?)", a.String(), p.String()))
	}
	return res
}
