package glv

import (
	"math/big"
)

func roundDiv(num, den *big.Int) *big.Int {
	quo, rem := new(big.Int).QuoRem(num, den, new(big.Int))
	absRem := new(big.Int).Abs(rem)
	halfDen := new(big.Int).Rsh(den, 1)
	if absRem.Cmp(halfDen) >= 0 {
		if num.Sign() >= 0 {
			quo.Add(quo, big.NewInt(1))
		} else {
			quo.Sub(quo, big.NewInt(1))
		}
	}
	return quo
}
