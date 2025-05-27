package glv

import (
	"math/big"
)

func FindOmega(p *big.Int) *big.Int {
	// Euler's theorem: ω = g^((p-1)/3) is a cube root of unity
	one := big.NewInt(1)
	exp := new(big.Int).Sub(p, one)
	exp.Div(exp, big.NewInt(3))

	// Try multiple generators (e.g., g = 2, 3, 5...) to find ω ≠ 1
	candidates := []*big.Int{big.NewInt(2), big.NewInt(3), big.NewInt(5)}
	for _, g := range candidates {
		omega := new(big.Int).Exp(g, exp, p)
		if omega.Cmp(one) != 0 {
			// Verify omega^3 ≡ 1
			cube := new(big.Int).Exp(omega, big.NewInt(3), p)
			if cube.Cmp(one) == 0 {
				return omega
			}
		}
	}
	panic("No valid ω found (cube root of 1 ≠ 1)")
}
