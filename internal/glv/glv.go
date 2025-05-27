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

func FindLambdaBN254(q *big.Int) *big.Int {
	one := big.NewInt(1)
	// λ = (-1 + sqrt(-3))/2 mod q
	negOne := new(big.Int).Neg(one)
	negOne.Mod(negOne, q)

	three := big.NewInt(3)
	minusThree := new(big.Int).Neg(three)
	minusThree.Mod(minusThree, q)

	// sqrt(-3) mod q
	sqrtMinusThree := new(big.Int).ModSqrt(minusThree, q)
	if sqrtMinusThree == nil {
		panic("sqrt(-3) does not exist mod q")
	}

	num := new(big.Int).Add(negOne, sqrtMinusThree)
	den := big.NewInt(2)
	denInv := new(big.Int).ModInverse(den, q)
	if denInv == nil {
		panic("No inverse for 2 mod q")
	}
	lambda := new(big.Int).Mul(num, denInv)
	lambda.Mod(lambda, q)
	return lambda
}

func DecomposeAlpha(alpha, lambda, q *big.Int) (*big.Int, *big.Int) {

	type triple struct{ r, s, t *big.Int }
	var seq []triple

	r0 := new(big.Int).Set(q)
	r1 := new(big.Int).Set(lambda)
	s0 := big.NewInt(1)
	s1 := big.NewInt(0)
	t0 := big.NewInt(0)
	t1 := big.NewInt(1)

	seq = append(seq, triple{new(big.Int).Set(r0), new(big.Int).Set(s0), new(big.Int).Set(t0)})
	seq = append(seq, triple{new(big.Int).Set(r1), new(big.Int).Set(s1), new(big.Int).Set(t1)})

	for {
		if r1.Sign() == 0 {
			break
		}
		quot := new(big.Int).Div(r0, r1)
		r2 := new(big.Int).Sub(r0, new(big.Int).Mul(quot, r1))
		s2 := new(big.Int).Sub(s0, new(big.Int).Mul(quot, s1))
		t2 := new(big.Int).Sub(t0, new(big.Int).Mul(quot, t1))

		seq = append(seq, triple{new(big.Int).Set(r2), new(big.Int).Set(s2), new(big.Int).Set(t2)})

		r0, r1 = r1, r2
		s0, s1 = s1, s2
		t0, t1 = t1, t2
	}

	// 2. Search l where rl >= sqrt(q)
	sqrtQ := new(big.Int).Sqrt(q)
	l := 0
	for i := 0; i < len(seq); i++ {
		if seq[i].r.Cmp(sqrtQ) >= 0 {
			l = i
		}
	}

	// 3. basis vectors
	a1 := new(big.Int).Set(seq[l+1].r)
	b1 := new(big.Int).Neg(seq[l+1].t)

	var a2, b2 *big.Int
	normL := new(big.Int).Add(new(big.Int).Mul(seq[l].r, seq[l].r), new(big.Int).Mul(seq[l].t, seq[l].t))
	normLp2 := new(big.Int).Add(new(big.Int).Mul(seq[l+2].r, seq[l+2].r), new(big.Int).Mul(seq[l+2].t, seq[l+2].t))
	if normL.Cmp(normLp2) <= 0 {
		a2 = new(big.Int).Set(seq[l].r)
		b2 = new(big.Int).Neg(seq[l].t)
	} else {
		a2 = new(big.Int).Set(seq[l+2].r)
		b2 = new(big.Int).Neg(seq[l+2].t)
	}

	// 4. Round
	c1 := roundDiv(new(big.Int).Mul(b2, alpha), q)
	c2 := roundDiv(new(big.Int).Neg(new(big.Int).Mul(b1, alpha)), q)

	// 5. k1 = alpha – c1*a1 – c2*a2
	tmp1 := new(big.Int).Mul(c1, a1)
	tmp2 := new(big.Int).Mul(c2, a2)
	alpha1 := new(big.Int).Sub(alpha, tmp1)
	alpha1.Sub(alpha1, tmp2)

	// 6. k2 = –(c1*b1 + c2*b2)
	tmp1 = new(big.Int).Mul(c1, b1)
	tmp2 = new(big.Int).Mul(c2, b2)
	alpha2 := new(big.Int).Add(tmp1, tmp2)
	alpha2.Neg(alpha2)

	return alpha1, alpha2
}
