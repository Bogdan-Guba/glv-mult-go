package glv

import (
	"glv-multiplication/internal/ecc"
	"math/big"
)

func FindOmega(p *big.Int) *big.Int {
	// Euler's theorem: w = g^((p-1)/3) is a cube root of unity
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
	panic("No valid w found (cube root of 1 ≠ 1)")
}

func FindLambdaBN254(q *big.Int) *big.Int {
	one := big.NewInt(1)
	// lambda = (-1 + sqrt(-3))/2 mod q
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

func SimultaneousMult(k1, k2 *big.Int, P1, P2 *ecc.Point, w int, params *ecc.CurveParams) *ecc.Point {

	table := make(map[[2]int]*ecc.Point)
	for i := 0; i < (1 << w); i++ {
		for j := 0; j < (1 << w); j++ {
			if i == 0 && j == 0 {
				continue
			}
			pi := ecc.ScalarMult(big.NewInt(int64(i)), P1, params)
			pj := ecc.ScalarMult(big.NewInt(int64(j)), P2, params)
			table[[2]int{i, j}] = ecc.AddPoints(pi, pj, params)
		}
	}

	// 2. k1, k2 crash with block lenght w
	bits := func(k *big.Int, w int) []int {
		res := []int{}
		kk := new(big.Int).Set(k)
		mask := big.NewInt(int64((1 << w) - 1))
		for kk.Sign() > 0 {
			res = append(res, int(new(big.Int).And(kk, mask).Int64()))
			kk.Rsh(kk, uint(w))
		}
		return res
	}
	K := bits(k1, w)
	L := bits(k2, w)
	d := len(K)
	if len(L) > d {
		d = len(L)
	}

	// 3. Main loop
	R := ecc.NewInfinity()
	for i := d - 1; i >= 0; i-- {
		// 4.1 R <- 2^w R
		for j := 0; j < w; j++ {
			R = ecc.DoublePoint(R, params)
		}
		// 4.2 R <- R + (Ki*P1 + Li*P2)
		ki, li := 0, 0
		if i < len(K) {
			ki = K[i]
		}
		if i < len(L) {
			li = L[i]
		}
		if ki != 0 || li != 0 {
			R = ecc.AddPoints(R, table[[2]int{ki, li}], params)
		}
	}
	return R
}
