package glv

import (
	"errors"
	"fmt"
	"glv-multiplication/internal/ecc"
	"math/big"
)

// Find omega and alpha for GLV multiplication but this metod is very slow and not for all!!!
// UniversalFindOmegaAndLambda search ω and λ, when:
// 1. w^3 = 1 mod p, ω != 1
// 2. lambda^2 + lambda + 1 = 0 mod q, lambda != 1
// 3. phi(P) = (omega·x, y) == [lambda]P
func FindOmegaAndLambda(P *ecc.Point, params *ecc.CurveParams) (*big.Int, *big.Int, error) {
	p := params.P
	q := params.Q
	one := big.NewInt(1)
	three := big.NewInt(3)
	exp := new(big.Int).Div(new(big.Int).Sub(p, one), three)

	// search w where w^3 = 1 mod p
	for g := int64(2); g < 100; g++ {
		base := big.NewInt(g)
		omega := new(big.Int).Exp(base, exp, p)
		if omega.Cmp(one) == 0 {
			continue
		}
		// check if omega^3 == 1 mod p
		check := new(big.Int).Exp(omega, three, p)
		if check.Cmp(one) != 0 {
			continue
		}
		// brootforce labda, by root of lambda^2 + lambda + 1 = 0 mod q
		for try := int64(1); try < 10000; try++ {
			lambda := big.NewInt(try)
			tmp := new(big.Int).Mul(lambda, lambda)
			tmp.Add(tmp, lambda)
			tmp.Add(tmp, one)
			tmp.Mod(tmp, q)
			if tmp.Sign() != 0 || lambda.Cmp(one) == 0 {
				continue
			}
			// check phi(P) == [lamda]P
			phiP := ecc.NewPoint(
				ecc.FieldMul(P.X, omega, p),
				P.Y,
			)
			res := ecc.ScalarMult(lambda, P, params)
			if res.IsEqual(phiP) {
				return omega, lambda, nil
			}
		}
	}
	return nil, nil, errors.New("failed to find harmony omega and lambda")
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

		seq = append(seq, triple{r2, s2, t2})

		r0, r1 = r1, r2
		s0, s1 = s1, s2
		t0, t1 = t1, t2
	}

	sqrtQ := new(big.Int).Sqrt(q)
	l := 0
	for i := 0; i < len(seq); i++ {
		if seq[i].r.Cmp(sqrtQ) >= 0 {
			l = i
		}
	}

	a1 := new(big.Int).Set(seq[l+1].r)
	b1 := new(big.Int).Neg(seq[l+1].t)

	var a2, b2 *big.Int
	if l+2 >= len(seq) {
		a2 = new(big.Int).Set(seq[l].r)
		b2 = new(big.Int).Neg(seq[l].t)
	} else {
		normL := new(big.Int).Add(new(big.Int).Mul(seq[l].r, seq[l].r), new(big.Int).Mul(seq[l].t, seq[l].t))
		normLp2 := new(big.Int).Add(new(big.Int).Mul(seq[l+2].r, seq[l+2].r), new(big.Int).Mul(seq[l+2].t, seq[l+2].t))
		if normL.Cmp(normLp2) <= 0 {
			a2 = new(big.Int).Set(seq[l].r)
			b2 = new(big.Int).Neg(seq[l].t)
		} else {
			a2 = new(big.Int).Set(seq[l+2].r)
			b2 = new(big.Int).Neg(seq[l+2].t)
		}
	}

	c1 := roundDiv(new(big.Int).Mul(b2, alpha), q)
	c2 := roundDiv(new(big.Int).Neg(new(big.Int).Mul(b1, alpha)), q)

	tmp1 := new(big.Int).Mul(c1, a1)
	tmp2 := new(big.Int).Mul(c2, a2)
	alpha1 := new(big.Int).Sub(alpha, tmp1)
	alpha1.Sub(alpha1, tmp2)

	tmp1 = new(big.Int).Mul(c1, b1)
	tmp2 = new(big.Int).Mul(c2, b2)
	alpha2 := new(big.Int).Add(tmp1, tmp2)
	alpha2.Neg(alpha2)

	// normalization [0, q)
	alpha1.Mod(alpha1, q)
	alpha2.Mod(alpha2, q)

	return alpha1, alpha2

}

func SimultaneousMult(k1, k2 *big.Int, P1, P2 *ecc.Point, w int, params *ecc.CurveParams) *ecc.Point {
	q := params.Q
	k1mod := new(big.Int).Mod(k1, q)
	if k1mod.Sign() < 0 {
		k1mod.Add(k1mod, q)
	}
	k2mod := new(big.Int).Mod(k2, q)
	if k2mod.Sign() < 0 {
		k2mod.Add(k2mod, q)
	}

	if k1mod.Sign() == 0 && k2mod.Sign() == 0 {
		return ecc.NewInfinity()
	}

	twoPowW := 1 << w
	precomputed := make(map[uint64]*ecc.Point)

	for i := 0; i < twoPowW; i++ {
		for j := 0; j < twoPowW; j++ {
			if i == 0 && j == 0 {
				continue
			}
			pi := ecc.ScalarMult(big.NewInt(int64(i)), P1, params)
			qj := ecc.ScalarMult(big.NewInt(int64(j)), P2, params)
			sum := ecc.AddPoints(pi, qj, params)
			key := (uint64(i) << uint(w)) | uint64(j)
			precomputed[key] = sum
		}
	}

	maxBits := k1mod.BitLen()
	if k2mod.BitLen() > maxBits {
		maxBits = k2mod.BitLen()
	}
	d := (maxBits + w - 1) / w
	K := make([]uint64, d)
	L := make([]uint64, d)
	mask := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), uint(w)), big.NewInt(1))

	tmpK := new(big.Int).Set(k1mod)
	tmpL := new(big.Int).Set(k2mod)

	for i := 0; i < d; i++ {
		Ki := new(big.Int).And(tmpK, mask)
		Li := new(big.Int).And(tmpL, mask)
		K[i] = Ki.Uint64()
		L[i] = Li.Uint64()
		tmpK.Rsh(tmpK, uint(w))
		tmpL.Rsh(tmpL, uint(w))
	}

	R := ecc.NewInfinity()
	for i := d - 1; i >= 0; i-- {
		for j := 0; j < w; j++ {
			R = ecc.DoublePoint(R, params)
		}

		Ki := K[i]
		Li := L[i]

		if Ki == 0 && Li == 0 {
			continue
		}
		if Ki >= uint64(twoPowW) || Li >= uint64(twoPowW) {
			panic(fmt.Sprintf("Ki or Li out of bounds: Ki=%d, Li=%d", Ki, Li))
		}
		key := (Ki << uint(w)) | Li
		R = ecc.AddPoints(R, precomputed[key], params)
	}

	// Check if the result needs to be negated
	neg := false
	if k1.Sign() < 0 {
		neg = !neg
	}
	if k2.Sign() < 0 {
		neg = !neg
	}
	if neg {
		R = R.Negation(params)
	}

	return R
}

// GLVMultiply
func GLVMultiply(alpha *big.Int, P *ecc.Point, params *ecc.CurveParams) *ecc.Point {
	w := 2
	// A1. Search omega
	// A2. Search lambda

	omega, lambda, err := FindOmegaAndLambda(P, params)
	if err != nil {
		fmt.Println(err)
		fmt.Println("Set standart values for omega and lambda")
		//hardcoded values for omega and lambda, because my method is very slow
		omega, _ = new(big.Int).SetString("0693b5d6484f0ab41464c9dfe60fd3f4c43d8a7f70c5b863d1c422a0e6e62b1e", 16)
		lambda, _ = new(big.Int).SetString("30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001", 16)
	}
	// Search phi(P)
	phiP := ecc.NewPoint(
		ecc.FieldMul(P.X, omega, params.P),
		P.Y,
	)
	// A3. Decompose alpha = a1 + a2*lambda
	a1, a2 := DecomposeAlpha(alpha, lambda, params.Q)
	// A4. multiple point multiplication
	return SimultaneousMult(a1, a2, P, phiP, w, params)
}
