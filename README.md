# GLV Scalar Multiplication on BN254 Curve

This project implements the [(GLV)] scalar multiplication algorithm for elliptic curves, demonstrated on the BN254 curve.

## Features

- **Elliptic Curve Arithmetic:** Basic operations (addition, doubling, negation, scalar multiplication) over the BN254 curve:  
   $ y^2 = x^3 + 3 $
  over a large prime field $\mathbb{F}_p$.
- **GLV Endomorphism:** Efficient scalar multiplication using the GLV method, including:
  - Endomorphism search (ω)
  - Eigenvalue search (λ)
  - Scalar decomposition
  - Simultaneous multi-scalar multiplication
- **Classic Scalar Multiplication:** Double-and-add method for comparison.
- **Performance Comparison:** Benchmarks between classic and GLV multiplication.

## Usage

1. **Build and Run:**
    ```bash
    go run main.go
    ```

2. **Input:**
    - Enter the X and Y coordinates of a point on the BN254 curve (in hexadecimal).
    - Enter the scalar α (in hexadecimal).

3. **Output:**
    - The program checks if the point lies on the curve.
    - Computes [α]P using both classic and GLV methods.
    - Prints the results and verifies if they match.
    - Shows timing and speedup.

## Example

```
BN254 curve (y² = x³ + 3):

Enter point coordinates for BN254 curve:
X coordinate (hex): 1
Y coordinate (hex): 2

Curve equation check (y² = x³ + 3):
y² = 4
x³ + 3 = 4

Point verified to be on curve.

Enter scalar α (hex): 123

Results:
Input point:  (1,
              2)
Scalar α:     123

Classical:    (....)
GLV:          (....)
Match:        true

Performance:
Classical:    1.23ms
GLV:          0.98ms
Speedup:      1.25x
```

## Project Structure

- `main.go` — CLI for user input and demonstration.
- `internal/ecc/` — Elliptic curve operations and field arithmetic.
- `internal/glv/` — GLV endomorphism, scalar decomposition, and fast multiplication.

## Notes

- The GLV method is only applicable for certain points (typically the generator or its multiples) due to the mathematical properties of the endomorphism.
- For arbitrary points, classic multiplication is always valid.
- The project is for educational and research purposes.

**Author:**  
Bogdan Guba  
2025