
import numpy as np

def naive_gaussian_elimination(A, B):
    n = len(B)

    # Create copies to preserve the original inputs
    A = A.astype(float).copy()
    B = B.astype(float).copy()

    # Display the original system
    print("=" * 65)
    print("          NAIVE GAUSSIAN ELIMINATION")
    print("=" * 65)

    print("\nSTEP 1: ORIGINAL SYSTEM")
    print("-" * 65)

    for i in range(n):
        equation = " + ".join(
            f"({A[i][j]:g})x{j+1}" for j in range(n)
        )
        print(f"{equation} = {B[i]:g}")

    def display_matrix():
        """Display the augmented matrix [A | B]."""
        for i in range(n):
            left = "  ".join(f"{A[i][j]:10.4f}" for j in range(n))
            print(f"[ {left} | {B[i]:10.4f} ]")

    # --------------------------------------------------
    # STEP 2: FORWARD ELIMINATION
    # --------------------------------------------------
    print("\nSTEP 2: FORWARD ELIMINATION")
    print("-" * 65)

    for k in range(n - 1):

        print(f"\nPivot column {k+1}:")
        print(f"Pivot element = A[{k+1},{k+1}] = {A[k][k]:.4f}")

        if abs(A[k][k]) < 1e-12:
            raise ValueError(
                "Zero or near-zero pivot encountered. "
                "Naive Gaussian elimination cannot continue "
                "without pivoting."
            )

        for i in range(k + 1, n):

            factor = A[i][k] / A[k][k]

            print(f"\nEliminating A[{i+1},{k+1}] = {A[i][k]:.4f}")
            print(
                f"Multiplier = A[{i+1},{k+1}] / "
                f"A[{k+1},{k+1}]"
            )
            print(
                f"           = {A[i][k]:.4f} / "
                f"{A[k][k]:.4f} = {factor:.4f}"
            )

            print(
                f"Row operation: R{i+1} = "
                f"R{i+1} - ({factor:.4f})R{k+1}"
            )

            for j in range(k, n):
                A[i][j] = A[i][j] - factor * A[k][j]

            B[i] = B[i] - factor * B[k]

            # Remove tiny floating-point residuals
            for j in range(k, n):
                if abs(A[i][j]) < 1e-12:
                    A[i][j] = 0.0

            print("\nUpdated augmented matrix [A | B]:")
            display_matrix()

    # Check the final pivot
    if abs(A[n-1][n-1]) < 1e-12:
        raise ValueError(
            "The final pivot is zero or near zero. "
            "The system may not have a unique solution."
        )

    print("\nForward elimination completed.")
    print("\nUpper triangular system:")
    display_matrix()

    # --------------------------------------------------
    # STEP 3: BACK SUBSTITUTION
    # --------------------------------------------------
    print("\nSTEP 3: BACK SUBSTITUTION")
    print("-" * 65)

    X = np.zeros(n)

    for i in range(n - 1, -1, -1):

        print(f"\nCalculating x{i+1}:")

        known_terms = 0.0

        for j in range(i + 1, n):
            term = A[i][j] * X[j]
            known_terms += term

            print(
                f"  A[{i+1},{j+1}] * x{j+1} "
                f"= {A[i][j]:.4f} * {X[j]:.4f} "
                f"= {term:.4f}"
            )

        print(f"  Right-hand side = {B[i]:.4f}")
        print(f"  Sum of known terms = {known_terms:.4f}")

        print(
            f"  x{i+1} = ({B[i]:.4f} - "
            f"{known_terms:.4f}) / {A[i][i]:.4f}"
        )

        X[i] = (B[i] - known_terms) / A[i][i]

        print(f"  Therefore, x{i+1} = {X[i]:.6f}")

    # --------------------------------------------------
    # STEP 4: FINAL SOLUTION
    # --------------------------------------------------
    print("\nSTEP 4: FINAL SOLUTION")
    print("-" * 65)

    for i in range(n):
        print(f"x{i+1} = {X[i]:.6f}")

    print("\n" + "=" * 65)
    print("Gaussian elimination completed successfully.")
    print("=" * 65)

    return X


# Example system:
# 2x + y - z = 8
# -3x - y + 2z = -11
# -2x + y + 2z = -3

A = np.array([
    [2.0,  1.0, -1.0],
    [-3.0, -1.0,  2.0],
    [-2.0,  1.0,  2.0]
])

B = np.array([8.0, -11.0, -3.0])

solution = naive_gaussian_elimination(A, B)

