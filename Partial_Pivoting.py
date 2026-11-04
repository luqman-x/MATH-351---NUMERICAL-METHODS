import numpy as np

def gaussian_elimination_partial_pivoting(A, B):
    # Convert inputs to floating-point arrays
    A = np.array(A, dtype=float).copy()
    B = np.array(B, dtype=float).copy()

    n = len(B)

    if A.shape != (n, n):
        raise ValueError("A must be a square matrix matching B.")

    def display_matrix():
        """Display the augmented matrix [A | B]."""
        for i in range(n):
            row = "  ".join(f"{A[i, j]:9.4f}" for j in range(n))
            print(f"[ {row} | {B[i]:9.4f} ]")

  
    print(" GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING")
    

    print("\nORIGINAL AUGMENTED MATRIX [A | B]")
    display_matrix()

    
    # STEP 2: FORWARD ELIMINATION WITH PARTIAL PIVOTING
    print("\nFORWARD ELIMINATION WITH PARTIAL PIVOTING")

    for k in range(n - 1):

        print(f"\n{'-' * 55}")
        print(f"Pivot column {k + 1}")

        # Find row with largest absolute entry in column k
        pivot_row = k + np.argmax(np.abs(A[k:n, k]))

        print(f"Current pivot: A[{k+1},{k+1}] = {A[k, k]:.4f}")
        print(f"Largest candidate is in row {pivot_row + 1}")
        print(f"Selected pivot = {A[pivot_row, k]:.4f}")

        # Check whether the column has a usable pivot
        if abs(A[pivot_row, k]) < 1e-12:
            raise ValueError(
                "No usable pivot found. The matrix may be singular."
            )

        # Swap rows if necessary
        if pivot_row != k:
            print(f"\nSwapping row {k+1} with row {pivot_row+1}")

            A[[k, pivot_row]] = A[[pivot_row, k]]
            B[[k, pivot_row]] = B[[pivot_row, k]]

            print("Augmented matrix after row swap:")
            display_matrix()
        else:
            print("No row swap needed.")

        # Eliminate entries below the pivot
        for i in range(k + 1, n):

            factor = A[i, k] / A[k, k]

            print(f"\nEliminating A[{i+1},{k+1}] = {A[i, k]:.4f}")
            print(
                f"Multiplier = {A[i, k]:.4f} / "
                f"{A[k, k]:.4f} = {factor:.4f}"
            )
            print(
                f"Row operation: R{i+1} = "
                f"R{i+1} - ({factor:.4f})R{k+1}"
            )

            A[i, k:] = A[i, k:] - factor * A[k, k:]
            B[i] = B[i] - factor * B[k]

            # Remove tiny floating-point residuals
            A[i, k] = 0.0

            print("\nUpdated augmented matrix:")
            display_matrix()

    # Check the final pivot
    if abs(A[n - 1, n - 1]) < 1e-12:
        raise ValueError(
            "The final pivot is zero or near zero. "
            "The system may not have a unique solution."
        )

    print("\nForward elimination completed.")
    print("\nUpper triangular augmented matrix:")
    display_matrix()

   
    # STEP 3: BACK SUBSTITUTION

    print("\nBACK SUBSTITUTION")

    X = np.zeros(n)

    for i in range(n - 1, -1, -1):

        print(f"\nCalculating variable {i+1}:")

        known_sum = 0.0

        for j in range(i + 1, n):
            term = A[i, j] * X[j]
            known_sum += term

            print(
                f"A[{i+1},{j+1}] * x{j+1} "
                f"= {A[i, j]:.4f} * {X[j]:.4f} "
                f"= {term:.4f}"
            )

        print(f"Right-hand side = {B[i]:.4f}")
        print(f"Sum of known terms = {known_sum:.4f}")

        print(
            f"x{i+1} = ({B[i]:.4f} - "
            f"{known_sum:.4f}) / {A[i, i]:.4f}"
        )

        X[i] = (B[i] - known_sum) / A[i, i]

        print(f"x{i+1} = {X[i]:.6f}")

   
    # STEP 4: FINAL SOLUTION
    print("\n FINAL SOLUTION")
    print("-" * 10)

    for i in range(n):
        print(f"x{i+1} = {X[i]:.6f}")

    print("COMPUTATION COMPLETED SUCCESSFULLY")

    return X


# Example:
# y + 2z = 5
# 2x + 3y + z = 10
# 4x + y + 3z = 13

A = np.array([
    [0, 1, 2],
    [2, 3, 1],
    [4, 1, 3]
])

B = np.array([5, 10, 13])

solution = gaussian_elimination_partial_pivoting(A, B)