import numpy as np

def lu_decomposition(A, B):

    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)

    n = len(B)

    if A.shape != (n, n):
        raise ValueError("A must be a square matrix.")

    # Preserve original matrices
    A_original = A.copy()
    B_original = B.copy()

    # Initialize matrices
    L = np.eye(n)
    U = A.copy()
    P = np.eye(n)

    print("=" * 65)
    print("       LU DECOMPOSITION WITH PARTIAL PIVOTING")
    print("=" * 65)

    def display(M, name):
        print(f"\n{name} =")
        print(np.array2string(
            M,
            precision=4,
            suppress_small=True
        ))

    display(A, "Original Matrix A")
    display(B.reshape(-1, 1), "Vector B")

    # ------------------------------------------------
    # STEP 1: LU DECOMPOSITION
    # ------------------------------------------------

    print("\nSTEP 1: LU DECOMPOSITION")
    print("-" * 65)

    for k in range(n):

        print(f"\nPIVOT COLUMN {k+1}")

        # Find largest pivot
        pivot_row = k + np.argmax(np.abs(U[k:n, k]))

        print(f"Selected pivot = U[{pivot_row+1},{k+1}]")
        print(f"Pivot value = {U[pivot_row, k]:.4f}")

        if abs(U[pivot_row, k]) < 1e-12:
            raise ValueError("Matrix is singular.")

        # Swap rows if necessary
        if pivot_row != k:

            print(f"\nSwapping row {k+1} with row {pivot_row+1}")

            U[[k, pivot_row]] = U[[pivot_row, k]]
            P[[k, pivot_row]] = P[[pivot_row, k]]

            # Swap previously calculated L entries
            if k > 0:
                L[[k, pivot_row], :k] = L[[pivot_row, k], :k]

            display(U, "Updated U")
            display(L, "Updated L")
            display(P, "Updated P")

        # Elimination
        for i in range(k + 1, n):

            factor = U[i, k] / U[k, k]

            print(f"\nMultiplier m{i+1}{k+1}")
            print(f"= {U[i,k]:.4f} / {U[k,k]:.4f}")
            print(f"= {factor:.4f}")

            L[i, k] = factor

            print(
                f"Row operation: R{i+1} = "
                f"R{i+1} - ({factor:.4f})R{k+1}"
            )

            U[i, k:] = U[i, k:] - factor * U[k, k:]

            U[i, k] = 0.0

            display(L, "Current L")
            display(U, "Current U")

    # ------------------------------------------------
    # STEP 2: VERIFY PA = LU
    # ------------------------------------------------

    print("\nSTEP 2: VERIFY LU DECOMPOSITION")
    print("-" * 65)

    display(P, "Permutation Matrix P")
    display(L, "Lower Triangular Matrix L")
    display(U, "Upper Triangular Matrix U")

    print("\nChecking whether PA = LU:")

    display(P @ A_original, "PA")
    display(L @ U, "LU")

    if np.allclose(P @ A_original, L @ U):
        print("\nVerification successful: PA = LU")
    else:
        raise ArithmeticError("LU decomposition verification failed.")

    # ------------------------------------------------
    # STEP 3: FORWARD SUBSTITUTION
    # ------------------------------------------------

    print("\nSTEP 3: FORWARD SUBSTITUTION")
    print("-" * 65)

    PB = P @ B_original

    display(PB.reshape(-1, 1), "PB")

    Y = np.zeros(n)

    for i in range(n):

        known_sum = 0

        for j in range(i):
            known_sum += L[i, j] * Y[j]

        print(f"\nCalculating Y{i+1}:")

        print(
            f"Y{i+1} = ({PB[i]:.4f} - "
            f"{known_sum:.4f}) / {L[i,i]:.4f}"
        )

        Y[i] = (PB[i] - known_sum) / L[i, i]

        print(f"Y{i+1} = {Y[i]:.6f}")

    display(Y.reshape(-1, 1), "Solution Vector Y")

    # ------------------------------------------------
    # STEP 4: BACKWARD SUBSTITUTION
    # ------------------------------------------------

    print("\nSTEP 4: BACKWARD SUBSTITUTION")
    print("-" * 65)

    X = np.zeros(n)

    for i in range(n - 1, -1, -1):

        known_sum = 0

        for j in range(i + 1, n):
            known_sum += U[i, j] * X[j]

        print(f"\nCalculating X{i+1}:")

        print(
            f"X{i+1} = ({Y[i]:.4f} - "
            f"{known_sum:.4f}) / {U[i,i]:.4f}"
        )

        X[i] = (Y[i] - known_sum) / U[i, i]

        print(f"X{i+1} = {X[i]:.6f}")

    # ------------------------------------------------
    # STEP 5: FINAL SOLUTION
    # ------------------------------------------------

    print("\nSTEP 5: FINAL SOLUTION")
    print("-" * 65)

    for i in range(n):
        print(f"X{i+1} = {X[i]:.6f}")

    print("\n" + "=" * 65)
    print("LU DECOMPOSITION COMPLETED SUCCESSFULLY")
    print("=" * 65)

    return X


# -----------------------------------------------
# EXAMPLE PROBLEM
#
# y + 2z = 5
# 2x + 3y + z = 10
# 4x + y + 3z = 13
# -----------------------------------------------

A = np.array([
    [0, 1, 2],
    [2, 3, 1],
    [4, 1, 3]
])

B = np.array([5, 10, 13])

solution = lu_decomposition(A, B)
