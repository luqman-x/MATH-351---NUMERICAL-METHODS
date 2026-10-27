import numpy as np

def naive_gaussian_elimination(A, B):
    n = len(B)

    # Forward elimination
    for k in range(n - 1):
        for i in range(k + 1, n):
            factor = A[i][k] / A[k][k]

            for j in range(k, n):
                A[i][j] = A[i][j] - factor * A[k][j]

            B[i] = B[i] - factor * B[k]

    # Back substitution
    X = np.zeros(n)

    X[n - 1] = B[n - 1] / A[n - 1][n - 1]

    for i in range(n - 2, -1, -1):
        sum_ax = 0

        for j in range(i + 1, n):
            sum_ax += A[i][j] * X[j]

        X[i] = (B[i] - sum_ax) / A[i][i]

    return X


# Example system:
# 2x + y - z = 8
# -3x - y + 2z = -11
# -2x + y + 2z = -3

A = np.array([
    [2.0, 1.0, -1.0],
    [-3.0, -1.0, 2.0],
    [-2.0, 1.0, 2.0]
])

B = np.array([8.0, -11.0, -3.0])

solution = naive_gaussian_elimination(A, B)

print("Solution:")
print(solution)