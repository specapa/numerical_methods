from lab1.matops import copy_mat, identity, mat_vec, zeros


def lu_decompose(A):
    n = len(A)
    U = copy_mat(A)
    L = zeros(n)
    P = identity(n)
    swaps = 0
    for k in range(n):
        piv = k
        mx = abs(U[k][k])
        for i in range(k + 1, n):
            v = abs(U[i][k])
            if v > mx:
                mx = v
                piv = i
        if mx == 0:
            raise ValueError("матрица вырождена, LU невозможно")
        if piv != k:
            U[k], U[piv] = U[piv], U[k]
            P[k], P[piv] = P[piv], P[k]
            L[k], L[piv] = L[piv], L[k]
            swaps += 1
        L[k][k] = 1.0
        akk = U[k][k]
        for i in range(k + 1, n):
            mult = U[i][k] / akk
            L[i][k] = mult
            U[i][k] = 0.0
            for j in range(k + 1, n):
                U[i][j] -= mult * U[k][j]
    return L, U, P, swaps


def forward_sub(L, b):
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        s = b[i]
        Li = L[i]
        for j in range(i):
            s -= Li[j] * y[j]
        y[i] = s / Li[i]
    return y


def back_sub(U, y):
    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        Ui = U[i]
        for j in range(i + 1, n):
            s -= Ui[j] * x[j]
        x[i] = s / Ui[i]
    return x


def lu_solve(L, U, P, b):
    pb = mat_vec(P, b)
    y = forward_sub(L, pb)
    return back_sub(U, y)


def lu_inverse(L, U, P):
    n = len(L)
    inv = zeros(n)
    e = [0.0] * n
    for j in range(n):
        e[j] = 1.0
        col = lu_solve(L, U, P, e)
        e[j] = 0.0
        for i in range(n):
            inv[i][j] = col[i]
    return inv


def lu_det(U, swaps):
    d = 1.0
    for i in range(len(U)):
        d *= U[i][i]
    if swaps % 2:
        d = -d
    return d
