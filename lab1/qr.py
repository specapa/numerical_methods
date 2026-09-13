import math

from lab1.matops import copy_mat, mat_mul, zeros


def qr_decomp(A):
    n = len(A)
    Q = zeros(n)
    R = zeros(n)
    cols = [[A[i][j] for i in range(n)] for j in range(n)]
    for j in range(n):
        v = cols[j][:]
        for i in range(j):
            qi = [Q[r][i] for r in range(n)]
            rij = 0.0
            for r in range(n):
                rij += qi[r] * v[r]
            R[i][j] = rij
            for r in range(n):
                v[r] -= rij * qi[r]
        nrm = 0.0
        for r in range(n):
            nrm += v[r] * v[r]
        nrm = math.sqrt(nrm)
        if nrm == 0:
            raise ValueError("нулевой столбец при QR")
        R[j][j] = nrm
        for r in range(n):
            Q[r][j] = v[r] / nrm
    return Q, R


def _subdiag_small(T, i, eps):
    return abs(T[i][i - 1]) <= eps * (abs(T[i - 1][i - 1]) + abs(T[i][i]))


def _qr_converged(T, eps):
    n = len(T)
    i = n - 1
    while i > 0:
        if _subdiag_small(T, i, eps):
            i -= 1
            continue
        if i >= 2 and _subdiag_small(T, i - 1, eps):
            i -= 2
            continue
        if i == 1:
            return True
        return False
    return True


def _eigs_from_schur(T, eps):
    n = len(T)
    eigs = []
    i = 0
    while i < n:
        if i == n - 1 or _subdiag_small(T, i + 1, eps):
            eigs.append(T[i][i])
            i += 1
            continue
        a = T[i][i]
        b = T[i][i + 1]
        c = T[i + 1][i]
        d = T[i + 1][i + 1]
        tr = a + d
        det = a * d - b * c
        disc = tr * tr - 4.0 * det
        if disc >= 0:
            sd = math.sqrt(disc)
            eigs.append(0.5 * (tr + sd))
            eigs.append(0.5 * (tr - sd))
        else:
            sd = math.sqrt(-disc)
            re = tr / 2.0
            im = sd / 2.0
            eigs.append(complex(re, im))
            eigs.append(complex(re, -im))
        i += 2
    return eigs


def qr_eigenvalues(A, eps, max_iter=5000):
    T = copy_mat(A)
    n = len(T)
    it = 0
    for it in range(1, max_iter + 1):
        Q, R = qr_decomp(T)
        T = mat_mul(R, Q)
        if _qr_converged(T, eps):
            break
    return _eigs_from_schur(T, eps), T, it
