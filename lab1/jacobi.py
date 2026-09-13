import math

from lab1.matops import copy_mat, identity


def jacobi_eigen(A, eps, max_sweeps=200):
    n = len(A)
    T = copy_mat(A)
    V = identity(n)
    it = 0
    while True:
        p, q = 0, 1
        mx = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                v = abs(T[i][j])
                if v > mx:
                    mx = v
                    p, q = i, j
        if mx < eps:
            break
        it += 1
        if it > max_sweeps * n * n:
            break
        app = T[p][p]
        aqq = T[q][q]
        apq = T[p][q]
        if abs(apq) < 1e-18:
            continue
        tau = (aqq - app) / (2.0 * apq)
        t = math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = c * t
        T[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        T[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        T[p][q] = 0.0
        T[q][p] = 0.0
        for i in range(n):
            if i == p or i == q:
                continue
            aip = T[i][p]
            aiq = T[i][q]
            T[i][p] = c * aip - s * aiq
            T[p][i] = T[i][p]
            T[i][q] = s * aip + c * aiq
            T[q][i] = T[i][q]
        for i in range(n):
            vip = V[i][p]
            viq = V[i][q]
            V[i][p] = c * vip - s * viq
            V[i][q] = s * vip + c * viq
    vals = [T[i][i] for i in range(n)]
    return vals, V, it
