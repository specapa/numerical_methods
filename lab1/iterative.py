from lab1.matops import vec_norm_inf, vec_sub


def _check_diag(A):
    n = len(A)
    for i in range(n):
        if A[i][i] == 0:
            raise ValueError("нуль на диагонали, итерации невозможны")


def simple_iteration(A, b, eps, x0=None, max_iter=10000):
    n = len(A)
    _check_diag(A)
    x = [0.0] * n if x0 is None else x0[:]
    hist = []
    for k in range(1, max_iter + 1):
        nxt = [0.0] * n
        for i in range(n):
            s = b[i]
            Ai = A[i]
            for j in range(n):
                if j != i:
                    s -= Ai[j] * x[j]
            nxt[i] = s / Ai[i]
        diff = vec_norm_inf(vec_sub(nxt, x))
        hist.append(diff)
        x = nxt
        if diff < eps:
            return x, k, hist
    return x, max_iter, hist


def seidel(A, b, eps, x0=None, max_iter=10000):
    n = len(A)
    _check_diag(A)
    x = [0.0] * n if x0 is None else x0[:]
    hist = []
    for k in range(1, max_iter + 1):
        old = x[:]
        for i in range(n):
            s = b[i]
            Ai = A[i]
            for j in range(n):
                if j != i:
                    s -= Ai[j] * x[j]
            x[i] = s / Ai[i]
        diff = vec_norm_inf(vec_sub(x, old))
        hist.append(diff)
        if diff < eps:
            return x, k, hist
    return x, max_iter, hist
