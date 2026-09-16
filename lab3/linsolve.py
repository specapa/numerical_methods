def gauss_solve(A, b):
    """Решение СЛАУ методом Гаусса с частичным выбором главного элемента."""
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for k in range(n):
        piv = k
        mx = abs(M[k][k])
        for i in range(k + 1, n):
            v = abs(M[i][k])
            if v > mx:
                mx = v
                piv = i
        if mx == 0:
            raise ValueError("вырожденная матрица")
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
        akk = M[k][k]
        for i in range(k + 1, n):
            factor = M[i][k] / akk
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]
    return x


def thomas(a, b, c, d):
    """Прогонка: a поддиагональ (n-1), b диагональ (n), c наддиагональ (n-1), d правая часть."""
    n = len(b)
    cp = [0.0] * n
    dp = [0.0] * n
    if b[0] == 0:
        raise ValueError("нулевой элемент на диагонали")
    cp[0] = c[0] / b[0] if n > 1 else 0.0
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        denom = b[i] - a[i - 1] * cp[i - 1]
        if denom == 0:
            raise ValueError("деление на ноль в прогонке")
        if i < n - 1:
            cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / denom
    x = [0.0] * n
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x
