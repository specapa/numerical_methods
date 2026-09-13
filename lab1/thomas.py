def thomas(a, b, c, d):
    n = len(b)
    alpha = [0.0] * n
    beta = [0.0] * n
    if b[0] == 0:
        raise ValueError("нулевой элемент на диагонали")
    alpha[0] = -c[0] / b[0]
    beta[0] = d[0] / b[0]
    for i in range(1, n):
        z = a[i - 1] * alpha[i - 1] + b[i]
        if z == 0:
            raise ValueError("деление на ноль в прогонке")
        if i < n - 1:
            alpha[i] = -c[i] / z
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / z
    x = [0.0] * n
    x[n - 1] = beta[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]
    return x, alpha, beta
