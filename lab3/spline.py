from lab3.linsolve import thomas


def natural_cubic_spline(xs, ys):
    """
    Кубический сплайн с нулевой кривизной на концах: S''(x0)=S''(xn)=0.
    На [x_i, x_{i+1}]:
      S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3
    """
    n = len(xs) - 1
    if n < 1:
        raise ValueError("нужно минимум 2 точки")
    if len(ys) != n + 1:
        raise ValueError("разные длины x и y")

    h = [xs[i + 1] - xs[i] for i in range(n)]
    for hi in h:
        if hi == 0:
            raise ValueError("совпадающие узлы")

    a = list(ys)
    c = [0.0] * (n + 1)

    if n >= 2:
        # система (3.13) относительно c_1..c_{n-1}; c_0 = c_n = 0
        m = n - 1
        sub = [0.0] * (m - 1)
        diag = [0.0] * m
        sup = [0.0] * (m - 1)
        rhs = [0.0] * m
        for i in range(1, n):
            row = i - 1
            diag[row] = 2.0 * (h[i - 1] + h[i])
            rhs[row] = 3.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])
            if row > 0:
                sub[row - 1] = h[i - 1]
            if row < m - 1:
                sup[row] = h[i]
        c_inner = thomas(sub, diag, sup, rhs)
        for i in range(1, n):
            c[i] = c_inner[i - 1]

    b = [0.0] * n
    d = [0.0] * n
    for i in range(n):
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i])
        b[i] = (ys[i + 1] - ys[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0

    return {"a": a, "b": b, "c": c, "d": d, "x": list(xs), "h": h}


def spline_eval(spline, x):
    xs = spline["x"]
    n = len(xs) - 1
    if x <= xs[0]:
        i = 0
    elif x >= xs[n]:
        i = n - 1
    else:
        i = 0
        for k in range(n):
            if xs[k] <= x <= xs[k + 1]:
                i = k
                break
    dx = x - xs[i]
    a = spline["a"][i]
    b = spline["b"][i]
    c = spline["c"][i]
    d = spline["d"][i]
    return a + b * dx + c * dx * dx + d * dx * dx * dx
