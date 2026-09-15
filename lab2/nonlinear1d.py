def poly_eval(coeffs, x):
    s = 0.0
    p = 1.0
    for c in coeffs:
        s += c * p
        p *= x
    return s


def poly_deriv_coeffs(coeffs):
    if len(coeffs) <= 1:
        return [0.0]
    return [(i + 1) * coeffs[i + 1] for i in range(len(coeffs) - 1)]


def make_polynomial(coeffs):
    dcoeffs = poly_deriv_coeffs(coeffs)

    def f(x):
        return poly_eval(coeffs, x)

    def df(x):
        return poly_eval(dcoeffs, x)

    return f, df


def make_phi_relaxation(f, lam):
    def phi(x):
        return x - lam * f(x)

    return phi


def simple_iteration(phi, x0, eps, max_iter=10000):
    x = float(x0)
    hist = []
    for k in range(1, max_iter + 1):
        nxt = phi(x)
        err = abs(nxt - x)
        hist.append(err)
        x = nxt
        if err < eps:
            return x, k, hist
    return x, max_iter, hist


def newton(f, df, x0, eps, max_iter=10000):
    x = float(x0)
    hist = []
    for k in range(1, max_iter + 1):
        d = df(x)
        if d == 0.0:
            raise ValueError("нулевая производная, метод Ньютона неприменим")
        nxt = x - f(x) / d
        err = abs(nxt - x)
        hist.append(err)
        x = nxt
        if err < eps:
            return x, k, hist
    return x, max_iter, hist
