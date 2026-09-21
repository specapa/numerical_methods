from lab4.helpers import copy_vec, gauss_solve, linspace
from lab4.ode_ivp import rk4, second_order_to_system


def secant_root(phi, eta0, eta1, tol=1e-10, max_iter=50):
    f0 = phi(eta0)
    f1 = phi(eta1)
    a, b = eta0, eta1
    fa, fb = f0, f1
    for _ in range(max_iter):
        if abs(fb - fa) < 1e-30:
            raise ValueError("метод секущих: нулевая разность φ")
        a_next = b - fb * (b - a) / (fb - fa)
        fa, a = fb, b
        b = a_next
        fb = phi(b)
        if abs(fb) < tol or abs(b - a) < tol:
            return b, fb
    return b, fb


def _ic_from_eta(alpha, beta, gamma, eta):
    if abs(beta) < 1e-15 and abs(alpha) > 1e-15:
        return [gamma / alpha, eta]
    if abs(alpha) < 1e-15 and abs(beta) > 1e-15:
        return [eta, gamma / beta]
    if abs(beta) > 1e-15:
        return [eta, (gamma - alpha * eta) / beta]
    raise ValueError("вырожденное граничное условие")


def shooting(
    rhs,
    a,
    b,
    h,
    bc_left,
    bc_right,
    eta0=0.0,
    eta1=1.0,
    tol=1e-10,
    ivp_solver=None,
    from_right=False,
    singular_left=False,
):
    """Метод стрельбы для краевой задачи y'' = rhs(x,y,y').

    bc_left / bc_right: (α, β, γ) для α y + β y' = γ.
    from_right=True — пристрелка с правого конца.
    singular_left=True — не интегрировать до x=a (вырождение a(x));
    условие слева проверяется экстраполяцией с узла a+h.
    """
    if ivp_solver is None:
        ivp_solver = rk4
    al, bl, gl = bc_left
    ar, br, gr = bc_right
    fsys = second_order_to_system(rhs)

    if not from_right:

        def phi(eta):
            y0 = _ic_from_eta(al, bl, gl, eta)
            xs, ys = ivp_solver(fsys, a, y0, h, b)
            yb, ypb = ys[-1][0], ys[-1][1]
            return ar * yb + br * ypb - gr

        eta, resid = secant_root(phi, eta0, eta1, tol=tol)
        xs, ys = ivp_solver(fsys, a, _ic_from_eta(al, bl, gl, eta), h, b)
        return xs, ys, eta, resid

    x_stop = a + h if singular_left else a

    def _left_residual(xs, ys):
        ya, ypa = ys[-1][0], ys[-1][1]
        if singular_left:
            x1 = xs[-1]
            ypp = rhs(x1, ya, ypa)
            ya_ext = ya - h * ypa + 0.5 * h * h * ypp
            ypa_ext = ypa - h * ypp
            return al * ya_ext + bl * ypa_ext - gl
        return al * ya + bl * ypa - gl

    def phi(eta):
        y0 = _ic_from_eta(ar, br, gr, eta)
        xs, ys = ivp_solver(fsys, b, y0, h, x_stop)
        return _left_residual(xs, ys)

    eta, resid = secant_root(phi, eta0, eta1, tol=tol)
    xs, ys = ivp_solver(fsys, b, _ic_from_eta(ar, br, gr, eta), h, x_stop)
    if singular_left:
        if abs(bl) < 1e-15 and abs(al) > 1e-15:
            ya = gl / al
            ypa = ys[-1][1] - h * rhs(xs[-1], ys[-1][0], ys[-1][1])
        else:
            ya = ys[-1][0] - h * ys[-1][1] + 0.5 * h * h * rhs(
                xs[-1], ys[-1][0], ys[-1][1]
            )
            ypa = ys[-1][1]
        xs = xs + [a]
        ys = ys + [[ya, ypa]]
    xs = list(reversed(xs))
    ys = list(reversed(ys))
    return xs, ys, eta, resid


def finite_difference(
    a_coef,
    b_coef,
    c_coef,
    d_coef,
    a,
    b,
    h,
    bc_left,
    bc_right,
    left_order=2,
    right_order=2,
):
    """Конечно-разностный метод для a(x)y'' + b(x)y' + c(x)y = d(x).

    a_coef,...,d_coef — callable(x) -> float.
    Граничные условия α y + β y' = γ; производные — односторонние O(h) или O(h^2).
    Размерность системы = число узлов (не зашита в код).
    """
    xs, n = linspace(a, b, h)
    N = n
    m = N + 1
    A = [[0.0] * m for _ in range(m)]
    rhs = [0.0] * m
    al, bl, gl = bc_left
    ar, br, gr = bc_right

    def fill_left():
        if abs(bl) < 1e-15:
            A[0][0] = al
            rhs[0] = gl
            return
        if left_order >= 2:
            A[0][0] = al + bl * (-3.0) / (2.0 * h)
            A[0][1] = bl * 4.0 / (2.0 * h)
            A[0][2] = bl * (-1.0) / (2.0 * h)
        else:
            A[0][0] = al + bl * (-1.0) / h
            A[0][1] = bl / h
        rhs[0] = gl

    def fill_right():
        if abs(br) < 1e-15:
            A[N][N] = ar
            rhs[N] = gr
            return
        if right_order >= 2:
            A[N][N - 2] = br * 1.0 / (2.0 * h)
            A[N][N - 1] = br * (-4.0) / (2.0 * h)
            A[N][N] = ar + br * 3.0 / (2.0 * h)
        else:
            A[N][N - 1] = br * (-1.0) / h
            A[N][N] = ar + br / h
        rhs[N] = gr

    fill_left()
    fill_right()

    for k in range(1, N):
        xk = xs[k]
        ak = a_coef(xk)
        bk = b_coef(xk)
        ck = c_coef(xk)
        dk = d_coef(xk)
        A[k][k - 1] = ak / (h * h) - bk / (2.0 * h)
        A[k][k] = -2.0 * ak / (h * h) + ck
        A[k][k + 1] = ak / (h * h) + bk / (2.0 * h)
        rhs[k] = dk

    y = gauss_solve(A, rhs)
    return xs, y


def linear_rhs_from_coefs(a_coef, b_coef, c_coef, d_coef):
    """y'' = (d - b y' - c y) / a  для использования в стрельбе."""

    def rhs(x, y, yp):
        ax = a_coef(x)
        if abs(ax) < 1e-14:
            raise ValueError(f"нулевой коэффициент при y'' в x={x}")
        return (d_coef(x) - b_coef(x) * yp - c_coef(x) * y) / ax

    return rhs
