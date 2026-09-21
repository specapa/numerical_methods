from lab4.helpers import add_scaled, copy_vec, linspace, scale, vec_add


def _step_nodes(x0, x_end, h):
    xs, n = linspace(x0, x_end, h)
    return xs, n


def euler(f, x0, y0, h, x_end):
    """Явный метод Эйлера для системы y' = f(x, y).

    f: (x, y_list) -> list той же размерности
    y0: начальный вектор (размерность не фиксирована)
    h > 0; при x_end < x0 интегрирование идёт назад.
    """
    xs, n = _step_nodes(x0, x_end, h)
    step = xs[1] - xs[0] if n > 0 else h
    y = copy_vec(y0)
    ys = [copy_vec(y)]
    for k in range(n):
        dy = f(xs[k], y)
        y = add_scaled(y, dy, step)
        ys.append(y)
    return xs, ys


def rk4_step(f, x, y, h):
    k1 = f(x, y)
    k2 = f(x + 0.5 * h, add_scaled(y, k1, 0.5 * h))
    k3 = f(x + 0.5 * h, add_scaled(y, k2, 0.5 * h))
    k4 = f(x + h, add_scaled(y, k3, h))
    incr = scale(vec_add(k1, scale(k2, 2.0), scale(k3, 2.0), k4), h / 6.0)
    return add_scaled(y, incr, 1.0)


def rk4(f, x0, y0, h, x_end):
    """Классический метод Рунге–Кутты 4-го порядка."""
    xs, n = _step_nodes(x0, x_end, h)
    step = xs[1] - xs[0] if n > 0 else h
    y = copy_vec(y0)
    ys = [copy_vec(y)]
    for k in range(n):
        y = rk4_step(f, xs[k], y, step)
        ys.append(y)
    return xs, ys


def adams4(f, x0, y0, h, x_end, starter=None):
    """Явный метод Адамса 4-го порядка (Адамса–Бэшфорта).

    Стартовые значения в первых четырёх узлах — через starter
    (по умолчанию RK4). Формула (4.25):
    y_{k+1} = y_k + (h/24)*(55 f_k - 59 f_{k-1} + 37 f_{k-2} - 9 f_{k-3})
    """
    if starter is None:
        starter = rk4
    xs, n = _step_nodes(x0, x_end, h)
    step = xs[1] - xs[0] if n > 0 else h
    if n <= 3:
        return starter(f, x0, y0, h, x_end)

    _, ys = starter(f, x0, y0, h, x0 + 3 * step)
    ys = [copy_vec(y) for y in ys]
    fs = [f(xs[i], ys[i]) for i in range(4)]

    y = copy_vec(ys[-1])
    for k in range(3, n):
        combo = vec_add(
            scale(fs[k], 55.0),
            scale(fs[k - 1], -59.0),
            scale(fs[k - 2], 37.0),
            scale(fs[k - 3], -9.0),
        )
        y = add_scaled(y, combo, step / 24.0)
        ys.append(copy_vec(y))
        fs.append(f(xs[k + 1], y))
    return xs, ys


def second_order_to_system(rhs):
    """Свести y'' = rhs(x, y, y') к системе [y, z], z=y'."""

    def f(x, state):
        y, z = state[0], state[1]
        return [z, rhs(x, y, z)]

    return f
