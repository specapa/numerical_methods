def find_index(xs, x_star, tol=1e-12):
    for i, x in enumerate(xs):
        if abs(x - x_star) <= tol:
            return i
    raise ValueError(f"точка x*={x_star} не найдена среди узлов")


def first_derivative_left(xs, ys, i):
    """Левосторонняя производная (3.18), O(h)."""
    if i < 1:
        raise ValueError("нет левого соседа")
    return (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1])


def first_derivative_right(xs, ys, i):
    """Правосторонняя производная (3.18), O(h)."""
    if i >= len(xs) - 1:
        raise ValueError("нет правого соседа")
    return (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])


def _triple_indices(n, i):
    if n < 3:
        raise ValueError("нужно минимум 3 точки")
    if i <= 0:
        return 0, 1, 2
    if i >= n - 1:
        return n - 3, n - 2, n - 1
    return i - 1, i, i + 1


def first_derivative_quadratic(xs, ys, i):
    """
    Первая производная по квадратичной аппроксимации (3.20), O(h^2).
    На равномерной сетке совпадает с полусуммой левой и правой.
    """
    i0, i1, i2 = _triple_indices(len(xs), i)
    x = xs[i]
    x0, x1, x2 = xs[i0], xs[i1], xs[i2]
    y0, y1, y2 = ys[i0], ys[i1], ys[i2]
    t01 = (y1 - y0) / (x1 - x0)
    t12 = (y2 - y1) / (x2 - x1)
    return (
        t01 * (2.0 * x - x1 - x2) / (x0 - x2)
        + t12 * (2.0 * x - x0 - x1) / (x2 - x0)
    )


def second_derivative_quadratic(xs, ys, i):
    """Вторая производная по квадратичной аппроксимации (3.21)."""
    i0, i1, i2 = _triple_indices(len(xs), i)
    x0, x1, x2 = xs[i0], xs[i1], xs[i2]
    y0, y1, y2 = ys[i0], ys[i1], ys[i2]
    t01 = (y1 - y0) / (x1 - x0)
    t12 = (y2 - y1) / (x2 - x1)
    return 2.0 * (t12 - t01) / (x2 - x0)


def differentiate_at(xs, ys, x_star):
    i = find_index(xs, x_star)
    left = first_derivative_left(xs, ys, i) if i >= 1 else None
    right = first_derivative_right(xs, ys, i) if i < len(xs) - 1 else None
    y1 = first_derivative_quadratic(xs, ys, i)
    y2 = second_derivative_quadratic(xs, ys, i)
    return {
        "index": i,
        "y_prime_left": left,
        "y_prime_right": right,
        "y_prime": y1,
        "y_double_prime": y2,
    }
