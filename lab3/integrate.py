def nodes(a, b, h):
    if h <= 0:
        raise ValueError("шаг должен быть положительным")
    n = int(round((b - a) / h))
    if abs(a + n * h - b) > 1e-10 * max(1.0, abs(b)):
        raise ValueError(f"отрезок [{a}, {b}] не делится на шаг h={h}")
    xs = [a + i * h for i in range(n + 1)]
    xs[-1] = b
    return xs


def rectangle_midpoint(f, a, b, h):
    """Метод средних прямоугольников."""
    xs = nodes(a, b, h)
    n = len(xs) - 1
    s = 0.0
    for i in range(n):
        mid = 0.5 * (xs[i] + xs[i + 1])
        s += f(mid)
    return s * h


def trapezoid(f, a, b, h):
    xs = nodes(a, b, h)
    n = len(xs) - 1
    ys = [f(x) for x in xs]
    s = 0.5 * (ys[0] + ys[n])
    for i in range(1, n):
        s += ys[i]
    return s * h


def simpson(f, a, b, h):
    xs = nodes(a, b, h)
    n = len(xs) - 1
    if n % 2 != 0:
        raise ValueError("для Симпсона число интервалов должно быть чётным")
    ys = [f(x) for x in xs]
    s = ys[0] + ys[n]
    for i in range(1, n):
        s += (4.0 if i % 2 else 2.0) * ys[i]
    return s * h / 3.0


def runge_romberg(Fh, Fkh, k, p):
    """
    Уточнение / оценка по Рунге–Ромбергу (3.30):
      F ≈ F_h + (F_h - F_{kh}) / (k^p - 1)
    F_h — результат с шагом h, F_kh — с шагом kh (более крупным).
    Возвращает уточнённое значение и оценку погрешности R ≈ (F_h - F_kh)/(k^p - 1).
    """
    denom = k ** p - 1.0
    if denom == 0:
        raise ValueError("некорректные k, p")
    R = (Fh - Fkh) / denom
    F = Fh + R
    return F, R


def integrate_all(f, a, b, h1, h2):
    """
    h1 — более крупный шаг, h2 — более мелкий (обычно h1 = k*h2).
    Порядки: прямоугольники (середины) p=2, трапеции p=2, Симпсон p=4.
    """
    k = h1 / h2
    if abs(k - round(k)) > 1e-10:
        raise ValueError("отношение шагов должно быть целым")
    k = int(round(k))

    methods = [
        ("rectangles", rectangle_midpoint, 2),
        ("trapezoid", trapezoid, 2),
        ("simpson", simpson, 4),
    ]
    out = {}
    for name, method, p in methods:
        F_coarse = method(f, a, b, h1)
        F_fine = method(f, a, b, h2)
        F_rr, R = runge_romberg(F_fine, F_coarse, k, p)
        out[name] = {
            "h1": F_coarse,
            "h2": F_fine,
            "runge_romberg": F_rr,
            "error_est": abs(R),
            "p": p,
            "k": k,
        }
    return out
