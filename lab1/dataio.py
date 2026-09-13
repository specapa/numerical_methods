from pathlib import Path


def _parse_line(line):
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    line = line.replace(",", " ").replace(";", " ")
    return [float(x) for x in line.split()]


def _read_number_rows(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            vals = _parse_line(line)
            if vals is not None:
                rows.append(vals)
    if not rows:
        raise ValueError(f"пустой файл: {path}")
    return rows


def read_matrix(path):
    rows = _read_number_rows(path)
    n = len(rows)
    for i, row in enumerate(rows):
        if len(row) != n:
            raise ValueError(
                f"{path}: строка {i + 1} длины {len(row)}, ожидалось {n}"
            )
    return rows


def read_vector(path):
    rows = _read_number_rows(path)
    if len(rows) == 1:
        return rows[0]
    vec = []
    for row in rows:
        if len(row) != 1:
            raise ValueError(f"{path}: вектор должен быть одним столбцом или одной строкой")
        vec.append(row[0])
    return vec


def read_scalar(path):
    rows = _read_number_rows(path)
    return rows[0][0]


def read_tridiag(path):
    rows = _read_number_rows(path)
    if len(rows) != 4:
        raise ValueError(
            f"{path}: для прогонки нужно 4 строки: поддиагональ, диагональ, наддиагональ, правая часть"
        )
    a, b, c, d = rows
    n = len(b)
    if len(a) != n - 1:
        raise ValueError(f"поддиагональ: {len(a)} чисел, нужно {n - 1}")
    if len(c) != n - 1:
        raise ValueError(f"наддиагональ: {len(c)} чисел, нужно {n - 1}")
    if len(d) != n:
        raise ValueError(f"правая часть: {len(d)} чисел, нужно {n}")
    return a, b, c, d


def data_dir():
    here = Path(__file__).resolve().parent.parent
    return here / "data"
