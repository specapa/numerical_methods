from pathlib import Path


def _parse_line(line):
    line = line.strip().lstrip("\ufeff")
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


def read_vector(path):
    rows = _read_number_rows(path)
    if len(rows) == 1:
        return rows[0]
    vec = []
    for row in rows:
        if len(row) != 1:
            raise ValueError(f"{path}: вектор должен быть одной строкой или столбцом")
        vec.append(row[0])
    return vec


def read_scalar(path):
    rows = _read_number_rows(path)
    return rows[0][0]


def read_interval(path):
    v = read_vector(path)
    if len(v) != 2:
        raise ValueError(f"{path}: ожидалось два числа a b")
    return v[0], v[1]


def read_bc(path):
    """Граничное условие α*y + β*y' = γ — три числа."""
    v = read_vector(path)
    if len(v) != 3:
        raise ValueError(f"{path}: ожидалось α β γ")
    return v[0], v[1], v[2]


def read_coeff_table(path):
    """Таблица коэффициентов линейного ОДУ: x  a  b  c  d  (a y''+b y'+c y = d).

    Можно задать только столбцы без x — тогда узлы берутся из сетки отдельно.
    """
    rows = _read_number_rows(path)
    return rows


def data_dir():
    return Path(__file__).resolve().parent.parent / "data"
