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
            raise ValueError(f"{path}: вектор должен быть одним столбцом или одной строкой")
        vec.append(row[0])
    return vec


def read_scalar(path):
    rows = _read_number_rows(path)
    return rows[0][0]


def read_points(path):
    """Читает таблицу (x, y): либо два столбца, либо две строки."""
    rows = _read_number_rows(path)
    if len(rows) == 2 and len(rows[0]) == len(rows[1]) and len(rows[0]) > 2:
        return rows[0], rows[1]
    xs, ys = [], []
    for row in rows:
        if len(row) != 2:
            raise ValueError(f"{path}: ожидаются пары x y или две строки")
        xs.append(row[0])
        ys.append(row[1])
    return xs, ys


def read_xy_files(x_path, y_path):
    x = read_vector(x_path)
    y = read_vector(y_path)
    if len(x) != len(y):
        raise ValueError(f"разные длины: {x_path} ({len(x)}) и {y_path} ({len(y)})")
    return x, y


def data_dir():
    here = Path(__file__).resolve().parent.parent
    return here / "data" / "lab3"
