from lab1.dataio import read_matrix, read_scalar, read_vector, data_dir

__all__ = ["read_matrix", "read_scalar", "read_vector", "data_dir", "read_coeffs"]


def read_coeffs(path):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            line = line.replace(",", " ").replace(";", " ")
            rows.extend(float(x) for x in line.split())
    if not rows:
        raise ValueError(f"пустой файл коэффициентов: {path}")
    return rows
