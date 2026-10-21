from pathlib import Path


def discover(input_dir: Path):
    """Clasifica los archivos de la carpeta. Devuelve (pdfs, dicoms[(dcm, csv|None)], ignorados)."""
    archivos = sorted(p for p in input_dir.iterdir() if p.is_file() and not p.name.startswith("."))
    pdfs = [p for p in archivos if p.suffix.lower() == ".pdf"]
    dcms = [p for p in archivos if p.suffix.lower() == ".dcm"]
    csvs = {p.stem: p for p in archivos if p.suffix.lower() == ".csv"}
    dicoms = [(d, csvs.pop(d.stem, None)) for d in dcms]
    ignorados = [p for p in archivos if p not in pdfs and p not in dcms and p.suffix.lower() != ".csv"]
    ignorados += list(csvs.values())  # CSV sin su DCM
    return pdfs, dicoms, ignorados
