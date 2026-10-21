import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pymupdf

from .common import FS_OUT, LEADS, ECGRecord, to_fixed_rate

MM = 72 / 25.4  # puntos PDF por mm
MIN_SEGMENTOS = 100  # una traza de señal tiene miles de segmentos; los marcos y etiquetas, pocos


def _trazas(page):
    ds = [
        d for d in page.get_drawings()
        if d["color"] == (0.0, 0.0, 0.0) and len(d["items"]) > MIN_SEGMENTOS
    ]
    return sorted(ds, key=lambda d: d["rect"].y0 + d["rect"].y1)


def is_vector_page(page):
    return len(_trazas(page)) == len(LEADS)


def _extraer_traza(d, velocidad, ganancia):
    # El primer segmento es la marca vertical de inicio de derivación, no señal.
    pts = [(i[1].x, i[1].y) for i in d["items"][1:]] + [(d["items"][-1][2].x, d["items"][-1][2].y)]
    x, y = np.array(pts).T
    x_u, inv = np.unique(x, return_inverse=True)
    y_u = np.bincount(inv, weights=y) / np.bincount(inv)
    t = (x_u - x_u[0]) / MM / velocidad
    mv = -(y_u - np.median(y_u)) / MM / ganancia  # Y del PDF crece hacia abajo; el cero absoluto se pierde
    return t, mv


def _num(patron, texto, grupos=1, flags=0):
    m = re.search(patron, texto, flags)
    if not m:
        return [np.nan] * grupos
    return [float(g) for g in m.groups()]


def parse_text(texto):
    """Extrae metadatos del texto de la página (formato GE MAC 2000)."""
    lineas = texto.split("\n")
    meta = {}
    m = re.search(r"(\d{2}\.\d{2}\.\d{4}\s+\d{2}:\d{2}:\d{2})", texto)
    if m:
        meta["acquisition_datetime"] = datetime.strptime(m.group(1), "%d.%m.%Y %H:%M:%S").isoformat()
        i = next(k for k, l in enumerate(lineas) if m.group(1) in l)
        if i + 1 < len(lineas):
            meta["institution"] = lineas[i + 1]
    # Tras la presión arterial vienen fabricante, modelo, versión y software.
    j = next((k for k, l in enumerate(lineas) if "mmHg" in l), None)
    if j is not None and j + 4 < len(lineas):
        meta["manufacturer"], meta["model"] = lineas[j + 1], lineas[j + 2]
        meta["software"] = f"{lineas[j + 4]} (v{lineas[j + 3]})"
    meta["heart_rate_bpm"], = _num(r"(\d+)\s*/min", texto)
    meta["bp_systolic_mmhg"], meta["bp_diastolic_mmhg"] = _num(r"(\d+)\s*/\s*(\d+)\s*mmHg", texto, 2)
    meta["paper_speed_mm_s"], = _num(r"(\d+(?:\.\d+)?)\s*mm/s", texto)
    meta["gain_mm_mv"], = _num(r"(\d+(?:\.\d+)?)\s*mm/mV", texto)
    meta["filter_low_hz"], meta["filter_high_hz"] = _num(r"(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)\s*Hz", texto, 2)
    meta["notch_hz"], = _num(r"^(\d+)\s*Hz$", texto, flags=re.M)
    meta["qrs_ms"], = _num(r"\nQRS\n(\d+)\s*ms", texto)
    meta["pr_ms"], = _num(r"\nPQ\n(\d+)\s*ms", texto)
    meta["p_ms"], = _num(r"\nP\n(\d+)\s*ms", texto)
    meta["rr_ms"], meta["pp_ms"] = _num(r"RR/PP\n(\d+)\s*/\s*(\d+)\s*ms", texto, 2)
    meta["p_axis_deg"], meta["qrs_axis_deg"], meta["t_axis_deg"] = _num(
        r"P/QRS/T\n(-?\d+)\s*/\s*(-?\d+)\s*/\s*(-?\d+)", texto, 3
    )
    q = re.search(r"QT\s*/\s*(QTc\w*)\n(\d+)\s*/\s*(\d+)\s*ms", texto)
    if q:
        meta["qtc_formula"], meta["qt_ms"], meta["qtc_ms"] = q.group(1), float(q.group(2)), float(q.group(3))
    l = re.search(r"\n(\d+x\d+x\d+_\d+)\n", texto)
    if l:
        meta["lead_layout"] = l.group(1)
    s = re.search(r"\n([^\n]*)\nQRS\n", texto)
    if s:
        meta["report_status"] = s.group(1)
    r = re.search(r"P/QRS/T\n.*?\n:\n(.*?)\n\d+/\d+\s*$", texto, re.S)
    if r:
        meta["interpretation"] = " | ".join(x.strip() for x in r.group(1).split("\n") if x.strip())
    return meta


def read_pdf(pdf_path: Path, avisar=print) -> list[ECGRecord]:
    """Una ECGRecord por página vectorial. Las páginas sin 12 trazas se omiten con aviso."""
    doc = pymupdf.open(pdf_path)
    registros = []
    for n, page in enumerate(doc, 1):
        if not is_vector_page(page):
            avisar(f"  page {n}: not vectorized (images: {len(page.get_images())}); skipped")
            continue
        texto = page.get_text()
        meta = parse_text(texto)
        vel = meta.get("paper_speed_mm_s", np.nan)
        gan = meta.get("gain_mm_mv", np.nan)
        vel = 25.0 if np.isnan(vel) else vel
        gan = 10.0 if np.isnan(gan) else gan
        meta["paper_speed_mm_s"], meta["gain_mm_mv"] = vel, gan

        trazas = [_extraer_traza(d, vel, gan) for d in _trazas(page)]
        t_fin = min(t[-1] for t, _ in trazas)
        t_new = np.arange(0, t_fin, 1 / FS_OUT)
        mv = np.column_stack([np.interp(t_new, t, y) for t, y in trazas])
        t_new, mv = to_fixed_rate(t_new, mv, FS_OUT)

        meta.update(
            source_type="pdf", page=n, calibration_pulse_trimmed=False,
            extra_json=json.dumps({"pdf_text": texto.split("\n")}, ensure_ascii=False),
        )
        registros.append(ECGRecord(name=f"{pdf_path.stem}_p{n}", t=t_new, mv=mv, metadata=meta))
    return registros
