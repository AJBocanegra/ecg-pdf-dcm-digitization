import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pymupdf

from .common import FS_OUT, LEADS, ECGRecord, center_by_median, to_fixed_rate

MM = 72 / 25.4  # puntos PDF por mm
MIN_SEGMENTOS = 100  # una traza de señal tiene miles de segmentos; los marcos y etiquetas, pocos


def _es_color_de_traza(color):
    return color is not None and max(color) <= 0.1


def _trazas(page):
    ds = [
        d for d in page.get_drawings()
        if _es_color_de_traza(d["color"]) and len(d["items"]) > MIN_SEGMENTOS
    ]
    # MAC VU360 places six leads in each of two columns. Other tested
    # layouts use one column, so this key preserves their vertical order.
    return sorted(ds, key=lambda d: (round(d["rect"].x0, 1), d["rect"].y0 + d["rect"].y1))


def is_vector_page(page):
    return len(_trazas(page)) == len(LEADS)


def _marca_inicial(d, ganancia):
    item = d["items"][0]
    if item[0] != "l" or len(d["items"]) < 2 or d["items"][1][0] != "l":
        return None
    a, b = item[1:]
    inicio = d["items"][1][1]
    altura = abs(b.y - a.y)
    if (abs(b.x - a.x) < 0.1 and abs(inicio.x - a.x) < 0.1
            and abs(inicio.y - b.y) > 0.25
            and 0.8 * ganancia * MM <= altura <= 1.2 * ganancia * MM):
        return (a.y + b.y) / 2
    return None


def _ordenar_por_etiquetas(page, dibujos, ganancia):
    etiquetas = []
    nombres = {lead.upper(): lead for lead in LEADS}
    for bloque in page.get_text("dict")["blocks"]:
        for linea in bloque.get("lines", []):
            texto = "".join(span["text"] for span in linea["spans"])
            nombre = nombres.get(re.sub(r"\s+", "", texto).upper())
            if nombre:
                x0, y0, _, y1 = linea["bbox"]
                etiquetas.append((nombre, x0, (y0 + y1) / 2))
    asignados = {}
    for d in dibujos:
        items = d["items"]
        marca = _marca_inicial(d, ganancia)
        inicio = items[1][1] if marca is not None else items[0][1]
        y = marca if marca is not None else inicio.y
        candidatos = sorted(
            (abs(ey - y), nombre) for nombre, ex, ey in etiquetas
            if abs(ex - inicio.x) <= 25 and abs(ey - y) <= 30
        )
        if not candidatos or (len(candidatos) > 1 and candidatos[1][0] - candidatos[0][0] < 2):
            return None
        nombre = candidatos[0][1]
        if nombre in asignados:
            return None
        asignados[nombre] = d
    if set(asignados) != set(LEADS):
        return None
    return [asignados[lead] for lead in LEADS]


def _pulsos(page):
    pulsos = []
    for d in page.get_drawings():
        if not _es_color_de_traza(d.get("color")):
            continue
        items = d["items"]
        for i in range(len(items) - 4):
            seg = items[i:i + 5]
            if any(item[0] != "l" for item in seg):
                continue
            p = [seg[0][1]] + [item[2] for item in seg]
            if any(abs(seg[j][2].x - seg[j + 1][1].x) > 0.2 or
                   abs(seg[j][2].y - seg[j + 1][1].y) > 0.2 for j in range(4)):
                continue
            height = p[0].y - p[2].y
            width = p[3].x - p[2].x
            if (5 <= height <= 60 and 3 <= width <= 30
                    and 0 <= p[1].x - p[0].x <= 8 and 0 <= p[5].x - p[4].x <= 8
                    and abs(p[0].y - p[1].y) < 0.3
                    and abs(p[1].x - p[2].x) < 0.3
                    and abs(p[2].y - p[3].y) < 0.3
                    and abs(p[3].x - p[4].x) < 0.3
                    and abs(p[4].y - p[5].y) < 0.3
                    and abs(p[0].y - p[4].y) < 0.3):
                pulsos.append({"x_fin": p[5].x, "y_cero": (p[0].y + p[4].y) / 2,
                               "altura": height, "valor_mv": 1.0, "valor_explicito": False})
    for bloque in page.get_text("dict")["blocks"]:
        for linea in bloque.get("lines", []):
            texto = "".join(span["text"] for span in linea["spans"]).strip()
            m = re.fullmatch(r"(\d+(?:[.,]\d+)?)\s*mV", texto, flags=re.I)
            if not m or not pulsos:
                continue
            x0, y0, x1, y1 = linea["bbox"]
            x, y = (x0 + x1) / 2, (y0 + y1) / 2
            cercano = min(pulsos, key=lambda p: (p["x_fin"] - x) ** 2 + (p["y_cero"] - y) ** 2)
            if abs(cercano["x_fin"] - x) <= 40 and abs(cercano["y_cero"] - y) <= 40:
                valor = float(m.group(1).replace(",", "."))
                if valor <= 0:
                    raise ValueError("calibration pulse value must be positive")
                cercano["valor_mv"] = valor
                cercano["valor_explicito"] = True
    return pulsos


def _reglas(page, avisar):
    lineas = []
    for d in page.get_drawings():
        color = d.get("color")
        if color is None or max(color) > 0.2 or len(d["items"]) > 12:
            continue
        lineas.extend((i[1], i[2]) for i in d["items"] if i[0] == "l")

    def es_marca(a, b, x, y, horizontal):
        dx, dy = abs(b.x - a.x), abs(b.y - a.y)
        if horizontal:
            return (dx < 0.5 and 1 <= dy <= 10 and abs(a.x - x) < 1
                    and min(a.y, b.y) - 1 <= y <= max(a.y, b.y) + 1)
        return (dy < 0.5 and 1 <= dx <= 10 and abs(a.y - y) < 1
                and min(a.x, b.x) - 1 <= x <= max(a.x, b.x) + 1)

    candidatos = {"time": [], "amplitude": []}
    for bloque in page.get_text("dict")["blocks"]:
        for linea in bloque.get("lines", []):
            texto = "".join(span["text"] for span in linea["spans"]).strip()
            match = re.fullmatch(r"(\d+(?:[.,]\d+)?)\s*(mV|s)", texto, re.I)
            if not match:
                continue
            valor = float(match.group(1).replace(",", "."))
            if valor <= 0:
                avisar(f"  invalid calibration ruler label: {texto}")
                continue
            horizontal = match.group(2).lower() == "s"
            tipo = "time" if horizontal else "amplitude"
            x0, y0, x1, y1 = linea["bbox"]
            x, y = (x0 + x1) / 2, (y0 + y1) / 2
            for a, b in lineas:
                dx, dy = abs(b.x - a.x), abs(b.y - a.y)
                length = dx if horizontal else dy
                if (not 12 <= length <= 200
                        or (dy if horizontal else dx) >= 0.5
                        or abs((a.x + b.x) / 2 - x) > 20
                        or abs((a.y + b.y) / 2 - y) > 20):
                    continue
                if all(any(es_marca(c, d, p.x, p.y, horizontal) for c, d in lineas)
                       for p in (a, b)):
                    candidatos[tipo].append(length / valor)
    resultado = {}
    for tipo, valores in candidatos.items():
        if valores and np.allclose(valores, valores[0], rtol=0.01):
            resultado[tipo] = float(np.mean(valores))
        elif valores:
            avisar(f"  conflicting {tipo} calibration rulers; ignoring them")
    return resultado


def _extraer_traza(d, velocidad, ganancia, escala=None, fin_pulso=None, escala_tiempo=None):
    marca = _marca_inicial(d, ganancia)
    items = d["items"][1:] if marca is not None else d["items"]
    if any(item[0] != "l" for item in items):
        raise ValueError("unsupported waveform path: expected line segments")
    pts = [(i[1].x, i[1].y) for i in items] + [(items[-1][2].x, items[-1][2].y)]
    x, y = np.array(pts).T
    if fin_pulso is not None:
        keep = x >= fin_pulso
        x, y = x[keep], y[keep]
    if len(x) < 2:
        raise ValueError("waveform has fewer than two points after calibration")
    x_u, inv = np.unique(x, return_inverse=True)
    y_u = np.bincount(inv, weights=y) / np.bincount(inv)
    t = (x_u - x_u[0]) / (escala_tiempo if escala_tiempo is not None else MM * velocidad)
    escala = escala if escala is not None else MM * ganancia
    mv = -y_u / escala
    return t, mv, fin_pulso is not None


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
        if np.isnan(vel):
            avisar(f"  page {n}: paper speed not found; assuming 25 mm/s")
        if np.isnan(gan):
            avisar(f"  page {n}: gain not found; assuming 10 mm/mV")
        vel = 25.0 if np.isnan(vel) else vel
        gan = 10.0 if np.isnan(gan) else gan
        if vel <= 0 or gan <= 0:
            raise ValueError(f"page {n}: paper speed and gain must be positive")
        meta["paper_speed_mm_s"], meta["gain_mm_mv"] = vel, gan

        dibujos = _trazas(page)
        ordenados = _ordenar_por_etiquetas(page, dibujos, gan)
        asignacion = "labels" if ordenados is not None else "fixed_order"
        if ordenados is None:
            avisar(f"  page {n}: lead labels incomplete or ambiguous; using fixed spatial order")
            ordenados = dibujos
        pulsos = _pulsos(page)
        reglas = _reglas(page, lambda mensaje: avisar(f"  page {n}: {mensaje}"))
        for tipo, esperado in (("time", MM * vel), ("amplitude", MM * gan)):
            if tipo in reglas and abs(reglas[tipo] / esperado - 1) > 0.1:
                avisar(f"  page {n}: {tipo} ruler disagrees with paper scale; using ruler")
        if not pulsos and "amplitude" not in reglas:
            avisar(f"  page {n}: no calibration pulse or amplitude ruler found; using paper gain")
        trazas = []
        recortados = []
        escalas = []
        for d in ordenados:
            marca = _marca_inicial(d, gan)
            inicio = d["items"][1][1] if marca is not None else d["items"][0][1]
            cercanos = [p for p in pulsos
                        if p["x_fin"] <= inicio.x + 20 * MM and abs(p["y_cero"] - inicio.y) < 20]
            pulso = min(cercanos, key=lambda p: abs(p["y_cero"] - inicio.y)) if cercanos else None
            referencia = pulso if pulso else pulsos[0] if len(pulsos) == 1 else None
            escala = reglas.get("amplitude")
            fuente = "ruler"
            if escala is None:
                escala = referencia["altura"] / referencia["valor_mv"] if referencia else None
                fuente = "pulse" if escala is not None else "paper_gain"
            fin = pulso["x_fin"] if pulso and inicio.x < pulso["x_fin"] else None
            t, mv, recortado = _extraer_traza(d, vel, gan, escala, fin, reglas.get("time"))
            trazas.append((t, mv))
            recortados.append(recortado)
            escalas.append(fuente)
        if pulsos and "paper_gain" in escalas:
            avisar(f"  page {n}: no matching pulse for some leads; using paper gain")
        if pulsos and any(abs(p["altura"] / p["valor_mv"] / MM - gan) > gan * 0.1 for p in pulsos):
            avisar(f"  page {n}: calibration pulse disagrees with paper gain; using ruler or pulse scale")
        if "amplitude" in reglas and any(
                abs(p["altura"] / p["valor_mv"] / reglas["amplitude"] - 1) > 0.1 for p in pulsos):
            avisar(f"  page {n}: calibration pulse disagrees with amplitude ruler; using ruler")
        t_fin = min(t[-1] for t, _ in trazas)
        t_new = np.arange(0, t_fin, 1 / FS_OUT)
        mv = np.column_stack([np.interp(t_new, t, y) for t, y in trazas])
        t_new, mv = to_fixed_rate(t_new, mv, FS_OUT)
        mv = center_by_median(mv)

        meta.update(
            source_type="pdf", page=n, calibration_pulse_trimmed=any(recortados),
            extra_json=json.dumps({
                "pdf_text": texto.split("\n"),
                "lead_assignment": asignacion,
                "calibration_pulses": len(pulsos),
                "pulse_values_mv": [p["valor_mv"] for p in pulsos],
                "pulse_values_assumed": [not p["valor_explicito"] for p in pulsos],
                "baseline": "per_lead_median",
                "ruler_points_per_second": reglas.get("time"),
                "ruler_points_per_mv": reglas.get("amplitude"),
                "time_scale": "ruler" if "time" in reglas else "paper_speed",
                "amplitude_scale": dict(zip(LEADS, escalas)),
            }, ensure_ascii=False),
        )
        registros.append(ECGRecord(name=f"{pdf_path.stem}_p{n}", t=t_new, mv=mv, metadata=meta))
    return registros
