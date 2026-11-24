import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pydicom

from .common import FS_OUT, LEADS, ECGRecord, center_by_median, to_fixed_rate

ANNOTATION_FIELDS = {
    "RR Interval": "rr_ms",
    "PR Interval": "pr_ms",
    "QRS Duration": "qrs_ms",
    "QT Interval": "qt_ms",
    "QTc Interval": "qtc_ms",
    "Ventricular Heart Rate": "heart_rate_bpm",
    "P Axis": "p_axis_deg",
    "QRS Axis": "qrs_axis_deg",
    "T Axis": "t_axis_deg",
}


def _lead_name(texto):
    return re.sub(r"^lead\s*", "", texto.strip(), flags=re.I).replace(" ", "")


def _leer_grupo(g):
    fs = float(g.SamplingFrequency)
    n_ch, n_s = g.NumberOfWaveformChannels, g.NumberOfWaveformSamples
    assert g.WaveformBitsAllocated == 16 and g.WaveformSampleInterpretation == "SS"
    crudo = np.frombuffer(g.WaveformData, dtype="<i2").reshape(n_s, n_ch)  # muestras intercaladas
    etiquetas, cols = [], []
    for j, c in enumerate(g.ChannelDefinitionSequence):
        etiquetas.append(_lead_name(c.ChannelSourceSequence[0].CodeMeaning))
        escala = float(c.ChannelSensitivity) * float(c.get("ChannelSensitivityCorrectionFactor", 1)) / 1000  # uV -> mV
        cols.append((crudo[:, j] - float(c.get("ChannelBaseline", 0))) * escala)
    por_lead = dict(zip(etiquetas, cols))
    mv = np.column_stack([por_lead.get(l, np.full(n_s, np.nan)) for l in LEADS])
    return fs, mv, g.ChannelDefinitionSequence[0]


def recortar_calibracion(t, mv):
    """Detecta un pulso rectangular final común, supuesto de 1 mV, y calibra cada canal."""
    if len(t) < 3 or not np.isfinite(mv).all():
        return t, mv, False
    dt = np.median(np.diff(t))
    n_min = max(3, int(round(0.02 / dt)))
    d = np.diff(mv, axis=0)
    subidas = np.flatnonzero((d > 0.5).all(axis=1)) + 1
    for subida in subidas[::-1]:
        if subida < len(mv) // 2 or subida < n_min:
            continue
        bajadas = np.flatnonzero((d[subida:] < -0.5).all(axis=1)) + subida + 1
        for bajada in bajadas:
            if bajada - subida < n_min or len(mv) - bajada < n_min:
                continue
            cero = np.median(mv[subida - n_min:subida], axis=0)
            alto = np.median(mv[subida:bajada], axis=0)
            altura = alto - cero
            if not (np.all((altura > 0.5) & (altura < 2.0))
                    and np.allclose(mv[subida - n_min:subida], cero, atol=0.01)
                    and np.allclose(mv[subida:bajada], alto, atol=0.01)
                    and np.allclose(mv[bajada:], cero, atol=0.01)):
                continue
            inicio = subida
            while inicio > 0 and np.allclose(mv[inicio - 1], cero, atol=0.01):
                inicio -= 1
            if inicio == 0:
                continue
            return t[:inicio], (mv[:inicio] - cero) / altura, True
    return t, mv, False


def _fecha(valor):
    try:
        return datetime.strptime(str(valor)[:14], "%Y%m%d%H%M%S").isoformat()
    except ValueError:
        return np.nan


def _anotaciones(ds):
    campos, otras, fiduciales, textos = {}, {}, {}, []
    for a in ds.get("WaveformAnnotationSequence", []):
        nombre = a.ConceptNameCodeSequence[0].CodeMeaning if "ConceptNameCodeSequence" in a else ""
        if "NumericValue" in a:
            valor = float(a.NumericValue)
            if nombre in ANNOTATION_FIELDS:
                campos[ANNOTATION_FIELDS[nombre]] = valor
            else:
                unidad = a.MeasurementUnitsCodeSequence[0].CodeValue if "MeasurementUnitsCodeSequence" in a else ""
                otras[nombre] = f"{valor:g} {unidad}".strip()
        elif "ReferencedSamplePositions" in a:
            fiduciales[nombre] = [int(x) for x in np.atleast_1d(a.ReferencedSamplePositions)]
        elif a.get("UnformattedTextValue", "").strip():
            textos.append(a.UnformattedTextValue.strip())
    return campos, otras, fiduciales, textos


def read_dicom(dcm_path: Path, csv_path: Path | None = None) -> list[ECGRecord]:
    ds = pydicom.dcmread(dcm_path)
    grupos = [_leer_grupo(g) + (g.get("MultiplexGroupLabel", ""),) for g in ds.WaveformSequence]
    fs, mv, canal, etiqueta = max(grupos, key=lambda x: len(x[1]))
    t = np.arange(len(mv)) / fs
    t, mv, recortado = recortar_calibracion(t, mv)
    if not recortado:
        print("  no reliable trailing calibration pulse found; using DICOM channel sensitivities")
    t, mv = to_fixed_rate(t, mv, FS_OUT)
    mv = center_by_median(mv)

    campos, otras, fiduciales, textos = _anotaciones(ds)
    extra = {"waveform_group": etiqueta, "other_measurements": otras,
             "fiducials_samples": fiduciales, "baseline": "per_lead_median"}
    if csv_path:
        extra["paired_dicom_dump_csv"] = csv_path.name
    meta = {
        "source_type": "dicom",
        "acquisition_datetime": _fecha(ds.get("AcquisitionDateTime", "")),
        "institution": ds.get("InstitutionName") or np.nan,
        "manufacturer": ds.get("Manufacturer") or np.nan,
        "model": ds.get("ManufacturerModelName") or np.nan,
        "software": str(ds.get("SoftwareVersions", "")) or np.nan,
        "patient_sex": ds.get("PatientSex") or np.nan,
        "source_sampling_hz": fs,
        "filter_low_hz": float(canal.get("FilterLowFrequency", np.nan)),
        "filter_high_hz": float(canal.get("FilterHighFrequency", np.nan)),
        "notch_hz": float(canal.get("NotchFilterFrequency", np.nan)),
        "calibration_pulse_trimmed": recortado,
        "interpretation": " | ".join(textos) or np.nan,
        "extra_json": json.dumps(extra, ensure_ascii=False),
        **campos,
    }
    return [ECGRecord(name=dcm_path.stem, t=t, mv=mv, metadata=meta)]
