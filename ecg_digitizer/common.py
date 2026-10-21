from dataclasses import dataclass, field

import numpy as np

LEADS = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
FS_OUT = 500  # Hz

# Esquema único de metadatos, independiente del origen (PDF o DICOM).
METADATA_COLUMNS = [
    "source_type", "source_file", "page", "signal_file",
    "acquisition_datetime", "institution", "manufacturer", "model", "software",
    "patient_sex",
    "source_sampling_hz", "output_sampling_hz", "n_samples", "duration_s",
    "paper_speed_mm_s", "gain_mm_mv", "filter_low_hz", "filter_high_hz", "notch_hz",
    "lead_layout", "report_status", "calibration_pulse_trimmed",
    "heart_rate_bpm", "bp_systolic_mmhg", "bp_diastolic_mmhg",
    "rr_ms", "pp_ms", "pr_ms", "p_ms", "qrs_ms", "qt_ms", "qtc_ms", "qtc_formula",
    "p_axis_deg", "qrs_axis_deg", "t_axis_deg",
    "interpretation", "extra_json",
]


@dataclass
class ECGRecord:
    name: str
    t: np.ndarray  # segundos
    mv: np.ndarray  # (n_muestras, 12) en mV, columnas en orden LEADS
    metadata: dict = field(default_factory=dict)


def to_fixed_rate(t, mv, fs=FS_OUT):
    """Remuestrea a `fs` Hz; si ya está a esa tasa devuelve los datos sin tocar."""
    if len(t) > 1 and np.allclose(np.diff(t), 1 / fs, atol=1e-6):
        return t - t[0], mv
    t = t - t[0]
    t_new = np.arange(0, t[-1], 1 / fs)
    out = np.column_stack([np.interp(t_new, t, mv[:, j]) for j in range(mv.shape[1])])
    return t_new, out
