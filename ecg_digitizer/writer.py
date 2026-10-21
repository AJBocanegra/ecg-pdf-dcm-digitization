from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import FS_OUT, LEADS, METADATA_COLUMNS, ECGRecord


def write_record(rec: ECGRecord, out_dir: Path, plot=False):
    out_dir.mkdir(parents=True, exist_ok=True)
    signal = pd.DataFrame(rec.mv, columns=LEADS)
    signal["time_stamp"] = rec.t
    signal_path = out_dir / f"{rec.name}.csv"
    signal.to_csv(signal_path, index=False, float_format="%.5f")

    meta = {c: np.nan for c in METADATA_COLUMNS}
    meta.update(rec.metadata)
    meta.update(
        source_file=meta["source_file"] if not pd.isna(meta["source_file"]) else rec.name,
        signal_file=signal_path.name, output_sampling_hz=FS_OUT,
        n_samples=len(rec.t), duration_s=round(float(rec.t[-1]), 3),
    )
    pd.DataFrame([meta], columns=METADATA_COLUMNS).to_csv(out_dir / f"{rec.name}_metadata.csv", index=False)

    if plot:
        fig, axs = plt.subplots(len(LEADS), 1, figsize=(14, 18), sharex=True)
        for ax, lead in zip(axs, LEADS):
            ax.plot(signal["time_stamp"], signal[lead], lw=0.6, color="k")
            ax.set_ylabel(lead, rotation=0, labelpad=20)
            ax.grid(alpha=0.3)
        axs[-1].set_xlabel("s")
        fig.tight_layout()
        fig.savefig(out_dir / f"{rec.name}.png", dpi=70)
        plt.close(fig)
    return signal_path
