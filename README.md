# ECG Digitizer

A small Python command-line tool for extracting 12-lead ECG waveforms from vector PDFs and DICOM waveform files. It writes each ECG as a standardized 500 Hz CSV and writes a separate metadata CSV using the same schema for both input formats.

## Features

- Reads DICOM 12-lead waveform data directly from `.dcm` files. A same-basename `.csv` metadata dump may be present, but is not used as the waveform source.
- Reads vectorized ECG PDFs with one ECG per page. Each qualifying page becomes a separate record.
- Resamples waveform data to 500 Hz.
- Exports signal values in millivolts and timestamps in seconds.
- Exports ECG measurements and annotations, including QT and QTc when available. Missing values are written as `NaN`.
- Optionally creates one validation plot per ECG.

Raster-only PDF pages are skipped with a warning. PDF waveform extraction is currently tuned to the GE MAC 2000 vector layout used during development; other PDF formats may need adjustments.

## Requirements

- Python 3.10 or newer
- pip

## Installation

Clone or download this repository, then run these commands from the project directory. On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Usage

Place input files in a folder. The tool scans that folder (not its subfolders) for PDFs and DICOM files. DICOM companion CSV files are matched by basename and are optional.

```bash
python main.py --input /path/to/input --output /path/to/output
```

Short options are also available:

```bash
python main.py -i /path/to/input -o /path/to/output
```

Add `--plot` to write validation PNGs alongside the CSV files:

```bash
python main.py -i /path/to/input -o /path/to/output --plot
```

Each ECG produces two CSV files:

- `<record_name>.csv`: one row per sample. Columns are `I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6, time_stamp`. Lead values are in mV and `time_stamp` is in seconds.
- `<record_name>_metadata.csv`: one row with a shared set of metadata fields for PDF and DICOM records. Measurements not available in the source are `NaN`. Additional annotations and source details are preserved in `extra_json` where possible.

For a multi-page PDF, output names include the page number, for example `<pdf_name>_p1.csv` and `<pdf_name>_p1_metadata.csv`.

## Input handling and limitations

- The input folder is scanned non-recursively. Unrecognized files and CSV files without a matching DICOM are ignored.
- DICOM waveform samples are scaled using the DICOM channel sensitivity and correction factor. A trailing calibration pulse is removed when detected.
- PDF traces are detected from their vector paths. The PDF calibration (paper speed and gain) is read from page text when available, with defaults of 25 mm/s and 10 mm/mV.
- PDF leads are median-centered because the vector trace does not preserve an absolute baseline. The PDF's original filtering has already been applied by the recording device.
- The metadata parser recognizes the tested GE MAC 2000 report layout. A different PDF manufacturer or template may leave some metadata fields empty or require parser changes.
- This tool is for waveform extraction and data handling, not clinical interpretation or diagnosis. Validate extracted traces before any downstream use.

## Privacy

Do not commit patient recordings or generated outputs. ECG PDFs and DICOM files may contain protected health information, and PDF metadata output preserves page text in `extra_json`. Output signal and metadata CSV files may also be sensitive. The repository `.gitignore` excludes common input and output file types, but always inspect `git status` before committing and do not use `git add -f` on patient data.

## Project layout

```text
ecg_digitizer/
  common.py       Shared record model, lead order, and resampling
  discovery.py    Input file discovery
  dicom_reader.py DICOM waveform and annotation extraction
  pdf_reader.py   Vector PDF signal and metadata extraction
  writer.py       Standardized signal/metadata CSV and optional plots
main.py           Command-line entry point
requirements.txt  Python dependencies
```
