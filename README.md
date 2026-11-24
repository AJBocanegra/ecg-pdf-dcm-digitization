# ECG Digitizer

A small Python command-line tool for extracting 12-lead ECG waveforms from vector PDFs and DICOM waveform files. It writes each ECG as a standardized 1000 Hz CSV and writes a separate metadata CSV using the same schema for both input formats.

## Features

- Reads DICOM 12-lead waveform data directly from `.dcm` files. A same-basename `.csv` metadata dump may be present, but is not used as the waveform source.
- Reads vectorized ECG PDFs with one ECG per page. Each qualifying page becomes a separate record.
- Resamples waveform data to 1000 Hz.
- Exports signal values in microvolts and timestamps in seconds, using the `example_csv.csv` column layout.
- Exports ECG measurements and annotations, including QT and QTc when available. Missing values are written as `NaN`.
- Optionally creates one validation plot per ECG.

Raster-only PDF pages are skipped with a warning. PDF waveform extraction supports the tested GE MAC 2000, GE MAC VU360, and Schiller vector layouts; other PDF formats may need adjustments.

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

- `<record_name>.csv`: one row per sample. Columns are `I, II, III, AVR, AVL, AVF, V1, V2, V3, V4, V5, V6, TIME IN SEC.`. Lead values are in microvolts (µV) and `TIME IN SEC.` is in seconds. The number of rows depends on the duration of the available traces.
- `<record_name>_metadata.csv`: one row with a shared set of metadata fields for PDF and DICOM records. Measurements not available in the source are `NaN`. Additional annotations and source details are preserved in `extra_json` where possible.

For a multi-page PDF, output names include the page number, for example `<pdf_name>_p1.csv` and `<pdf_name>_p1_metadata.csv`.

## Input handling and limitations

- The input folder is scanned non-recursively. Unrecognized files and CSV files without a matching DICOM are ignored.
- PDF traces are detected from dark vector paths; the page must contain exactly 12 candidate paths. Lead names in the PDF text are matched to traces by position. If all 12 names cannot be matched unambiguously, the program warns and falls back to fixed spatial order (left to right, then top to bottom): `I, II, III, aVR, aVL, aVF, V1–V6`. A different layout can therefore be mislabeled; check the source and output.
- PDF vector coordinates are measured in points (72 points per inch, approximately 2.835 points/mm). Paper speed and gain are read from the text. If absent, the program warns and assumes 25 mm/s and 10 mm/mV respectively.
- The program looks for labeled, capped horizontal time and vertical voltage rulers on every PDF page. Where confidently detected, the ruler's measured points per second or points per mV takes precedence over the corresponding paper scale. A discrepancy greater than 10% produces a warning; conflicting rulers are ignored with a warning. Detection is geometric and cannot recognize every template.
- If no voltage ruler is found, the program uses a detected rectangular calibration pulse, or otherwise the printed/default gain. Pulses are assumed to represent **1 mV** if their value is not explicitly identified. A pulse inside a waveform path is removed before timing the signal. A single detected pulse supplies a common amplitude scale; multiple pulses are associated with nearby rows where possible. Disagreements with the ruler or paper gain produce warnings. The initial vertical markers in tested GE MAC2000 traces are removed from the signal, not used as the zero.
- In the tested Schiller PDF, the labeled 1 s and 1 mV rulers set the time and voltage scales. No calibration pulse is present in its waveform paths, so the signal start is not trimmed.
- DICOM channel sensitivity and correction factor first convert samples to mV. A recognizable, common, rectangular trailing pulse is assumed to be 1 mV; its height is used to scale each channel, then the pulse and preceding flat portion are removed. If no reliable pulse is detected, the program warns and retains the DICOM-scaled signal.
- After calibration, trimming, and resampling, **every PDF and DICOM lead is centered on its own median**. This changes the absolute baseline and does not guarantee that the median is the true physiological zero.
- In PDFs with several columns, each extracted lead is timed from the start of its own displayed trace. The CSV places them on the same time axis but does **not** establish that the columns represent simultaneous acquisition; do not assume synchrony. All leads are truncated to the shortest extracted trace. The output sampling rate is interpolated and does not imply the original acquisition rate was 1000 Hz.
- The PDF's original filtering has already been applied by the recording device. The metadata parser recognizes the tested GE MAC 2000 report layout; other templates may leave fields empty or require parser changes. PDF `extra_json` records the lead assignment, pulse count, detected ruler measurements, and chosen time and amplitude scale sources. For both PDF and DICOM it records the median baseline policy.
- This tool is for waveform extraction and data handling, not clinical interpretation or diagnosis. Validate extracted traces before any downstream use.

## Privacy

Do not commit patient recordings or generated outputs. ECG PDFs and DICOM files may contain protected health information, and PDF metadata output preserves page text in `extra_json`. Output signal and metadata CSV files may also be sensitive. The repository `.gitignore` excludes common input and output file types, but always inspect `git status` before committing and do not use `git add -f` on patient data.

## License

This project is licensed under the [MIT License](./LICENSE).

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
