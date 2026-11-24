import argparse
import sys
from pathlib import Path

from ecg_digitizer.discovery import discover
from ecg_digitizer.dicom_reader import read_dicom
from ecg_digitizer.pdf_reader import read_pdf
from ecg_digitizer.writer import write_record


def main(argv=None):
    ap = argparse.ArgumentParser(description="Extract ECG waveforms from vector PDFs or DICOM files to 1000 Hz CSV.")
    ap.add_argument("-i", "--input", required=True, type=Path, help="folder with PDFs or DICOM files")
    ap.add_argument("-o", "--output", required=True, type=Path, help="output folder")
    ap.add_argument("--plot", action="store_true", help="save one validation PNG per ECG")
    args = ap.parse_args(argv)

    if not args.input.is_dir():
        ap.error(f"{args.input} is not a folder")

    pdfs, dicoms, ignorados = discover(args.input)
    print(f"Found: {len(pdfs)} PDF(s), {len(dicoms)} DICOM(s); ignored: {len(ignorados)}")
    for p in ignorados:
        print(f"  ignored: {p.name}")

    ok, fallos = 0, 0
    for dcm, csv in dicoms:
        print(f"[DICOM] {dcm.name}" + ("" if csv else "  (no companion CSV)"))
        try:
            for rec in read_dicom(dcm, csv):
                rec.metadata["source_file"] = dcm.name
                write_record(rec, args.output, args.plot)
                ok += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            fallos += 1

    for pdf in pdfs:
        print(f"[PDF] {pdf.name}")
        try:
            recs = read_pdf(pdf)
            for rec in recs:
                rec.metadata["source_file"] = pdf.name
                write_record(rec, args.output, args.plot)
            print(f"  extracted {len(recs)} ECG(s)")
            ok += len(recs)
        except Exception as e:
            print(f"  ERROR: {e}")
            fallos += 1

    print(f"ECGs exported: {ok}; files with errors: {fallos}; output: {args.output}")
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
