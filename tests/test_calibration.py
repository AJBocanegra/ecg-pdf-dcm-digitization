import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf

from ecg_digitizer.common import LEADS, ECGRecord, center_by_median
from ecg_digitizer.dicom_reader import recortar_calibracion
from ecg_digitizer.pdf_reader import MM, _extraer_traza, _ordenar_por_etiquetas, _pulsos, _reglas
from ecg_digitizer.writer import write_record


def line(a, b):
    return ("l", pymupdf.Point(*a), pymupdf.Point(*b))


class Page:
    def __init__(self, drawings, labels=()):
        self.drawings = drawings
        self.labels = labels

    def get_drawings(self):
        return self.drawings

    def get_text(self, kind):
        return {"blocks": [{"lines": [
            {"spans": [{"text": name}], "bbox": (x, y, x + 5, y + 4)}
            for name, x, y in self.labels
        ]}]}


def drawing(items):
    return {"items": items, "color": (0.0, 0.0, 0.0)}


class PdfCalibrationTests(unittest.TestCase):
    def test_reordered_labels_override_spatial_order(self):
        names = LEADS[::-1]
        drawings = [drawing([line((10, 50 * i), (11, 50 * i))])
                    for i in range(12)]
        labels = [(name, 10, 50 * i - 2) for i, name in enumerate(names)]
        ordered = _ordenar_por_etiquetas(Page(drawings, labels), drawings, 10)
        self.assertEqual(ordered, drawings[::-1])

    def test_partial_or_ambiguous_labels_fall_back(self):
        drawings = [drawing([line((10, 50 * i), (11, 50 * i))])
                    for i in range(12)]
        labels = [(name, 10, 50 * i - 2) for i, name in enumerate(LEADS)]
        self.assertIsNone(_ordenar_por_etiquetas(Page(drawings, labels[:-1]), drawings, 10))
        labels[1] = ("I", 10, 48)
        self.assertIsNone(_ordenar_por_etiquetas(Page(drawings, labels), drawings, 10))

    def test_pulse_inside_path_is_removed_and_used_for_scale(self):
        height = 10 * MM
        pts = [(0, 100), (2, 100), (2, 100 - height),
               (15, 100 - height), (15, 100), (17, 100),
               (17 + MM, 100 - height / 2), (17 + 2 * MM, 100)]
        d = drawing([line(a, b) for a, b in zip(pts, pts[1:])])
        pulse, = _pulsos(Page([d]))
        t, mv, trimmed = _extraer_traza(d, 25, 10, pulse["altura"], pulse["x_fin"])
        np.testing.assert_allclose(t, [0, 0.04, 0.08])
        np.testing.assert_allclose(mv - mv[0], [0, 0.5, 0])
        self.assertTrue(trimmed)

    def test_explicit_pulse_value_is_read_from_nearby_text(self):
        height = 20 * MM
        pts = [(0, 100), (2, 100), (2, 100 - height),
               (15, 100 - height), (15, 100), (17, 100)]
        d = drawing([line(a, b) for a, b in zip(pts, pts[1:])])
        pulse, = _pulsos(Page([d], [("2 mV", 18, 98)]))
        self.assertEqual(pulse["valor_mv"], 2)
        self.assertTrue(pulse["valor_explicito"])
        self.assertAlmostEqual(pulse["altura"] / pulse["valor_mv"], 10 * MM)

    def test_waveform_first_segment_is_not_removed_without_marker(self):
        d = drawing([line((0, 100), (MM, 100 - MM * 10)),
                     line((MM, 100 - MM * 10), (2 * MM, 100))])
        t, mv, trimmed = _extraer_traza(d, 25, 10)
        np.testing.assert_allclose(t, [0, 0.04, 0.08])
        self.assertFalse(trimmed)
        self.assertEqual(len(mv), 3)

    def test_vertical_marker_is_removed_without_dropping_signal(self):
        d = drawing([line((0, 100 + 5 * MM), (0, 100 - 5 * MM)),
                     line((0, 100), (MM, 100 - 5 * MM))])
        t, mv, trimmed = _extraer_traza(d, 25, 10)
        np.testing.assert_allclose(t, [0, 0.04])
        np.testing.assert_allclose(mv - mv[0], [0, 0.5])
        self.assertFalse(trimmed)

    def test_labeled_rulers_override_both_paper_scales(self):
        rulers = [
            drawing([line((0, 100), (100, 100))]),
            drawing([line((0, 98), (0, 102))]),
            drawing([line((100, 98), (100, 102))]),
            drawing([line((0, 40), (0, 80))]),
            drawing([line((-2, 40), (2, 40))]),
            drawing([line((-2, 80), (2, 80))]),
        ]
        page = Page(rulers, [("2 s", 47.5, 93), ("2 mV", 4, 58)])
        scales = _reglas(page, self.fail)
        self.assertEqual(scales, {"time": 50, "amplitude": 20})
        trace = drawing([line((0, 50), (50, 30)), line((50, 30), (100, 50))])
        t, mv, _ = _extraer_traza(trace, 25, 10, scales["amplitude"],
                                   escala_tiempo=scales["time"])
        np.testing.assert_allclose(t, [0, 1, 2])
        np.testing.assert_allclose(mv - mv[0], [0, 1, 0])

    def test_unmarked_line_is_not_a_ruler(self):
        page = Page([drawing([line((0, 100), (100, 100))])], [("1 s", 47.5, 93)])
        self.assertEqual(_reglas(page, self.fail), {})

    def test_conflicting_rulers_are_ignored_with_warning(self):
        lines = []
        for y, length in [(100, 60), (200, 80)]:
            lines.extend([drawing([line((0, y), (length, y))]),
                          drawing([line((0, y - 2), (0, y + 2))]),
                          drawing([line((length, y - 2), (length, y + 2))])])
        warnings = []
        page = Page(lines, [("1 s", 27.5, 93), ("1 s", 37.5, 193)])
        self.assertEqual(_reglas(page, warnings.append), {})
        self.assertEqual(len(warnings), 1)


class MedianTests(unittest.TestCase):
    def test_per_lead_centering_preserves_nan_and_morphology(self):
        mv = np.array([[10., -4., np.nan], [11., 0., np.nan], [14., 2., np.nan]])
        centered = center_by_median(mv)
        np.testing.assert_allclose(centered[:, :2], [[-1., -4.], [0., 0.], [3., 2.]])
        self.assertTrue(np.isnan(centered[:, 2]).all())
        np.testing.assert_allclose(np.diff(centered[:, :2], axis=0), np.diff(mv[:, :2], axis=0))


class DicomCalibrationTests(unittest.TestCase):
    def test_rectangular_pulse_scales_each_channel_and_removes_tail(self):
        t = np.arange(400) / 1000
        signal = np.tile(np.linspace(-0.2, 0.2, 300)[:, None], (1, 12))
        baseline = np.arange(12) * 0.01
        height = np.linspace(0.7, 1.3, 12)
        mv = np.vstack([signal, np.tile(baseline, (30, 1)),
                        np.tile(baseline + height, (40, 1)),
                        np.tile(baseline, (30, 1))])
        t_out, mv_out, trimmed = recortar_calibracion(t, mv)
        self.assertTrue(trimmed)
        self.assertEqual(len(t_out), 300)
        np.testing.assert_allclose(mv_out, (signal - baseline) / height)

    def test_nonrectangular_spike_is_not_treated_as_calibration(self):
        t = np.arange(400) / 1000
        mv = np.zeros((400, 12))
        mv[350] = 1
        t_out, mv_out, trimmed = recortar_calibracion(t, mv)
        self.assertFalse(trimmed)
        np.testing.assert_array_equal(t_out, t)
        np.testing.assert_array_equal(mv_out, mv)


class OutputTests(unittest.TestCase):
    def test_template_header_and_microvolts(self):
        rec = ECGRecord("test", np.array([0.0, 0.001]),
                        np.array([[0.5] * 12, [-0.25] * 12]))
        with tempfile.TemporaryDirectory() as tmp:
            path = write_record(rec, Path(tmp))
            out = pd.read_csv(path)
        self.assertEqual(list(out.columns),
                         ["I", "II", "III", "AVR", "AVL", "AVF",
                          "V1", "V2", "V3", "V4", "V5", "V6", "TIME IN SEC."])
        np.testing.assert_allclose(out["I"], [500, -250])
        np.testing.assert_allclose(out["TIME IN SEC."], [0, 0.001])


if __name__ == "__main__":
    unittest.main()
