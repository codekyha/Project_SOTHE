"""p2_write_csv / p2_write_json ports: re-writing the MATLAB R2025b products reproduces them byte for byte."""
import json
import os
import tempfile
import unittest

import _util  # noqa: F401
import numpy as np

from sothe_p2 import phase2


class TestWriters(unittest.TestCase):
    def test_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            for name in phase2.CSV_FILES:
                src = os.path.join(_util.REF, name)
                raw = _util.text(src)
                head, M = _util.read_csv(src)
                out = os.path.join(d, name)
                phase2.write_csv(out, head, M)
                self.assertEqual(_util.text(out), raw, name)

    def test_json_roundtrip(self):
        raw = _util.text(os.path.join(_util.REF, "phase2_summary.json"))
        self.assertEqual(phase2.json_text(json.loads(raw)), raw)

    def test_json_rules(self):
        s = phase2.json_text({"a": 2.0, "b": 0.1, "c": [1, 2.5], "d": "q\"x", "e": True, "f": float("nan"), "g": {"h": [7]},
                              "i": np.array([[1, 2], [3, 4]])})
        self.assertEqual(s, '{\n  "a": 2,\n  "b": 0.10000000000000001,\n  "c": [1, 2.5],\n  "d": "q\\"x",\n  "e": true,\n'
                            '  "f": null,\n  "g": {\n    "h": 7\n  },\n  "i": [1, 3, 2, 4]\n}\n')


if __name__ == "__main__":
    unittest.main()
