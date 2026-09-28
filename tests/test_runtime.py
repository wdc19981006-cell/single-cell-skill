"""Runtime failure classification and process isolation without crashing R."""
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "runtime/r45"))
import run_r45 as launcher


class RuntimeTests(unittest.TestCase):
    def test_native_crash_even_after_pass_and_original_error_code(self):
        for code in (139, -1073741819, 3221225477, 7):
            with self.subTest(code=code), tempfile.TemporaryDirectory() as temporary:
                def process(command, **kwargs):
                    kwargs["stdout"].write(b"PASS\n")
                    class Child:
                        def wait(self): return code
                    return Child()
                with patch.object(launcher, "RUNTIME_LOG", Path(temporary)), \
                     patch.object(launcher.subprocess, "Popen", side_effect=process):
                    record = launcher.run(expression="1+1", stage="build")
                self.assertEqual(record["exit_code"], code)
                self.assertEqual(record["stage"], "build")
                self.assertEqual(record["status"], "R_SCRIPT_ERROR" if code == 7 else "R_NATIVE_CRASH")

    def test_missing_runtime_never_falls_back(self):
        with tempfile.TemporaryDirectory() as temporary, \
             patch.object(launcher, "RUNTIME_LOG", Path(temporary)), \
             patch.object(launcher, "RSCRIPT", Path(temporary) / "absent.exe"), \
             patch.object(launcher.subprocess, "Popen") as process:
            record = launcher.run(expression="1+1")
        self.assertEqual(record["status"], "R_RUNTIME_MISSING")
        self.assertNotEqual(record["exit_code"], 0)
        process.assert_not_called()

    def test_child_environment_is_fixed_without_mutating_parent(self):
        captured = {}
        def process(command, **kwargs):
            captured.update(kwargs["env"])
            self.assertIn("--no-environ", command)
            self.assertIn("--no-site-file", command)
            class Child:
                def wait(self): return 0
            return Child()
        with tempfile.TemporaryDirectory() as temporary, \
             patch.dict(os.environ, {"R_LIBS": "wrong-library", "LC_ALL": "C.UTF-8"}), \
             patch.object(launcher, "RUNTIME_LOG", Path(temporary)), \
             patch.object(launcher.subprocess, "Popen", side_effect=process):
            launcher.run(expression="1+1")
            self.assertEqual(os.environ["R_LIBS"], "wrong-library")
            self.assertEqual(os.environ["LC_ALL"], "C.UTF-8")
        self.assertEqual(captured["CLI_NO_THREAD"], "1")
        self.assertEqual(captured["LC_ALL"], "C")
        self.assertEqual(captured["R_LIBS"], str(launcher.R_LIBRARY))
        self.assertEqual(captured["R_LIBS_USER"], str(launcher.R_LIBRARY))


if __name__ == "__main__":
    unittest.main()
