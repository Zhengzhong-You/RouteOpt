import contextlib
import io
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

import ci_runtime


class SmokeTest(unittest.TestCase):
    def run_case(self, outputs):
        with patch.object(ci_runtime.subprocess, "run", side_effect=outputs) as run, \
                patch.object(ci_runtime.time, "sleep") as sleep, \
                contextlib.redirect_stdout(io.StringIO()):
            ci_runtime.run_smoke(Path("/runtime"), Path("/gurobi.lic"))
        return run, sleep

    def test_success(self):
        run, sleep = self.run_case([subprocess.CompletedProcess([], 0, "<LB= 216>\n<UB= 216>\n")])
        self.assertEqual(run.call_count, 1)
        sleep.assert_not_called()

    def test_retry_license_startup_only(self):
        failure = "error in solver.loadEnv(nullptr)\nERROR CODE = 10030"
        run, sleep = self.run_case([subprocess.CompletedProcess([], 1, failure),
                                   subprocess.CompletedProcess([], 0, "<LB= 216>\n<UB= 216>")])
        self.assertEqual(run.call_count, 2)
        sleep.assert_called_once_with(310)

    def test_does_not_retry_other_failures(self):
        for output in ("error in solver.loadEnv(nullptr)\nERROR CODE = 10009",
                       "error in solver.optimize()\nERROR CODE = 10030", "heap corruption"):
            with self.subTest(output=output), \
                    patch.object(ci_runtime.subprocess, "run", return_value=subprocess.CompletedProcess([], 1, output)) as run, \
                    patch.object(ci_runtime.time, "sleep") as sleep, \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(subprocess.CalledProcessError):
                    ci_runtime.run_smoke(Path("/runtime"), Path("/gurobi.lic"))
                self.assertEqual(run.call_count, 1)
                sleep.assert_not_called()

    def test_retry_is_bounded(self):
        failure = subprocess.CompletedProcess([], 1, "error in solver.loadEnv(nullptr)\nERROR CODE = 10030")
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_case([failure, failure])

    def test_wrong_objective_is_failure(self):
        with self.assertRaises(RuntimeError):
            self.run_case([subprocess.CompletedProcess([], 0, "<LB= 217>\n<UB= 217>")])


if __name__ == "__main__":
    unittest.main()
