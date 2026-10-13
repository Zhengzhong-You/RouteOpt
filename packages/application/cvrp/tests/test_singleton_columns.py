"""Column lifecycle regressions (#23 and #25); requires a built cvrp and solver license.

Run: python3 tests/test_singleton_columns.py
Set ROUTEOPT_TEST_BINARY to test a different build.
Fixtures: https://github.com/Zhengzhong-You/RouteOpt/issues/23
          https://github.com/Zhengzhong-You/RouteOpt/pull/25
"""

import math
import os
from pathlib import Path
import re
import subprocess
import unittest


APP = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("ROUTEOPT_TEST_BINARY", APP / "bin/cvrp")).resolve()


class SingletonColumnsTest(unittest.TestCase):
    def check_solution(self, instance, expected, *arguments, timeout=180):
        run = subprocess.run(
            [str(BINARY), str(instance), *arguments],
            cwd=APP, env={**os.environ, "OMP_NUM_THREADS": "1"},
            capture_output=True, text=True, timeout=timeout,
        )
        self.assertEqual(run.returncode, 0, run.stdout[-3000:] + run.stderr)
        output = run.stdout
        for bound in ("LB", "UB"):
            match = re.search(rf"<{bound}=\s*([\d.]+)>", output)
            self.assertIsNotNone(match, output[-3000:])
            self.assertEqual(float(match[1]), expected, output[-3000:])
        routes = [list(map(int, line.split("-"))) for line in
                  output.split("<Solution>\n")[-1].split("<UB=")[0].splitlines()]

        coordinates, demands = {}, {}
        section = None
        for line in instance.read_text().splitlines():
            if line.startswith("CAPACITY"):
                capacity = int(line.split(":")[1])
            elif line in ("NODE_COORD_SECTION", "DEMAND_SECTION", "DEPOT_SECTION"):
                section = line
            elif section == "NODE_COORD_SECTION":
                node, x, y = map(int, line.split())
                coordinates[node - 1] = (x, y)
            elif section == "DEMAND_SECTION":
                node, demand = map(int, line.split())
                demands[node - 1] = demand

        self.assertEqual(sorted(node for route in routes for node in route),
                         sorted(node for node in demands if node != 0))
        total = 0
        for route in routes:
            self.assertTrue(route)
            self.assertLessEqual(sum(demands[node] for node in route), capacity)
            sequence = [0, *route, 0]
            total += sum(math.floor(math.dist(coordinates[a], coordinates[b]) + 0.5)
                         for a, b in zip(sequence, sequence[1:]))
        self.assertEqual(total, expected)
        return output, routes

    def test_issue23(self):
        instance = APP / "tests/data/XML50_1261_24.vrp"
        for arguments in (("-u", "15021"), ()):
            with self.subTest(arguments=arguments):
                output, routes = self.check_solution(instance, 15020, *arguments)
                self.assertIn("Process Enumeration Matrix", output)
                self.assertTrue(any(len(route) == 1 for route in routes))

    def test_existing_small_instance(self):
        self.check_solution(APP / "instance/P-n20-k2.vrp", 216)

    def test_enumeration_cut_buffers(self):
        self.check_solution(APP / "tests/data/XML20_3253_21.vrp", 4074,
                            "-u", "4075", timeout=900)


if __name__ == "__main__":
    unittest.main()
