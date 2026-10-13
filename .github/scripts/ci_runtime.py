"""Package a platform build and run it without rebuilding dependencies."""

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time


def package(workspace, gurobi_root, archive):
    app = workspace / "packages/application/cvrp"
    with tempfile.TemporaryDirectory(dir=archive.parent) as temporary:
        runtime = Path(temporary) / "runtime"
        (runtime / "bin").mkdir(parents=True)
        (runtime / "lib").mkdir()
        (runtime / "instance").mkdir()
        binary = "cvrp.exe" if os.name == "nt" else "cvrp"
        shutil.copy2(app / "bin" / binary, runtime / "bin" / binary)
        shutil.copytree(app / "model", runtime / "model")
        shutil.copy2(app / "instance/P-n20-k2.vrp", runtime / "instance")
        library_dirs = [workspace / "packages/external/xgb/lib", app / "lib/hgs/lib",
                        gurobi_root / "lib", gurobi_root / "bin", app / "bin"]
        if os.name == "nt":
            library_dirs.append(Path("C:/vcpkg/installed/x64-windows/bin"))
        for directory in library_dirs:
            for path in directory.iterdir():
                if path.is_file() and re.search(r"\.(dll|dylib|so(?:\.\d+)*)$", path.name) and "_light" not in path.name:
                    shutil.copy2(path, runtime / "lib" / path.name)
        with tarfile.open(archive, "w:gz") as output:
            output.add(runtime, arcname="runtime")


def run_smoke(runtime, license_file):
    env = dict(os.environ, GRB_LICENSE_FILE=str(license_file), OMP_NUM_THREADS="1")
    library_dir = str(runtime / "lib")
    for key in ("PATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        env[key] = library_dir + os.pathsep + env.get(key, "")
    binary = "cvrp.exe" if os.name == "nt" else "cvrp"
    command = [str(runtime / "bin" / binary), "instance/P-n20-k2.vrp", "-u", "217"]
    for attempt in range(2):
        result = subprocess.run(command, cwd=runtime, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        print(result.stdout, flush=True)
        if result.returncode == 0:
            bounds = dict(re.findall(r"<(LB|UB)=\s*([\d.]+)>", result.stdout))
            if set(bounds) != {"LB", "UB"} or any(float(value) != 216 for value in bounds.values()):
                raise RuntimeError("CVRP smoke test did not prove LB = UB = 216")
            return
        # Finished processes can retain WLS tokens for five minutes. Retry only environment startup.
        if (attempt == 0 and "error in solver.loadEnv(nullptr)" in result.stdout
                and re.search(r"ERROR CODE = 10030\b", result.stdout)):
            print("Gurobi rejected environment startup (10030); waiting 310 seconds before one retry.", flush=True)
            time.sleep(310)
            continue
        raise subprocess.CalledProcessError(result.returncode, command)


def main():
    temporary = Path(os.environ["RUNNER_TEMP"])
    archive = temporary / "runtime.tar.gz"
    if sys.argv[1] == "package":
        package(Path.cwd(), Path(os.environ["GUROBI_ROOT"]), archive)
    elif sys.argv[1] == "smoke":
        with tarfile.open(archive) as source:
            source.extractall(temporary, filter="data")
        keys = ("WLSACCESSID", "WLSSECRET", "LICENSEID")
        values = {key: os.environ.get("GRB_" + key, "") for key in keys}
        if not all(values.values()):
            raise RuntimeError("Missing Gurobi WLS secrets")
        # The license is created only in the test job, outside all caches and artifacts.
        with tempfile.TemporaryDirectory(dir=temporary) as directory:
            license_file = Path(directory) / "gurobi.lic"
            license_file.write_text("".join(f"{key}={value}\n" for key, value in values.items())
                                    + "WLSTOKENDURATION=5\n", encoding="ascii")
            license_file.chmod(0o600)
            run_smoke(temporary / "runtime", license_file)
    else:
        raise SystemExit("Use package or smoke")


if __name__ == "__main__":
    main()
