#!/usr/bin/env python3
import os
import sys
import glob
import re
import subprocess
import multiprocessing
import shutil


def run_cmd(cmd, cwd=None):
    proc = subprocess.run(
        cmd, shell=True, cwd=cwd,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if proc.returncode != 0:
        print(f"Error in: {cmd}\n{proc.stderr}")
        sys.exit(1)
    return proc.stdout


def cpu_jobs():
    try:
        n = multiprocessing.cpu_count()
        return n if n and n > 0 else 4
    except Exception:
        return 4


def update_cmake(gurobi_path):
    cmake_file = os.path.join("packages", "external", "cmake_modules", "FindGUROBI.cmake")
    if not os.path.exists(cmake_file):
        print(f"Missing {cmake_file}")
        sys.exit(1)

    with open(cmake_file, "r") as f:
        content = f.read()

    new_root_line = f'set(GUROBI_ROOT "{gurobi_path}")'
    content = re.sub(r'set\s*\(\s*GUROBI_ROOT\s*".*?"\s*\)', new_root_line, content)

    lib_dir = os.path.join(gurobi_path, "lib")
    if not os.path.exists(lib_dir):
        print(f"Library dir not found: {lib_dir}")
        sys.exit(1)

    libs = glob.glob(os.path.join(lib_dir, "libgurobi*.so"))
    libs = [lib for lib in libs if "_light" not in os.path.basename(lib)]
    if not libs:
        print(f"No suitable libgurobi*.so found in {lib_dir}")
        sys.exit(1)
    libs.sort()
    lib_file = os.path.basename(libs[-1])

    content = re.sub(
        r'(find_library\s*\(\s*GUROBI_LIBRARY\s*\n\s*NAMES\s+)[^\n]+',
        r'\1' + lib_file,
        content,
        flags=re.MULTILINE,
    )
    with open(cmake_file, "w") as f:
        f.write(content)
    print("Updated FindGUROBI.cmake.")


def ensure_xgboost_version(ext_dir, version_tag):
    xgb_dir = os.path.join(ext_dir, "xgb")
    repo_url = "https://github.com/dmlc/xgboost"

    if not os.path.exists(xgb_dir):
        run_cmd(f'git clone --recursive --branch {version_tag} {repo_url} xgb', cwd=ext_dir)
        run_cmd('git submodule update --init --recursive', cwd=xgb_dir)
        return xgb_dir

    if not os.path.exists(os.path.join(xgb_dir, ".git")):
        run_cmd('rm -rf xgb', cwd=ext_dir)
        run_cmd(f'git clone --recursive --branch {version_tag} {repo_url} xgb', cwd=ext_dir)
        run_cmd('git submodule update --init --recursive', cwd=xgb_dir)
        return xgb_dir

    run_cmd('git fetch --tags --prune', cwd=xgb_dir)
    run_cmd(f'git checkout {version_tag}', cwd=xgb_dir)
    run_cmd('git submodule sync --recursive', cwd=xgb_dir)
    run_cmd('git submodule update --init --recursive', cwd=xgb_dir)
    return xgb_dir


def choose_generator():
    # Prefer Ninja if available; otherwise use Unix Makefiles
    if shutil.which("ninja"):
        return 'Ninja'
    return 'Unix Makefiles'


def cmake_configure(src_dir, build_dir, extra_flags=""):
    os.makedirs(build_dir, exist_ok=True)
    gen = choose_generator()
    # CMake 4.x can fail on projects that set an old minimum; raise the policy floor.
    policy_flag = "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
    cmd = f'cmake -B "{build_dir}" -S "{src_dir}" -DCMAKE_BUILD_TYPE=RelWithDebInfo {policy_flag} -G "{gen}" {extra_flags}'
    run_cmd(cmd)


def cmake_build(build_dir, target=None, jobs=None):
    j = jobs or cpu_jobs()
    tgt = f' --target {target}' if target else ''
    run_cmd(f'cmake --build "{build_dir}"{tgt} -j {j}')


def build_xgboost(xgb_dir):
    # Clean and configure
    run_cmd('rm -rf build bin', cwd=xgb_dir)
    cmake_configure(xgb_dir, os.path.join(xgb_dir, "build"))
    # Build (generator-agnostic)
    cmake_build(os.path.join(xgb_dir, "build"))


def build_hgs():
    hgs_dir = os.path.join("packages", "application", "cvrp", "lib", "hgs")
    run_cmd("rm -rf lib build", cwd=hgs_dir)
    os.makedirs(os.path.join(hgs_dir, "build"), exist_ok=True)
    run_cmd('cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..', cwd=os.path.join(hgs_dir, "build"))
    # Use cmake --build to be generator-agnostic; install target where available
    cmake_build(os.path.join(hgs_dir, "build"))
    cmake_build(os.path.join(hgs_dir, "build"), target="install")


def build_cvrpsep():
    cvrpsep_dir = os.path.join("packages", "external", "cvrpsep")
    run_cmd("rm -rf dep obj", cwd=cvrpsep_dir)
    run_cmd("make -j {}".format(cpu_jobs()), cwd=cvrpsep_dir)


def main():
    gurobi_path = input("Enter Gurobi installation path: ").strip()
    if not os.path.exists(gurobi_path):
        print("Gurobi not found. Exiting.")
        sys.exit(1)
    update_cmake(gurobi_path)

    XGB_VERSION = "v2.0.0"
    ext_dir = os.path.join("packages", "external")
    os.makedirs(ext_dir, exist_ok=True)

    xgb_dir = ensure_xgboost_version(ext_dir, XGB_VERSION)
    build_xgboost(xgb_dir)
    build_hgs()
    build_cvrpsep()

    print("Build process completed successfully.")


if __name__ == '__main__':
    main()
