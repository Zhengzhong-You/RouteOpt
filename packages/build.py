#!/usr/bin/env python3
import os
import sys
import glob
import re
import subprocess

def run_cmd(cmd, cwd=None):
    proc = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        print(f"Error in: {cmd}\n{proc.stderr}")
        sys.exit(1)
    return proc.stdout

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
    if not libs:
        print(f"No libgurobi*.so found in {lib_dir}")
        sys.exit(1)
    libs.sort()
    lib_file = os.path.basename(libs[-1])
    content = re.sub(r'(find_library\s*\(\s*GUROBI_LIBRARY\s*\n\s*NAMES\s+)[^\n]+', r'\1' + lib_file, content, flags=re.MULTILINE)
    with open(cmake_file, "w") as f:
        f.write(content)
    print("Updated FindGUROBI.cmake.")

def main():
    gurobi_path = input("Enter Gurobi installation path: ").strip()
    if not os.path.exists(gurobi_path):
        print("Gurobi not found. Exiting.")
        sys.exit(1)
    update_cmake(gurobi_path)

    ext_dir = os.path.join("packages", "external")
    xgb_dir = os.path.join(ext_dir, "xgb")
    if not os.path.exists(xgb_dir):
        run_cmd("git clone --recursive https://github.com/dmlc/xgboost", cwd=ext_dir)
    run_cmd("rm -rf build bin", cwd=xgb_dir)
    run_cmd('cmake -B build -S . -DCMAKE_BUILD_TYPE=RelWithDebInfo -GNinja', cwd=xgb_dir)
    run_cmd("ninja -j $(nproc)", cwd=os.path.join(xgb_dir, "build"))

    hgs_dir = os.path.join("packages", "application", "cvrp", "lib", "hgs")
    run_cmd("rm -rf lib build", cwd=hgs_dir)
    run_cmd("mkdir build", cwd=hgs_dir)
    run_cmd("cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..", cwd=os.path.join(hgs_dir, "build"))
    run_cmd("make intall -j $(nproc)", cwd=os.path.join(hgs_dir, "build"))

    cvrpsep_dir = os.path.join("packages", "external", "cvrpsep")
    run_cmd("rm -rf dep obj", cwd=cvrpsep_dir)
    run_cmd("make", cwd=cvrpsep_dir)

    print("Build process completed successfully.")

if __name__ == '__main__':
    main()
