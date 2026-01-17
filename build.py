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


def is_macos():
    return sys.platform == "darwin"


def is_windows():
    return os.name == "nt"


def rm_rf(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def rm_rf_many(paths, cwd=None):
    base = cwd or os.getcwd()
    for rel in paths:
        rm_rf(os.path.join(base, rel))


def has_gurobi_libs(path, exts):
    lib_dir = os.path.join(path, "lib")
    if not os.path.isdir(lib_dir):
        return False
    patterns = []
    for ext in exts:
        if ext == ".lib":
            patterns.append(f"gurobi[0-9]*{ext}")
            patterns.append(f"libgurobi[0-9]*{ext}")
        else:
            patterns.append(f"libgurobi[0-9]*{ext}")
    for pattern in patterns:
        if glob.glob(os.path.join(lib_dir, pattern)):
            return True
    return False


def is_gurobi_root(path, exts=None):
    if not (
        os.path.isdir(path)
        and os.path.isdir(os.path.join(path, "lib"))
        and os.path.isfile(os.path.join(path, "include", "gurobi_c.h"))
    ):
        return False
    if exts is None:
        return True
    return has_gurobi_libs(path, exts)


def resolve_gurobi_root(gurobi_path):
    if is_macos():
        exts = (".dylib", ".a")
        preferred = ("macos_universal2", "macos_arm64", "macosx_arm64", "mac64")
    elif is_windows():
        exts = (".lib",)
        preferred = ("win64",)
    else:
        exts = (".so", ".a")
        preferred = ("linux64", "linux_arm64", "win64")

    if is_gurobi_root(gurobi_path, exts):
        return gurobi_path

    # Common platform subfolders under a versioned root.
    for sub in preferred:
        candidate = os.path.join(gurobi_path, sub)
        if is_gurobi_root(candidate, exts):
            return candidate

    # Fallback: scan one level deep for a folder that looks like a Gurobi root.
    if os.path.isdir(gurobi_path):
        for entry in os.listdir(gurobi_path):
            candidate = os.path.join(gurobi_path, entry)
            if is_gurobi_root(candidate, exts):
                return candidate

    return None


def find_gurobi_libs(lib_dir):
    patterns = (
        "libgurobi[0-9]*.so",
        "libgurobi[0-9]*.dylib",
        "libgurobi[0-9]*.a",
        "gurobi[0-9]*.lib",
        "libgurobi[0-9]*.lib",
    )
    libs = []
    for pattern in patterns:
        libs.extend(glob.glob(os.path.join(lib_dir, pattern)))
    return [lib for lib in libs if "_light" not in os.path.basename(lib)]


def update_cmake(gurobi_path):
    cmake_file = os.path.join("packages", "external", "cmake_modules", "FindGUROBI.cmake")
    if not os.path.exists(cmake_file):
        print(f"Missing {cmake_file}")
        sys.exit(1)

    with open(cmake_file, "r") as f:
        content = f.read()

    # Use forward slashes for CMake and avoid backslash escapes in regex replacement.
    cmake_root = gurobi_path.replace("\\", "/")
    new_root_line = f'set(GUROBI_ROOT "{cmake_root}")'
    content = re.sub(
        r'set\s*\(\s*GUROBI_ROOT\s*".*?"\s*\)',
        lambda _m: new_root_line,
        content,
    )

    lib_dir = os.path.join(gurobi_path, "lib")
    if not os.path.exists(lib_dir):
        print(f"Library dir not found: {lib_dir}")
        sys.exit(1)

    libs = find_gurobi_libs(lib_dir)
    if not libs:
        print(f"No suitable libgurobi* found in {lib_dir}")
        sys.exit(1)

    # Best-effort update if the FindGUROBI module still uses a fixed library name.
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
    if is_windows():
        if shutil.which("nmake"):
            return 'NMake Makefiles'
        # Fall back to the default VS generator if available
        return 'Visual Studio 17 2022'
    return 'Unix Makefiles'


def cmake_configure(src_dir, build_dir, extra_flags=""):
    os.makedirs(build_dir, exist_ok=True)
    gen = choose_generator()
    # CMake 4.x can fail on projects that set an old minimum; raise the policy floor.
    policy_flag = "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
    toolchain = os.environ.get("CMAKE_TOOLCHAIN_FILE")
    toolchain_flag = f'-DCMAKE_TOOLCHAIN_FILE="{toolchain}"' if toolchain else ""
    cmd = (
        f'cmake -B "{build_dir}" -S "{src_dir}" -DCMAKE_BUILD_TYPE=RelWithDebInfo '
        f'{policy_flag} {toolchain_flag} -G "{gen}" {extra_flags}'
    )
    run_cmd(cmd)


def cmake_build(build_dir, target=None, jobs=None):
    j = jobs or cpu_jobs()
    tgt = f' --target {target}' if target else ''
    run_cmd(f'cmake --build "{build_dir}"{tgt} -j {j}')


def build_xgboost(xgb_dir):
    # Clean and configure
    rm_rf_many(["build", "bin"], cwd=xgb_dir)
    xgb_flags = []
    if is_macos():
        # OpenMP is mainly for XGBoost training; this build does not train models,
        # so we do not require libomp to be installed here.
        xgb_flags.append("-DUSE_OPENMP=OFF")
    cmake_configure(xgb_dir, os.path.join(xgb_dir, "build"), " ".join(xgb_flags))
    # Build (generator-agnostic)
    cmake_build(os.path.join(xgb_dir, "build"))


def build_hgs():
    hgs_dir = os.path.join("packages", "application", "cvrp", "lib", "hgs")
    rm_rf_many(["lib", "build"], cwd=hgs_dir)
    os.makedirs(os.path.join(hgs_dir, "build"), exist_ok=True)
    cmake_configure(hgs_dir, os.path.join(hgs_dir, "build"))
    # Use cmake --build to be generator-agnostic; install target where available
    cmake_build(os.path.join(hgs_dir, "build"))
    cmake_build(os.path.join(hgs_dir, "build"), target="install")


def build_cvrpsep():
    cvrpsep_dir = os.path.join("packages", "external", "cvrpsep")
    if shutil.which("make") and not is_windows():
        run_cmd("rm -rf dep obj", cwd=cvrpsep_dir)
        run_cmd("make -j {}".format(cpu_jobs()), cwd=cvrpsep_dir)
        return

    rm_rf_many(["build", "obj"], cwd=cvrpsep_dir)
    build_dir = os.path.join(cvrpsep_dir, "build")
    cmake_configure(cvrpsep_dir, build_dir)
    cmake_build(build_dir)


def main():
    gurobi_path = input("Enter Gurobi installation path: ").strip()
    if not os.path.exists(gurobi_path):
        print("Gurobi path does not exist. Exiting.")
        sys.exit(1)
    gurobi_root = resolve_gurobi_root(gurobi_path)
    if not gurobi_root:
        print("Gurobi root not found. Please provide the platform folder (e.g., linux64 or macos_universal2).")
        sys.exit(1)
    update_cmake(gurobi_root)

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
