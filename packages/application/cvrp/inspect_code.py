# !/usr/bin/env python3
import os
import re
import sys
import subprocess
import argparse


def parse_template(template_path):
    parameters = {}
    with open(template_path, "r") as file:
        lines = file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # Skip empty or comment lines.
        if not line or line.startswith("#"):
            i += 1
            continue
        # Check if the line defines an enum class.
        enum_match = re.match(r'^enum\s+class\s+(\w+)', line)
        if enum_match:
            enum_name = enum_match.group(1)
            # Begin collecting the enum definition
            enum_def = line
            # Continue reading until we find a semicolon, indicating the end of the enum block.
            while ';' not in enum_def and i < len(lines) - 1:
                i += 1
                enum_def += "\n" + lines[i].rstrip("\n")
            # Store the enum definition under its name
            parameters[enum_name] = enum_def
        else:
            # Process regular parameter lines: <variable> <value>
            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                i += 1
                continue
            parameters[parts[0]] = parts[1]
        i += 1
    return parameters


def update_files(base_dir, parameters):
    import os, sys, re
    src_dir = os.path.join(base_dir, "src")
    changed = set()
    for root, _, files in os.walk(src_dir):
        for filename in files:
            if filename.endswith("_macro.hpp"):
                file_path = os.path.join(root, filename)
                with open(file_path, "r") as f:
                    file_content = f.read()
                original_content = file_content
                file_modified = False

                for param_name, new_value in parameters.items():
                    # Check if the new_value is an enum definition by detecting "enum class"
                    if new_value.lstrip().startswith("enum class"):
                        # Build a regex pattern to match the entire enum block by its name.
                        pattern = r"(enum\s+class\s+" + re.escape(param_name) + r"\s*\{.*?\};)"
                        if re.search(pattern, file_content, flags=re.DOTALL):
                            file_content, count = re.subn(pattern, new_value, file_content, flags=re.DOTALL)
                            if count > 0:
                                print(f"Updated enum '{param_name}' in file '{file_path}'")
                                file_modified = True
                                changed.add(param_name)
                            else:
                                print(f"Error: Enum '{param_name}' in file '{file_path}' could not be updated.")
                                sys.exit(1)
                        else:
                            continue
                    else:
                        # For non-enum parameters, use the regex-based replacement.
                        pattern = (
                                r"((?:(?:static\s+)?constexpr)\s+[\w:]+\s+"
                                + re.escape(param_name) +
                                r"\s*\{\s*)([^}]+?)(\s*\}\s*;)"
                        )
                        match = re.search(pattern, file_content)
                        if match:
                            new_line = match.group(1) + new_value + match.group(3)
                            file_content = re.sub(pattern, new_line, file_content, count=1)
                            print(f"Updated parameter '{param_name}' in file '{file_path}'")
                            file_modified = True
                            changed.add(param_name)
                        else:
                            continue

                if file_modified and file_content != original_content:
                    with open(file_path, "w") as f:
                        f.write(file_content)
                    print(f"Updated file: {file_path}")
    for param_name in parameters:
        if param_name not in changed:
            print(f"Warning: Parameter '{param_name}' not found in any file.")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Update application based on template file.")
    parser.add_argument("template_name", help="Template file name, e.g. cvrp_train_model_1.txt")
    args = parser.parse_args()
    template_name = args.template_name

    base_dir = os.getcwd()

    # Step 1: Check if the template exists in the "templates" folder
    template_path = os.path.join(base_dir, "templates", template_name)
    if not os.path.exists(template_path):
        print(f"Error: Template file '{template_name}' not found in the templates folder.")
        sys.exit(1)

    # Step 2: Parse the template file to extract parameters
    # Expected format per line: <variable> <value>
    enum_parameters = parse_template(template_path)
    update_files(base_dir, enum_parameters)

    # Step 3: Update CMakeLists.txt with new project name (template name without extension)
    new_project_name = os.path.splitext(template_name)[0]
    cmake_path = os.path.join(base_dir, "CMakeLists.txt")
    if not os.path.exists(cmake_path):
        print("Error: CMakeLists.txt not found.")
        sys.exit(1)

    with open(cmake_path, "r") as f:
        cmake_lines = f.readlines()

    project_line_found = False
    updated_cmake_lines = []
    for line in cmake_lines:
        if line.strip().lower().startswith("project("):
            project_line_found = True
            line = f"project({new_project_name})\n"
        updated_cmake_lines.append(line)

    if not project_line_found:
        print("Error: No project() line found in CMakeLists.txt.")
        sys.exit(1)

    with open(cmake_path, "w") as f:
        f.writelines(updated_cmake_lines)
    print(f"CMakeLists.txt updated with project name: {new_project_name}")

    print("Compiling code by running 'sh build.sh'...")
    compile_proc = subprocess.run(["sh", "build.sh"], capture_output=True, text=True)

    if compile_proc.returncode != 0 or "error" in compile_proc.stderr.lower():
        print("Compilation failed. Output:")
        print(compile_proc.stdout)
        print(compile_proc.stderr)
        sys.exit(1)
    else:
        print("Compilation succeeded.")


if __name__ == "__main__":
    main()
