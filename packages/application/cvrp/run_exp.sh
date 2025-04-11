#!/bin/bash
# run_experiments.sh
# This script runs multiple CVRP experiments over a list of instance files.
# For each .vrp instance, it runs three variants:
#   1. Standard variant (using 'cvrp')
#   2. Root cut only variant (using 'cvrp_rootcutonly')
#   3. Robust child variant (using 'cvrp_robustchild')
#
# The script automatically names the output logs using the instance file’s basename.
# For example, if the file is "instance/P-n20-k2.vrp", the corresponding log
# for the standard run is "./logs/output_P-n20-k2.txt" etc.

# Array of VRP instance files.
vrp_files=(
    "instance/P-n20-k2.vrp"
    "instance/P-n60-k10.vrp"
    "instance/P-n65-k10.vrp"
    "instance/P-n60-k15.vrp"
    "instance/P-n70-k10.vrp"
    "instance/P-n76-k5.vrp"
    "instance/X-n101-k25.vrp"
    #"instance/X-n120-k6.vrp"
    "instance/X-n115-k10.vrp"
    "instance/X-n110-k13.vrp"
)

# u_values=(217 745 793 969 828 628 27592 13333 12748 14972)
# u_values=(217 745 793 969 828 628 27592 12748 14972)

# add 1%
u_values=(219 752 801 979 836 634 27868 12875 15122)

# Array of corresponding -u parameters.
# u_values=(13340)

# Total number of runs (3 experiments per instance)
n_instances=${#vrp_files[@]}
total_runs=$(( n_instances * 3 ))

echo "Starting experiments: Total runs = $total_runs"

# Function to display a simple progress bar.
show_progress() {
    current=$1
    percent=$(( (current * 100) / total_runs ))
    bar_length=30
    filled_length=$(( (bar_length * current) / total_runs ))
    bar=$(printf "%${filled_length}s" | tr ' ' '#')
    empty=$(printf "%$((bar_length - filled_length))s")
    echo -ne "\rProgress: [${bar}${empty}] ${percent}% (${current}/${total_runs})"
}

# Initialize run counter.
run_count=0

# Loop over each instance file.
for i in "${!vrp_files[@]}"; do
    file="${vrp_files[$i]}"
    u="${u_values[$i]}"
    # Extract the base name (e.g., "P-n20-k2") without the .vrp extension.
    base=$(basename "$file" .vrp)
    
    # 1. Standard variant using 'cvrp'
    run_count=$((run_count+1))
    show_progress "$run_count"
    ./asaved_exe/cvrp "$file" -u "$u" > "/home/haoran/solver/DataCollect/logs/outputUB1_${base}.txt" 2>&1

    # 2. Root cut only variant using 'cvrp_rootcutonly'
    run_count=$((run_count+1))
    show_progress "$run_count"
    ./asaved_exe/cvrp_rootcutonly "$file" -u "$u" > "/home/haoran/solver/DataCollect/logs/outputUB1_${base}_onlyr.txt" 2>&1

    # 3. Robust child variant using 'cvrp_robustchild'
    run_count=$((run_count+1))
    show_progress "$run_count"
    ./asaved_exe/cvrp_robustchild "$file" -u "$u" > "/home/haoran/solver/DataCollect/logs/outputUB1_${base}_robustc.txt" 2>&1
done

# Final progress update
show_progress "$total_runs"
echo -e "\nAll runs completed."

