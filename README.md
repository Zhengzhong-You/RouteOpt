# RouteOpt

---

## Note:

- This is a demo of the CVRP (Capacitated Vehicle Routing Problem) module. The machine learning module is not yet enabled.

- The VRPTW (Vehicle Routing Problem with Time Windows) module will be launched soon. Stay tuned!

- This executable (compiled using Linux) is intended *for academic use only*.

## Requirements

---

- [CMake](https://cmake.org/download/) version 3.16 or higher.
- [GUROBI](https://www.gurobi.com/downloads/gurobi-software/) version 10.0
- [CVRPSEP](https://econ.au.dk/research/researcher-websites/jens-lysgaard/cvrpsep/)
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) version 3.4.0

Please ensure these requirements are met before proceeding with the usage of the RouteOpt.

## Link Depencices

---

**Step 1:** Clone the RouteOpt repository:

```
https://github.com/Zhengzhong-You/RouteOpt-CVRP.git
```

**Step 2:** Construct the project directory tree:

```
cd RouteOpt-CVRP && mkdir Zips Dependency && cd Dependency
```

**Step 3:** Build the cvrpsep library:

```
(obtain the code) && cd cvrpsep
make
```

**Step 4:** Download Eigen:

```
cd ../
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
unzip eigen-3.4.0.zip && mv eigen-3.4.0.zip ../Zips
```

**Step 5:** Set the `GUROBI_HOME` environment variable:

In Linux, `vim ~/.bashrc`. In MacOS, `vim ~/.bash_profile`. After adding the following command line, `source ~/.bashrc` or `source ~/.bash_profile`

```
export GUROBI_HOME=<gurobi_root>/<gurobi_version>/<OS>
```

Example: `export GUROBI_HOME=${HOME}/gurobi1000/linux64`

**Step 6:** Revise `FindGUROBI.cmake`:

It's in `<path to solver>/RouteOpt/Application/CVRP/CVRP7/package`

```
find_library(<package>_LIBRARY
        NAMES <lib>
        PATHS "$ENV{<package>_HOME}/lib"
        )
```

Replace the `<lib>` with the right `lib`.

For example, if we use `gurobi1000`, the lib will be `libgurobi100.so` in Linux and `libgurobi100.dylib` in MacOS.

## Useage

---

The library accepts two kinds of parameters:

1. Utilize the `idx/<ins_file.ins>` (Recommended)
   
   ```
   ./RouteOpt_CVRP-d idx/<ins_file.ins> -n -u
   ```
   
   In this command, `-d` denotes the `.ins` file and `-n` stands for the instance corresponding to the `n`-th line of the file. Note that `-u` is an optional parameter indicating the initial upper bound, and it should be written without a space following `-u`.
   
   For instance, you can enter:
   
   ```
   ./RouteOpt_CVRP -d idx/OldIns.ins -n 33 -u27592
   ```
   
   This command will read the instance from the `33`-rd line of the `OldIns.ins` file.
   
   Alternatively, you can just enter:
   
   ```
   ./RouteOpt_CVRP -d idx/OldIns.ins -n 33 
   ```
   
   In this case, if the `.ins` file contains an Upper Bound (UB), it will be automatically read from the file. If no UB is provided, it will be initialized as 1e9.

2. Specify the complete path
   
   ```
   ./RouteOpt_CVRP ./../../../DataForCVRP/OldIns/X-n101-k25.vrp -u27592
   ```

## Output:

---

- **Branching**:
  
  - *nd_ind*: Node index
  
  - *nd_col*: Number of columns in the node's Linear Programming (LP)
  
  - *nd_val*: The value assigned to the node (used by the **priority queue** during branching)
  
  - *nd_dep*: The depth of the node in the branch-and-bound tree (BBT) (tree level)
  
  - *nd_rmn*: The remaining nodes in the BBT
  
  - *et/s*: Elapsed time up until the current point
  
  - *lb*: Global lower bound
  
  - *ub*: Global upper bound

- **Cutting**:
  
  - *rcc*: Number of rounded capacity cuts generated
  
  - *r3*: Number of rank1-3 cuts
  
  - *r1*: Number of rank1-1 cuts

- **Pricing**:
  
  - *it*: Number of iterations
  
  - *chgcol*: Number of columns generated in Column Generation (CG)
  
  - *ncol*: Number of columns in LP
  
  - *ncstr*: Number of rows in LP
  
  - *mt/s*: Time interval between two consecutive information printouts during **LP re-optimization**
  
  - *spt/s*: Time interval between two consecutive information printouts during **column generation**
  
  - *lpval*: Objective value of LP
  
  - *lb*: Global lower bound
  
  - *ub*: Global upper bound

- **ArcElimination**:
  
  - *half_forwarding/s*: Time spent in half_forward labeling
  
  - *concatenate/s*: Time spent in concatenation, used for the extension of a few labels
  
  - *left_half/s:* Time spent in left_half_forward labeling
  
  - *BucketArcs*: Number of remaining bucket arcs
  
  - *prev.*: Ratio of current bucket arcs number to the number from the last ArcElimination
  
  - *max.*: Ratio of current bucket arcs number to the total bucket arcs
  
  - *JumpArcs*: Number of jump arcs generated

- **EnumerateRoutes**:
  
  - *num_labels*: Number of labels generated in half_forward labeling
  
  - *num_routes*: Number of routes generated in half_forward labeling
  
  - *num_routes_now*: Total number of routes
  
  - *concatenate*: Time spent in concatenation
  
  - *all_routes*: Total number of revised routes (used in generating the column pool)
