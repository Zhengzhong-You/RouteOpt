# RouteOpt

---

## Note:

- This solver is designed for solving the CVRP (Capacitated Vehicle Routing Problem) and the VRPTW (Vehicle Routing
  Problem with Time Windows).

- This is intended *for academic use only*.

1. The solver consists of the following modules:
    1. **Default Setting (Setting I)**: All mentioned features are enabled, offering the highest computational performance.
    2. **Setting II**: Used exclusively for comparing different branching strategies.
    3. **Setting for open instances**: Requires the assistance of supercomputing. Detailed explanations of these
       features can be found
       in [Two-Stage Learning](https://www.researchgate.net/publication/374553305_Two-Stage_Learning_to_Branch_in_Branch-Price-and-Cut_Algorithms_for_Solving_Vehicle_Routing_Problems_Exactly).
       For specific usage, please refer to the [User Manual](https://github.com/Zhengzhong-You/RouteOpt-usermanual).

2. Currently, the solver is configured to assemble with the two-stage learning to branch module for CVRP and VRPTW. This
   solver's current parameter configuration is suitable for most instances, except for those instances where pricing is
   extremely challenging (for example, CVRP instances with long routes). For solving instances with long routes, the
   parameters that need to be adjusted can be found in
   the [User Manual](https://github.com/Zhengzhong-You/RouteOpt-usermanual).

3. The VRPTW module of RouteOpt is currently optimized only for the quick solution of type-2 instances, meaning those
   with very loose capacity constraints. This is due to the adaptive strategies mentioned
   in [Bucket Graph](https://pubsonline.informs.org/doi/abs/10.1287/trsc.2020.0985). Please refrain from using RouteOpt
   for solving type-1 instances, as the performance will be subpar! If you must solve the type-1 instances, please turn
   to the [VRPSolver](https://vrpsolver.math.u-bordeaux.fr/) with configuration one.

4. The design workflow, techniques and even outputs of RouteOpt pay a high tribute
   to [VRPSolver](https://vrpsolver.math.u-bordeaux.fr/). The scientific workflow settings of VRPSolver, combined with
   its ingenious computational techniques, taught me invaluable knowledge. I am deeply grateful for contributions to the
   VRPSolver by scholars like Artur Pessoa, Ruslan Sadykov, Eduardo Uchoa, and François Vanderbeck and others. Special
   thanks to Dr. Ruslan Sadykov for generously granting the rights to access the VRPSolver and assisting me with tuning
   VRPSolver to find research directions. Without the VRPSolver, RouteOpt would not have been possible.

## Requirements

---

- [CMake](https://cmake.org/download/) version 3.16 or higher (required).
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) version 3.4.0 or higher (required).
- [XGBoost](https://xgboost.readthedocs.io/en/latest/build.html) version 1.4.2 or higher (required).
- [GUROBI](https://www.gurobi.com/downloads/gurobi-software/) version 10.0 or higher (recommended).
- [CVRPSEP](https://econ.au.dk/research/researcher-websites/jens-lysgaard/cvrpsep/)

Please ensure these requirements are met before proceeding with the usage of the RouteOpt.

## Link Depencices

---

**Step 1:** Clone the RouteOpt repository:

```
https://github.com/Zhengzhong-You/pub-RouteOpt.git
```

**Step 2:** Construct the project directory tree:

```
cd pub-RouteOpt && mkdir Zips Dependency && cd Dependency
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

In Linux, `vim ~/.bashrc`. In MacOS, `vim ~/.bash_profile`. After adding the following command line, `source ~/.bashrc`
or `source ~/.bash_profile`

```
export GUROBI_HOME=<gurobi_root>/<gurobi_version>/<OS>
```

Example: `export GUROBI_HOME=${HOME}/gurobi1000/linux64`

**Step 6:** Revise `FindGUROBI.cmake`:

It's in `<path to solver>/pub-RouteOpt/Application/CVRP/CVRP/package`

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

The library accepts two kinds of parameters: (RouteOpt_VRPTW for VRPTW)

1. Utilize the `idx/<ins_file.ins>` (Recommended)

   ```
   ./CVRP -d idx/<ins_file.ins> -n -u
   ```

   In this command, `-d` denotes the `.ins` file and `-n` stands for the instance corresponding to the `n`-th line of
   the file. Note that `-u` is an optional parameter indicating the initial upper bound, and it should be written
   without a space following `-u`.

   For instance, you can enter:

   ```
   ./CVRP -d idx/OldIns.ins -n 33 -u27592
   ```

   This command will read the instance from the `33`-rd line of the `OldIns.ins` file.

   Alternatively, you can just enter:

   ```
   ./CVRP -d idx/OldIns.ins -n 33 
   ```

   In this case, if the `.ins` file contains an Upper Bound (UB), it will be automatically read from the file. If no UB
   is provided, it will be initialized as 1e9.

2. Specify the complete path

   ```
   ./CVRP ./../../../DataForCVRP/OldIns/X-n101-k25.vrp -u27592
   ```

## Additional Information

For details on:
- Log outputs
- Parameter adjustments
- Toggling techniques,
- Branching strategies

please refer to the [User Manual](https://github.com/Zhengzhong-You/RouteOpt-usermanual).
