# RouteOpt

---

## Note:

- RouteOpt solves the capacitated vehicle routing problem (CVRP) and the vehicle routing problem with time windows (VRPTW). Please cite our paper [_Two-Stage Learning to Branch in Branch-Price-and-Cut
Algorithms for Solving Vehicle Routing Problems Exactly_](https://www.researchgate.net/publication/374553305_Two-Stage_Learning_to_Branch_in_Branch-Price-and-Cut_Algorithms_for_Solving_Vehicle_Routing_Problems_Exactly) if you use RouteOpt for your research.

- RouteOpt for *academic* use only.

1. The solver consists of the following modules:
    1. **Default Setting (Setting I)**: Full-fledged version with all functionalities enabled, offering the best overall performance.
    2. **Setting II**: For comparing different branching strategies.
    3. **Setting for open instances**: Requires access to supercomputing resources. Detailed explanations of these features can be found in our [paper](https://www.researchgate.net/publication/374553305_Two-Stage_Learning_to_Branch_in_Branch-Price-and-Cut_Algorithms_for_Solving_Vehicle_Routing_Problems_Exactly). For specific usage, please refer to the [User Manual](https://github.com/Zhengzhong-You/RouteOpt-usermanual).

2. RouteOpt employs a novel two-stage learning-to-branch module. The current parameter configuration is suitable for most instances, except for ones where pricing is extremely challenging (e.g., CVRP instances with long routes). For such instances, the parameters need to be adjusted according to our recommendations in the [User Manual](https://github.com/Zhengzhong-You/RouteOpt-usermanual).

3. The VRPTW module of RouteOpt is currently optimized for solving type-2 instances (ones with very loose capacity constraints) fast, and is not suitable for type-1 instances. This is due to the adaptive strategies mentioned in [Bucket Graph](https://pubsonline.informs.org/doi/abs/10.1287/trsc.2020.0985).  

4. The workflow design, computational techniques, and output style of RouteOpt pay a high tribute
   to [VRPSolver](https://vrpsolver.math.u-bordeaux.fr/). We are deeply grateful for Drs. Artur Pessoa, Ruslan Sadykov, Eduardo Uchoa, François Vanderbeck, and others who have contributed to the VRPSolver. Special thanks go to Dr. Ruslan Sadykov for providing us access to the VRPSolver and helping us tune it. Without the VRPSolver, RouteOpt would not have been possible.

## Requirements

---

- [CMake](https://cmake.org/download/) version 3.16 or higher (required).
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) version 3.4.0 or higher (required).
- [XGBoost](https://xgboost.readthedocs.io/en/latest/build.html) version 1.4.2 or higher (required).
- [GUROBI](https://www.gurobi.com/downloads/gurobi-software/) version 10.0 or higher (recommended).
- [CVRPSEP](https://econ.au.dk/research/researcher-websites/jens-lysgaard/cvrpsep/)

Please ensure these requirements are met before proceeding with the following steps.

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

On Linux, use `vim ~/.bashrc`. On MacOS, use `vim ~/.bash_profile` to open the corresponding file. Then add the following line to it.
```
export GUROBI_HOME=<gurobi_root>/<gurobi_version>/<OS>
```

Example: `export GUROBI_HOME=${HOME}/gurobi1000/linux64`

After that, use `source ~/.bashrc`
or `source ~/.bash_profile` to reload the update file

**Step 6:** Revise `FindGUROBI.cmake`:

It's located in `<path to solver>/pub-RouteOpt/Application/CVRP/CVRP/package`

```
find_library(<package>_LIBRARY
        NAMES <lib>
        PATHS "$ENV{<package>_HOME}/lib"
        )
```

Replace the `<lib>` with the right `lib`.

For example, if you use `gurobi1000`, the lib will be `libgurobi100.so` on Linux and `libgurobi100.dylib` on MacOS.

## Useage

---

The library accepts two kinds of parameters: (RouteOpt_VRPTW for VRPTW)

1. Use the `idx/<ins_file.ins>` (Recommended)

   ```
   ./CVRP -d idx/<ins_file.ins> -n -u
   ```

   In this command, `-d` is followed by `.ins` file, and `-n` gets the instance name in the `n`-th line of
   the file. Note that `-u` is an optional parameter to provide an initial upper bound (UB) if a valid one is available. It should directly follow `-u`  without a space.

   For instance, you can type in:

   ```
   ./CVRP -d idx/OldIns.ins -n 33 -u27592
   ```

   This command will read the instance from the `33`-rd line of the `OldIns.ins` file.

   Alternatively, you can just enter:

   ```
   ./CVRP -d idx/OldIns.ins -n 33 
   ```

   In this case, if the `.ins` file contains a UB, it will be automatically read from the file. If no UB is provided, it will be initialized to 1e9.

2. Specify the complete path

   ```
   ./CVRP ./../../../DataForCVRP/OldIns/X-n101-k25.vrp -u27592
   ```

## Additional Information

For details on adjusting:
- Log files/outputs,
- Parameters,
- Branching strategies,

please refer to the [User Manual](https://github.com/Zhengzhong-You/RouteOpt-usermanual).
