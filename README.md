# RouteOpt

---

## Note:

- This is a demo of the CVRP (Capacitated Vehicle Routing Problem) module. The machine learning module is not yet enabled.

- The VRPTW (Vehicle Routing Problem with Time Windows) module will be launched soon. Stay tuned!

- This library is intended *for academic use only*.

## Requirements

---

- [CMake](https://cmake.org/download/) version 3.16 or higher.
- [GUROBI](https://www.gurobi.com/downloads/gurobi-software/) - please follow the link to install GUROBI.

Please ensure these requirements are met before proceeding with the usage of the RouteOpt library.

## Link RouteOpt Library

---

### 1. Linking the Library in Your Project

### Using CMake

If you're using CMake for your project, you can easily link to the RouteOpt library. First, in your CMakeLists.txt, add the following line:

```cmake
find_library(ROUTEOPT_LIBRARY RouteOpt PATHS /path/to/your/library)
```

Replace "/path/to/your/library" with the actual path where your built RouteOpt library is located.

Then, when defining your executable or another library that depends on RouteOpt, you need to link to it:

```cmake
add_executable(your_executable ${YOUR_SOURCES})
target_link_libraries(your_executable ${ROUTEOPT_LIBRARY})
```

Replace "your_executable" with the name of your executable and "${YOUR_SOURCES}" with the variable that contains your source files.

### Without CMake

If you're not using CMake, you would typically pass the path to the library directly to the compiler/linker.

For example, if you're using the g++ compiler, you can do something like this:

```bash
g++ your_source_files.cpp -L/path/to/your/library -lRouteOpt -o your_executable
```

Again, replace "/path/to/your/library" with the actual path where your built RouteOpt library is located, "your_source_files.cpp" with your source files, and "your_executable" with the desired name of your executable.

### 2. Including Headers

In your source files, you need to include the headers of the RouteOpt library. You can do this with the following line:

```cpp
#include <RouteOpt/YourHeader.h>
```

Replace "YourHeader.h" with the actual name of the header you need.

---

Please refer to the specific documentation and examples of the RouteOpt library for more details on how to use its functions and classes in your code. If you encounter any problems or need further assistance, please submit an issue in the RouteOpt repository.

## Useage

---

The library accepts two kinds of parameters:

1. Utilize the `idx/<ins_file.ins>` (Recommended)
   
   ```
   ./your_executable -d idx/<ins_file.ins> -n -u
   ```
   
   In this command, `-d` denotes the `.ins` file and `-n` stands for the instance corresponding to the `n`-th line of the file. Note that `-u` is an optional parameter indicating the initial upper bound, and it should be written without a space following `-u`.
   
   For instance, you can enter:
   
   ```
   ./your_executable -d idx/OldIns.ins -n 33 -u27592
   ```
   
   This command will read the instance from the `33`-rd line of the `OldIns.ins` file.
   
   Alternatively, you can just enter:
   
   ```
   ./your_executable -d idx/OldIns.ins -n 33 
   ```
   
   In this case, if the `.ins` file contains an Upper Bound (UB), it will be automatically read from the file. If no UB is provided, it will be initialized as 1e9.

2. Specify the complete path
   
   ```
   ./your_executable ./../../../DataForCVRP/OldIns/X-n101-k25.vrp -u27592
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
