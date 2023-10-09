# RouteOpt

---

Note:

- This package is *still under development*.
- This package is *only for academic use*.


## Requirements

---


- [GUROBI](https://www.gurobi.com/downloads/gurobi-software/) 
- [XGBoost](https://xgboost.readthedocs.io/en/stable/build.html)
- [CMake](https://cmake.org/download/) versions 3.16 and higher
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) versions 3.4.0 and higher


## Installation:

---

 Obtain RouteOpt by typing:

```
git clone https://github.com/Zhengzhong-You/RouteOpt.git
```

Construct project tree:

```
cd RouteOpt && mkdir Zips Dependency && cd Dependency
```

Build cvrpsep library:

```
git clone https://github.com/Zhengzhong-You/cvrpsep.git && cd cvrpsep
make
```

Download eigen:

```
cd ../
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
unzip eigen-3.4.0.zip && mv eigen-3.4.0.zip ../Zips
```

Download xgboost (for ML purpose):

If building on MacOS:

Obtain `libomp` from [Homebrew](https://brew.sh/) first.

```
brew install libomp
```

Rest is the same as building on Linux.

```
git clone --recursive https://github.com/dmlc/xgboost
cd xgboost && mkdir build && cd build
cmake ..
make -j$(nproc)
```

Set `GUROBI_HOME` :

In Linux, `vim ~/.bashrc`. In MacOS, `vim ~/.bash_profile`. After adding the following command line, `source ~/.bashrc` or `source ~/.bash_profile`

```
export GUROBI_HOME=<gurobi_root>/<gurobi_version>/<OS>
```

example: `export GUROBI_HOME=${HOME}/gurobi952/linux64`

Revise `FindGUROBI.cmake` and `FindXGB.cmake`:

They are in `<path to solver>/RouteOpt/Application/CVRP/CVRP7/package`

```
find_library(<package>_LIBRARY
        NAMES <lib>
        PATHS "$ENV{<package>_HOME}/lib"
        )
```

Replace the `<lib>` with the right `lib`.

For example, for `FindGUROBI.cmake` and we use `gurobi952`, the lib will be `libgurobi95.so` in Linux and `libgurobi95.dylib` in MacOS; for `FindXGB.cmake` , the lib will be `libxgboost.so` in Linux and `libxgboost.dylib` in MacOS.

Build:

```
cd ../../../Application/CVRP/CVRP7/ && cmake .
make -j$(nproc)
```

## Usage:

---

Two ways of running CVRP Application:

- Take the advantage of `idx/<ins_file.ins>` (Recommend)
  
  ```
  ./CVRP -d idx/<ins_file.ins> -n -u
  ```
  
  where `-d` specifies the `.ins` file and `-n` specifies the instance corresponding to the `n`-th line of the file. Note that `-u` is an optinal parameter indicating the initial upper bound and there is no space after `-u`.
  
  For example, you can type
  
  ```
  ./CVRP -d idx/OldIns.ins -n 33 -u27592
  ```
  
  , which means we read instance from `33`-rd line of `OldIns.ins` file.
  
  Or just type
  
  ```
  ./CVRP -d idx/OldIns.ins -n 33 
  ```
  
  , if the .ins file contains UB, it will read the UB automatically, otherwise, the UB will be initialized as 1e9.

- Specify the full path
  
  ```
  ./CVRP ./../../../DataForCVRP/OldIns/X-n101-k25.vrp -u27592
  ```

## Output:

---

- **Branching**:
  
  - *nd_ind*: index of node
  
  - *nd_col*: # of columns in node's LP
  
  - *nd_val*: the assigned value of node (used by **priority queue** when branching)
  
  - *nd_dep*: the depth of the node in the branch-and-bound tree (BBT) (tree level)
  
  - *nd_rmn*: the remaining nodes of the BBT
  
  - *et/s*: elapsed time so far
  
  - *lb*: global lower bound
  
  - *ub*: global upper bound

- **Cutting**:
  
  - *rcc*: # of rounded capacity cuts gnerated
  
  - *r3*: # of rank1-3 cuts
  
  - *r1*: # of rank1-1 cuts

- **Pricing**:
  
  - *it*: # of iterations
  
  - *chgcol*: # of columns generated in column generation (CG)
  
  - *ncol*: # of columns in LP
  
  - *ncstr*: # of rows in LP
  
  - *mt/s*: time interval of two continus infomation printing for **re-optimizing LP**
  
  - *spt/s*: time interval of two continus infomation printing for **column generation**
  
  - *lpval*: objective value of LP
  
  - *lb*: global lower bound
  
  - *ub*: global upper bound

- **ArcElimination**:
  
  - *half_forwarding/s*: time spent in half_forward labeling
  
  - *concatenate/s*: time spent in concatenation used for the extension of few labels
  
  - *left_half/s:* time spent in left_half_forward labeling
  
  - *BucketArcs*: # of bucket arcs remained
  
  - *prev.*: ratio of # of current bucket arcs compared to these of last ArcElimination
  
  - *max.*: ratio of # of current bucket arcs compared to all bucket arcs
  
  - *JumpArcs*: # of jump arcs generated

- **EnumerateRoutes**:
  
  - *num_labels*: # of labels generated in half_forward labeling
  
  - *num_routes*: # of routes generated in half_forward labeling
  
  - *num_routes_now*: # of routes in total
  
  - *concatenate*: time spent in concatenation
  
  - *all_routes*: # of revised routes in total (used in generating the column pool)
