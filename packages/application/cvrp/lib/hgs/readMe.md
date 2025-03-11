
# About This Fork

This repository is my personal fork of the project, "featuring a modern implementation of the hybrid genetic search (HGS) algorithm, specifically tailored for the capacitated vehicle routing problem (CVRP). In addition to the standard features, this codebase includes an extra neighborhood, known as SWAP*."
For more details on the original project, visit [HGS-CVRP by vidalt](https://github.com/vidalt/HGS-CVRP).

To seamlessly integrate this repository with RouteOpt, I have designed a custom interface.

## Installation Instructions

Follow these steps to set up the project and generate the `libhgscvrp.so` file:

1. **Navigate to the HGS Directory**
   ```bash
   cd hgs
   ```

2. **Create a Library Directory**
   ```bash
   mkdir lib
   ```

3. **Generate Build Configuration**
   ```bash
   cmake .
   ```

4. **Build and Install**
   ```bash
   make install
   ```

After completing these steps, you will find the compiled library file named `libhgscvrp.so` in the `lib` directory.