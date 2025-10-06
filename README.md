# RouteOpt 2.0

Welcome to **RouteOpt 2.0** – the latest, state-of-the-art optimization framework that is now free for use! We are excited to introduce this new version, designed to deliver enhanced performance, flexibility, and ease-of-use for solving complex routing and optimization problems.

### Important Notice

If you have any questions about using this repository — not just bug reports — feel free to open an issue. We typically respond within one day. However, before submitting, please check whether a similar issue has already been raised or resolved.

Additionally, using AI-assisted coding tools such as **Claude Code** can significantly improve your understanding and navigation of this repository. Even I use such tools to facilitate code exploration (this codebase is quite extensive, and my memory isn’t perfect sometimes!).  

Have fun exploring **RouteOpt 2.0**!

## Introduction

RouteOpt 2.0 has been meticulously developed to address a wide range of optimization challenges, including the Capacitated Vehicle Routing Problem (CVRP) and the Vehicle Routing Problem with Time Windows (VRPTW). With an improved branching module, modular architecture, and advanced solver integration, RouteOpt 2.0 empowers you to build custom solvers with unparalleled efficiency.

### Key Features

- **Enhanced Branching Module:** Leverage sophisticated branching strategies and unified interfaces that seamlessly integrate multiple optimization solvers.
- **Modular Architecture:** Enjoy a flexible framework where modules are easily extended or replaced to suit your specific problem requirements.
- **Detailed Documentation:** Access comprehensive and powerful documentation to guide you through every aspect of RouteOpt 2.0.
- **Improved Performance:** Experience significant performance improvements and scalability for large, complex problem instances.

## Getting Started

To get started with RouteOpt 2.0, clone the repository and follow our installation instructions provided in the documentation. Whether you are a researcher or a practitioner, our guides and examples will help you quickly harness the full potential of RouteOpt 2.0.


## Updates

## 1) Update on XGBoost Compatibility

The latest XGBoost (v3.0.4) does not support models built with older versions such as v2.0.0. To avoid compatibility issues, please reinstall XGBoost. The simplest way to rebuild is to remove the old folder and run the build script again:

```bash
cd /RouteOpt/packages/external
rm -rf xgb
cd /RouteOpt
python3 build.py
```

For new users, simply run `build.py`. The script will automatically pull XGBoost v2.0.0 from the official repository. We apologize for any inconvenience.

## 2) RouteOpt 2.0 Parameter Management

All key parameters for RouteOpt 2.0 are stored in the `xxx_macro.hpp` files. However, modifying these parameters manually may be inconvenient. Therefore, we now provide a Python script that helps you easily update these parameters.

We have added a `templates` folder under `packages/application/cvrp/`, which contains three `.txt` template files. You can modify these templates according to your needs. A sample template file looks like this:

```angular2html
app_type APPLICATION_TYPE::CVRP            // Specifies CVRP application type
ml_type  ML_TYPE::ML_USE_MODEL             // Use 2LBB method
IF_WRITE_NODE_OUT false                    // Do not write node output
NODE_MEMORY_ROUTE_LENGTH 20                // Use arc-memory-based rank-1 cuts if avg. route length > 20
IF_USE_STAB false                          // Do not use stabilization techniques
TIME_LIMIT 3600                            // Time limit set to 3600 seconds
````

If you encounter memory-related issues, refer to the troubleshooting section in our documentation and manually adjust parameters. Alternatively, you can add:

```angular2html
LABEL_ASSIGN 800000
```

in your template file to set a smaller memory request.

After editing your template file (`xxx.txt`), navigate to the CVRP application directory:

```bash
cd RouteOpt/packages/application/cvrp/
```

Then run:

```bash
python3 inspect_code.py xxx.txt
```
Your code will now compile with the updated parameters.


**Note:**  
This Python script can only modify key parameters initialized with braces `{}` in the code, for example:

```cpp
constexpr double ArcEliminationTimeIncreased = 20;
constexpr double HardTimeThresholdInAllEnumeration = 50;
constexpr int LABEL_ASSIGN{1600000};
````

Parameters initialized using an equals sign `=` (such as `constexpr double parameter = value;`) are not modifiable by this script.

If you find parameters that should be adjustable but currently aren't, please feel free to open an issue or submit a pull request.

## Documentation

For powerful and detailed documentation, please visit our [RouteOpt Documentation](https://Zhengzhong-You.github.io/RouteOpt-Docs/
). Here, you will find step-by-step guides, API references, and tutorials covering every module of RouteOpt 2.0.

## Migration Notice

**Important:** RouteOpt 1.0 has been removed. All features and improvements have been integrated into RouteOpt 2.0.  
For users previously working with RouteOpt 1.0, please note that the repository has been moved to a new location:

- **RouteOpt1.0 Repository:** [RouteOpt1.0 Repository](https://github.com/Zhengzhong-You/RouteOpt1).

We strongly encourage all users to migrate to RouteOpt 2.0 to take advantage of the latest advancements.

## Community and Support

We welcome contributions, feedback, and community involvement. If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request. Join our community discussions to share your experiences and collaborate with fellow users.

---

Thank you for choosing **RouteOpt 2.0**. We hope you find it as powerful and intuitive as we do. Happy optimizing!
