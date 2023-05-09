<h1 align="center">
<img src="qand_logo.jpg" width="200" style="border-radius:20%">
</h1>

# QAND
Quantitative and Algebraic Nonlinear Dynamics
## Dependencies

* Python (Preferred with a [conda](https://docs.conda.io/en/latest/miniconda.html) environment)
    ```Python
    conda create -n qand_env python=3.10 # optional
    conda activate qand_env # optional
    pip install --upgrade pip setuptools wheel # optional
    ```
* PyJulia ([Installation Guide](https://pyjulia.readthedocs.io/en/latest/installation.html))
* diffeqpy ([Installation Guide](https://pypi.org/project/diffeqpy/))

Now you can Install QAND. Change directory to the package parent folder and run this:
* `pip install -e .`

## Features

So far, these items have been developed

* System of Ordinary Differential Equations
    * Defining in Julia script.
    * Time Series/Trajectory from an initial point.
    * Bifurcation analysis with multiple nonlinear features (local min/max peaks, min/max intervals)
    * Save/Load for Differential Equation/Trajectory/Bifurcation Objects.
    * Plotting

TODO:
* Plotting Module
* System of Ordinary Differential Equations
    * Algebraic Definition
    * Solving algebraically for equilibrium points
    * Lyapunov exponent analysis
    * Stability Analysis
    * Basin of Attraction Analysis
    * Predefining Famous Systems
    * Save/Load for new objects
* Maps
    * Definition in Numba
    * Sequence series
    * Bifurcation
    * Save/Load
    * Predefining Famous Maps
* Complex Networks
    * Definition
    * Coupling Matrix and Topology analysis
    * Predefining Famous Topologies
    * Synchronization / Chimera states detection
    * MSF analysis
    * Save/Load

## Name
QAND (/qænd/) Means sugar cube in Persian. You can enjoy them with some hot saffron tea 🍵.
