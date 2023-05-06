# QAND
Quantitative and Algebraic Nonlinear Dynamics

## Dependencies

* Python (Preferred with a [conda](https://docs.conda.io/en/latest/miniconda.html) environment)
    ```Python
    conda create -n qand_env python=3.10 # optional
    conda activate qand_env # optional
    pip install --upgrade pip setuptools wheel # optional
    pip install numba # optional. if you install jupyter, install this first for numpy compatibility
    pip install jupyter # optional for notebook use
    ```
* numpy and matplotlib
    ```Python
    pip install matplotlib scipy # automatically installs numpy as well
    ```
* PyJulia ([Installation Guide](https://pyjulia.readthedocs.io/en/latest/installation.html))
* diffeqpy ([Installation Guide](https://pypi.org/project/diffeqpy/))

Now you can Install QAND. Change directory to the package parent folder and run this:
* `pip install -e .`

## Name
QAND (/qænd/) Means Sugar Cubes in Persian.
![QAND Logo](qand_logo.jpg)
