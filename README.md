# Normal form analysis

This code analyzes finite normal form games.

**Author: Vitor Farinha Luz\
Email: [vitor.farinhaluz@ubc.ca](mailto:vitor.farinhaluz@ubc.ca)**


The files contained here are the following:
1. "Feasible sets.ipynb": Jupyter notebook generating plots for feasible payoff sets
1. "Correlated equilibria.ipynb": Jupyter notebook generating plots for correlated equilibria
1. Manifest and Project.toml: auxiliary files needed to run Julia code
1. README.md: you are reading it
1. src folder: contains actual code doing calculations and plots:
    1. Games_src.jl: contains the main code doing the work in the background
    1. CorrelatedEqPolytopes.jl: adds to module calculation of correlated equilibria using linear algebra.
1. .gitignore: ignore it

The code can be run on a virtual machine using the following link:
https://mybinder.org/v2/gh/vfarinhaluz/Julia-games/HEAD

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vfarinhaluz/Julia-games/HEAD)

If you find any bugs, please [create an issue entry](https://github.com/vfarinhaluz/Julia-games/issues) explaining the problem and the code leading to it, or submit a pull request.