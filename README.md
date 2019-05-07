# gapsplit

Efficient sampling for COBRA models. `gapsplit` samples both convex (LP) and non-convex (MILP) models by targeting unexplored regions of the solution space.

`gapsplit` is available for Matlab and Python.

## Matlab

### Installation

Download the `gapsplit` repository and add the `matlab` directory to Matlab's path. `gapsplit` requires the [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) with a QP solver (or MIQP solver to sample MILP models).

Test your `gapsplit` installation with the `test_gapsplit` script. (Be sure to initialize the COBRA Toolbox first.)
```matlab
>> initCobraToolbox
>> test_gapsplit
Calculating feasible ranges for variables. (0.45 seconds)
Targeting 1 primary and 1 secondary variables.
Sampling LP model with 10/10 unblocked variables (100.00%).

Samples   Coverage   MinGap   Median   MaxGap     Elapsed     Remaining   Infeasible
 10/100     77.47%   0.1879   0.2499   0.2502        0.04          0.37            0
 20/100     88.10%   0.0943   0.1250   0.1258        0.07          0.28            0
 30/100     89.45%   0.0624   0.1234   0.1241        0.10          0.24            0
 40/100     94.04%   0.0476   0.0626   0.0629        0.13          0.20            0
 50/100     94.38%   0.0463   0.0623   0.0627        0.16          0.16            0
 60/100     94.70%   0.0331   0.0607   0.0623        0.19          0.13            0
 70/100     95.77%   0.0312   0.0407   0.0612        0.23          0.10            0
 80/100     96.61%   0.0278   0.0319   0.0461        0.26          0.06            0
 90/100     96.72%   0.0238   0.0314   0.0455        0.29          0.03            0
100/100     97.03%   0.0236   0.0312   0.0313        0.32
```

### Using `gapsplit`

`gapsplit` is a single function called `gapsplit`. For information on usage, see `help gapsplit`.

## Python

## Citation

Keaty TC, Jensen PA. **gapsplit: Efficient random sampling for non-convex constraint-based models.** *bioRxiv*.
