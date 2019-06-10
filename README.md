# gapsplit

Efficient random sampling for constraint-based (COBRA) models. `gapsplit` samples both convex (LP) and non-convex (MILP) models by targeting unexplored regions of the solution space.

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
100/100     97.03%   0.0236   0.0312   0.0313        0.32          0.00            0
```

### Using `gapsplit`

`gapsplit` is a single function called `gapsplit`. For more information, see `help gapsplit`.

### Sampling MILP models

`gapsplit` can sample non-convex models. The `gapsplit` function recognizes MILP models by the `vartype` field and uses the COBRA Toolbox MIQP solver.

To sample [TIGER](https://github.com/pauljensen/tiger), pass the model to `gapsplit` with the `'tiger'` option set to `true`. The TIGER model will be converted to a COBRA MILP model. TIGER needs to be installed for this conversion.

## Python

The `gapsplit.py` module is compatible with [cobrapy](https://opencobra.github.io/cobrapy). To test your installation, run the module as a script.
```sh
$ python3 gapsplit.py
Calculating feasible ranges using FVA.
Targeting 87/95 unblocked primary variables.
Targeting 4 secondary variables.

 Sample   Coverage   MinGap   Median   MaxGap     Elapsed     Remaining   Infeasible
 10/100     38.47%   0.2500   0.6096   0.9434        0.33          2.96            0
 20/100     59.19%   0.1292   0.3997   0.6743        0.62          2.50            0
 30/100     67.75%   0.1219   0.3227   0.5010        0.92          2.14            0
 40/100     72.82%   0.1112   0.2635   0.4019        1.20          1.80            0
 50/100     77.77%   0.1022   0.2369   0.3352        1.50          1.50            0
 60/100     80.07%   0.0625   0.1948   0.2906        1.80          1.20            0
 70/100     80.79%   0.0625   0.1874   0.2505        2.09          0.90            0
 80/100     83.56%   0.0550   0.1677   0.2420        2.39          0.60            0
 90/100     84.79%   0.0540   0.1574   0.2094        2.68          0.30            0
100/100     85.59%   0.0540   0.1498   0.1948        2.96          0.00            0
```

`gapsplit` can use the Gurobi interface directly to significantly reduce sampling times. Use the `gurobi_direct=True` parameter to enable this feature.

## Citation

Keaty TC, Jensen PA. **gapsplit: Efficient random sampling for non-convex constraint-based models.** *bioRxiv*.
