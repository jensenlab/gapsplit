
import sys
import time

import numpy as np
import pandas as pd
import cobra
import cobra.test

def gapsplit(
        model, n, max_tries=None,
        primary='max', primary_tol=0.001,
        secondary_frac=0.05,
        fva=None,
        min_range=1e-5,
        enforce_range=True,
        report_interval=0.1,
        quiet=False):
    """Randomly sample a COBRA model.

    Parameters
    ----------
    model: cobra.Model
        The model to sample. The model will not be modified during sampling.
    n: integer
        Number of samples to generate
    max_tries: integer, optional, default=None
        Sampling attempts that return infeasible or unbounded solutions are
        discarded. Thus the total number of optimizations may exceed `n` for
        difficult models. `max_tries` limits the total number of attempts. If
        None (default), gapsplit will continue until `n` feasible samples are
        found.
    primary: str, optional, default='max'
        Strategy for selection the primary target. Targets are chosen
        sequentially ('seq'), randomly ('random'), or by always targeting the
        variable with the largest relative gap ('max', default).
    primary_tol: float, optional, default=0.001
        The primary target is split by setting the upper and lower bounds to
        the midway point of the max gap. The bounds are set to within +/-
        `primary_tol` times the width of the gap to avoid infeasible solutions
        due to numerical issues.
    secondary_frac: float, optional, default=0.05
        Fraction of model variables randomly chosen as secondary targets during
        each iteration. Default is 0.05 (5% of reactions). If 0, no secondary
        targeting is used; this may decrease coverage but improves runtime for
        numerically difficult models.
    fva: pandas.DataFrame, optional, default=None
        gapsplit uses flux variability analysis (FVA) to find the feasible
        ranges for each variable. The user can supply the output of a previous
        `cobra.flux_analysis.flux_variability_analysis` run to avoid re-running
        FVA. If None (default), gapsplit will run FVA.
    min_range: float, optional, default=1e-5
        Variables are targeted only if their feasible range is larger than
        this value.
    enforce_range: boolean, optional, default=True
        If true (default), round solutions to fall within the feasible range.
        This prevents small deviations outside the feasible range from causing
        small decreases in coverage.
    report_interval: float or int, optional, default=0.1
        Show the coverage and gap statistics at this interval. If a number
        between 0.0 and 1.0 is given, gapsplit reports when that fraction of
        `n` samples is finished (i.e. if N=1000 and reportInterval=0.1, reports
        are printed every 100 samples.) To turn off reporting, set to 0.
    quiet: boolean, optional, default=True
        Set to false to keep gapslit from printing status updates.

    Returns
    -------
    pandas.DataFrame
        A data frame with rows = samples and columns = reactions. This is the
        same format as the other cobrapy samplers.
    """
    # output has rows = samples, columns = variables
    # cobrapy returns a pandas DF

    if quiet:
        report = lambda s: None
    else:
        report = lambda s: print(s)

    reactions = model.reactions

    if fva is None:
        report("Calculating feasible ranges using FVA.")
        fva = cobra.flux_analysis.flux_variability_analysis(
                model, reactions, fraction_of_optimum=0.0)
    else:
        report("Using supplied FVA ranges.")

    if secondary_frac >= 1.0:
        n_secondary = secondary_frac
    else:
        n_secondary = np.floor(secondary_frac * len(model.reactions)).astype(int)

    # only split reactions with feasible range >= min_range
    idxs = (fva.maximum - fva.minimum >= min_range).to_numpy().nonzero()[0]
    weights = (1/(fva.maximum - fva.minimum)**2).to_numpy()

    report("Targeting {}/{} unblocked primary variables.".format(len(idxs), len(model.reactions)))
    report("Targeting {} secondary variables.".format(n_secondary))

    report_header, report_format = _make_report_header(n)
    report("\n" + report_header)
    if report_interval < 1.0:
        report_interval = np.floor(report_interval * n).astype(int)

    samples = np.zeros((n, len(model.reactions)))
    k = 0
    infeasible_count = 0
    if primary == 'sequential':
        # primary_var will increment
        primary_var = -1
    if max_tries is None:
        # no limit on number of tries; pick a really big number
        max_tries = sys.maxsize
    start_time = time.time()
    for try_ in range(max_tries):
        relative, target, width = maxgap(samples[0:k,idxs], fva.iloc[idxs,:])
        if primary == 'max':
            primary_var = np.argmax(relative)
        elif primary == 'random':
            primary_var = np.random.choice(len(idxs), 1).astype(int)[0]
        elif primary == 'sequential':
            primary_var += 1
            if primary_var >= len(idxs):
                primary_var = 0

        primary_target = target[primary_var]
        primary_lb = primary_target - primary_tol*width[primary_var]
        primary_ub = primary_target + primary_tol*width[primary_var]

        secondary_vars = np.random.choice(len(idxs), n_secondary, replace=False)
        secondary_targets = target[secondary_vars]
        secondary_weights = weights[idxs[secondary_vars]]

        new_sample = _generate_sample(
            model, idxs[primary_var], primary_lb, primary_ub,
            idxs[secondary_vars], secondary_targets, secondary_weights)
        if new_sample is not None:
            if enforce_range:
                new_sample[new_sample > fva.maximum] = fva.maximum[new_sample > fva.maximum]
                new_sample[new_sample < fva.minimum] = fva.minimum[new_sample < fva.minimum]

            samples[k,:] = new_sample
            k += 1
            if k % report_interval == 0:
                elapsed = time.time() - start_time
                remaining = elapsed / k * (n - k)
                report(report_format.format(
                        i=k, n=n, cov=100*(1-np.mean(relative)),
                        min=np.min(relative), med=np.median(relative),
                        max=np.max(relative), ela=elapsed, rem=remaining,
                        inf=infeasible_count))
        else:
            infeasible_count += 1

        if k >= n: break

    if k < n:
        # max_tries reached; return fewer than n samples
        samples = samples[:k,:]

    return pd.DataFrame(data=samples,columns=fva.maximum.index)


def _generate_sample(
        model, primary_var, primary_lb, primary_ub,
        secondary_vars=None, secondary_targets=None, secondary_weights=None):
    """Formulate a [MI]QP to find a single solution."""
    with model:
        model.reactions[primary_var].lower_bound = primary_lb
        model.reactions[primary_var].upper_bound = primary_ub

        if secondary_vars is not None:
            quad_exp = 0
            for i, sec in enumerate(secondary_vars):
                diff = model.problem.Variable('difference_{}'.format(sec))
                cons = model.problem.Constraint(
                    model.reactions[sec].flux_expression - diff,
                    lb=secondary_targets[i], ub=secondary_targets[i])
                model.add_cons_vars([diff, cons])
                quad_exp += secondary_weights[i] * diff**2
            quad_obj = model.problem.Objective(quad_exp, direction='min')
            model.objective = quad_obj
        else:
            model.objective = model.problem.Objective(0)

        solution = model.optimize()
        if solution.status != 'optimal':
            return None
        else:
            return solution.fluxes


def _maxgap(points, fva=None):
    # points has rows = samples, columns = variables

    # make a copy because we're going to sort the columns
    points = points.copy()
    if fva is not None:
        points = np.vstack((fva.minimum, points, fva.maximum))
    points.sort(0)

    gaps = points[1:,:] - points[0:-1,:]
    width = gaps.max(0)
    loc = gaps.argmax(0)
    left = np.zeros(width.size)
    for i in range(width.size):
        left[i] = points[loc[i],i]
    relative = width / (points[-1,:] - points[0,:])
    target = left + width/2

    return relative, target, width


def _make_report_header(maxN):
    """Return the header and format string for reporting coverage."""
    nw = len(str(maxN))
    frac_width = 2*nw + 1  # width of 300/1000
    frac_header = 'Sample'
    frac_format = '{i:' + str(nw) + 'd}/{n:' + str(nw) + 'd}'
    if frac_width < len(frac_header):
        pad = ''.join([' ' for _ in range(len(frac_header) - frac_width)])
        frac_format = pad + frac_format
    elif len(frac_header) < frac_width:
        pad = ''.join([' ' for _ in range(frac_width - len(frac_header))])
        frac_header = pad + frac_header
    hdr = frac_header + "   Coverage   MinGap   Median   MaxGap     Elapsed     Remaining   Infeasible"
    fmt = frac_format + "    {cov:6.2f}%   {min:6.4f}   {med:6.4f}   {max:6.4f}   {ela:9.2f}     {rem:9.2f}   {inf:10d}"
    return hdr, fmt


if __name__ == '__main__':
    model = cobra.test.create_test_model('textbook')
    print(gapsplit(model,100,primary="max"))
