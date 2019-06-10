
import time

import numpy as np
import pandas as pd
import cobra
import cobra.test


def gapsplit(
        model, n, max_tries=None,
        primary='sequential', primary_tol=0.001,
        secondary_frac=0.05,
        fva=None,
        min_range=1e-5,
        enforce_range=True,
        report_interval=0.1,
        quiet=False,
        gurobi_direct=False):
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
    primary: str, optional, default='sequential'
        Strategy for selection the primary target. Targets are chosen
        sequentially ('sequential', default), randomly ('random'), or by always
        targeting the variable with the largest relative gap ('max').
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
    gurobi_direct: boolean, optional, default=False
        Use the gurobipy interface directly to sample the model. This can
        significantly reduce sampling times and gives identical results.
        Requires the gurobipy model and model.solver be set to 'gurobi'.

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

    if gurobi_direct:
        grb_model = _reduce_gurobi(model)

    samples = np.zeros((n, len(model.reactions)))
    k = 0
    infeasible_count = 0
    if primary == 'sequential':
        # primary_var will increment
        primary_var = -1
    try_ = 0
    start_time = time.time()
    while True:
        if max_tries is not None and try_ >= max_tries:
            break
        try_ += 1
        relative, target, width = _maxgap(samples[0:k,idxs], fva.iloc[idxs,:])
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

        if gurobi_direct:
            new_sample = _generate_sample_gurobi_direct(
                grb_model, idxs[primary_var], primary_lb, primary_ub,
                idxs[secondary_vars], secondary_targets, secondary_weights)
        else:
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


def _generate_sample_gurobi_direct(
        model, primary_var, primary_lb, primary_ub,
        secondary_vars=None, secondary_targets=None, secondary_weights=None):
    """Solve the model directly with the gurobipy interface.

    Cobrapy models have several features that makes them slow to sample. We can
    apply some Gurobi-specific optimizations to improve gapsplit runtimes.
        - "Unsplit" all the variables from forward and reverse to a single
          reaction that has positive and negative fluxes.
        - Set the objectives directly to skip the sympy interface.
        - Collection solutions for all variables simultaneously.

    The first optimization is handled by _reduce_gurobi. The latter tricks
    are used in this function.

    Inputs
    ------
    model: grb model object
        A copy of the Gurobi model object from cobra_mode._solver.problem.
        This model should be reduced using _reduce_gurobi.
    <all other inputs from _generate_sample>

    Returns
    -------
    Either a NumPy array with the new solution or None (if infeasible).
    """

    import gurobipy as grb
    # model is a copy of the grb model (from cobra_model._solver.problem)
    vars = model.getVars()
    primary = vars[primary_var]
    prev_lb = primary.getAttr("lb")
    prev_ub = primary.getAttr("ub")
    primary.setAttr("lb", primary_lb)
    primary.setAttr("ub", primary_ub)

    if secondary_vars is not None:
        qobj = grb.QuadExpr(0)
        sec_vars = [vars[i] for i in secondary_vars]
        qobj.addTerms(secondary_weights, sec_vars, sec_vars)
        qobj.addTerms(-2*secondary_weights*secondary_targets, sec_vars)
        model.setObjective(qobj, sense=grb.GRB.MINIMIZE)
    else:
        model.setObjective(grb.LinExpr(0))

    model.optimize()
    #if model.status == 3:
    #    print(model.status, primary_var, prev_lb, prev_ub, primary_lb, primary_ub)
    if model.status == 2:
        solution = np.array(model.getAttr("X", vars))
    else:
        solution = None

    # reset the bounds on the primary target
    primary.setAttr("lb", prev_lb)
    primary.setAttr("ub", prev_ub)
    model.update()

    return solution


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


def _reduce_gurobi(cobra):
    """Modify the gurobi model object to improve sampling efficiency."""
    grb = cobra._solver.problem.copy()
    varnames = [var.VarName for var in grb.getVars()]
    for rxn in cobra.reactions:
        if rxn.reverse_id in varnames:
            for_var = grb.getVarByName(rxn.id)
            rev_var = grb.getVarByName(rxn.reverse_id)
            ub = for_var.getAttr("ub") - rev_var.getAttr("lb")
            lb = for_var.getAttr("lb") - rev_var.getAttr("ub")
            for_var.setAttr("ub", ub)
            for_var.setAttr("lb", lb)
            grb.remove(rev_var)

    # For some reason FVA leaves a variable 'fva_old_objective' in the grb
    # model object. We need to remove it so the number of variables match the
    # number of reactions.
    if 'fva_old_objective' in varnames:
        grb.remove(grb.getVarByName('fva_old_objective'))

    grb.update()
    return grb


if __name__ == '__main__':
    model = cobra.test.create_test_model('textbook')
    gapsplit(model,100,primary="sequential")
