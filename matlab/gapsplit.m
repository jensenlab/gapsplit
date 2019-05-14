function [sampling] = gapsplit(model,n,varargin)
% GAPSPLIT  Random sampling for COBRA models
%
%   [SAMPLING] = GAPSPLIT(MODEL,N,...parameters...)
%
%   GAPSPLIT uniformly samples LP, QP, and MILP models by targeting
%   "gaps" in the sample space. GAPSPLIT can give better coverage with
%   with fewer samples than ACHR-family solvers, although each solution
%   requires more computation.
%
% INPUTS
%   model   A COBRA model structure.
%   N       Number of samples to generate.
%
% PARAMETERS
%   'primary'        Strategy for selection the primary target. Targets are
%                    chosen sequentially ('seq', default), randomly 
%                    ('random'), or by always targeting the variable with
%                    the largest relative gap ('max').
%   'secondaryFrac'  Fraction of model variables randomly chosen as
%                    secondary targets during each iteration. Default is
%                    0.05 (5% of reactions). If 0, no secondary targeting
%                    is used; this may decrease coverage but improves
%                    runtime for numerically difficult models.
%   'vars'           Indices of variables to target during sampling. By
%                    default all variables are targeted. Regardless of
%                    targeting, all variables are sampled.
%   'maxTries'       Sampling attempts that return infeasible or unbounded
%                    solutions are discarded. Thus the total number of
%                    samples may exceed N for difficult models. 'maxTries'
%                    limits the number of attempts. Gapsplit will terminate
%                    before N samples are found in 'maxTries' is reached.
%                    The default is Inf, so GAPSPLIT will try until all N
%                    samples are found.
%   'minRange'       Variables are split only if their feasible range
%                    is larger than this value (default 1e-5).
%   'primaryTol'     The primary target is split by setting the upper and
%                    lower bounds to the midway point of the max gap.
%                    The bounds are set to within +/- primaryTol of the
%                    midway point to avoid infeasible solutions due to
%                    numerical issues.
%   'minval',        minval and maxval are vectors holding the feasible
%      'maxval'      range for each variable. If these are not provided,
%                    GAPSPLIT will use FLUXVARIABILITY to find the bounds.
%   'enforceRange'   If true (default), round solutions to fall within the
%                    feasible range. This prevents small deviations outside
%                    the feasible range from causing small decreases in
%                    coverage.
%   'reportInterval' Show the coverage and gap statistics at this interval.
%                    If a number between 0.0 and 1.0 is given, GAPSPLIT
%                    reports when that fraction of N samples is finished
%                    (i.e. if N=1000 and reportInterval=0.1, reports are
%                    printed every 100 samples.) To turn off reporting, set
%                    to 0. Default is 0.1.
%   'debug'          Turn on debugging information. Default is false.
%   'solverParams'   Cell array of parameters passed to solver.
%
% OUTPUTS
%   A sampling structure with the following fields:
%   samples           Array of sample points. Columns are variables in the
%                     model. Each row is a sample.
%   coverage          Fractional coverage after each sample. Coverage is
%                     defined as 1 - mean(relative maxgap) across all 
%                     variables. The relative maxgap for a variable is the
%                     size of the maxgap divided by the variable's feasible
%                     range.
%   minGap,medianGap, Minimum, median, and maximum gap across all
%      maxGap         variables; calculated after each sample.
%   elapsed           Time (in seconds) since sampling started until each
%                     sample was found.
%   minval,maxval     Feasible range for each variable.

global CBT_QP_SOLVER;

param = inputParser;
param.addParameter('minval',[],@isnumeric);
param.addParameter('maxval',[],@isnumeric);
param.addParameter('enforceRange',true,@islogical);
param.addParameter('secondaryFrac',0.05,@(x) 0.0 <= x && x <= 1.0);
param.addParameter('reportInterval',0.1);
param.addParameter('vars',[]);
param.addParameter('minRange',1e-5);
param.addParameter('primaryTol',1e-3);
param.addParameter('debug',false);
param.addParameter('maxTries',Inf);
param.addParameter('solverParams',{});

errMsg = "Primary selector must be 'seq', 'max', or 'random'.";
isValidPrimary = @(x) assert(ismember(x,{'seq','max','random'}),errMsg);
param.addParameter('primary','seq',isValidPrimary);

param.parse(varargin{:});

debug = @(s) param.Results.debug && fprintf(s);
primaryTol = param.Results.primaryTol;
maxTries = param.Results.maxTries;

reportInterval = param.Results.reportInterval;
if reportInterval < 1.0
    % fraction of total samples
    reportInterval = floor(reportInterval * n);
end

minval = param.Results.minval(:);
maxval = param.Results.maxval(:);
enforceRange = param.Results.enforceRange;
if isempty(minval) || isempty(maxval)
    if reportInterval > 0
        fprintf('Calculating feasible ranges for variables.');
    end
    tstart = tic;
    [minval,maxval] = fluxVariability(model,0);
    if reportInterval > 0
        fprintf(' (%.2f seconds)\n', toc(tstart));
    end
elseif reportInterval > 0
    fprintf('Using supplied feasible ranges for variables.\n');
end

model = buildLPproblemFromModel(model);
[m,p] = size(model.A);

vars = param.Results.vars;
if isempty(vars)
    vars = 1:p;
end
fvaRange = maxval - minval;
vars = vars(fvaRange(vars) >= param.Results.minRange);
Nvars = length(vars);

quadWeights = 1 ./ fvaRange.^2;

% add constraints for secondary target
ksec = floor(param.Results.secondaryFrac * Nvars);
useSec = ksec >= 1;
kcidx = [false(p,1); true(ksec,1)];
kridx = [false(m,1); true(ksec,1)];
model.A(kridx,kcidx) = -eye(ksec);
model.b(kridx) = 0;
model.c = zeros(p+ksec,1);
model.lb(kcidx) = -Inf;
model.ub(kcidx) = Inf;
model.csense(kridx) = 'E';

model.F = sparse(p+ksec,p+ksec);
model.F(kcidx,kcidx) = 1;

sampling.samples = zeros(n,p);
sampling.coverage = zeros(n,1);
sampling.minGap = zeros(n,1);
sampling.medianGap = zeros(n,1);
sampling.maxGap = zeros(n,1);
sampling.elapsed = zeros(n,1);
sampling.minval = minval;
sampling.maxval = maxval;

[hdr,fmt] = makeFormatStrs(n);
if reportInterval > 0
    fprintf('Targeting 1 primary and %i secondary variables.\n',ksec);
    fprintf('Sampling LP model with %i/%i unblocked variables (%6.2f%%).\n',Nvars,p,Nvars/p*100);
    fprintf(hdr);
end

tstart = tic;
[~,l,u,r] = maxgap([minval'; maxval']);
tries = 0;
i = 0;
while tries < maxTries && i < n
    prevModel = model;
    
    % fix the primary target to split the max gap
    switch param.Results.primary
        case 'seq'
            primary = rem(i,Nvars);
            if primary == 0
                primary = vars(Nvars);
            else
                primary = vars(primary);
            end
        case 'max'
            [~,idx] = max(r(vars));
            primary = vars(idx);
        case 'random'
            primary = vars(randi(Nvars,1));
    end
    primaryTarget = mean([l(primary) u(primary)]);
    targetRange = sort([1-primaryTol 1+primaryTol]*primaryTarget);
    model.lb(primary) = targetRange(1);
    model.ub(primary) = targetRange(2);
    debug(sprintf('Primary var %i set to [%f,%f] from range [%f,%f].\n', ...
        primary, model.lb(primary), model.ub(primary), minval(primary), maxval(primary)));
    
    if useSec
        secondary = vars(randi(Nvars,ksec,1));
        model.A(kridx,secondary) = eye(ksec);
        model.b(kridx) = mean([l(secondary); u(secondary)]);
        model.F(kcidx,kcidx) = diag(quadWeights(secondary));
    end
    
    try
        if strcmpi(CBT_QP_SOLVER,'gurobi')
            % barrier or automatic methods lead to residual problems in 
            % Gurobi. Use primal/dual simplex instead.
            sol = solveCobraQP(model,'Method',1,param.Results.solverParams{:});
        else
            sol = solveCobraQP(model,param.Results.solverParams{:});
        end
        debug(sprintf('   Solution was %s (%i) in %f seconds with objective %f.\n', ...
            sol.origStat, sol.stat, sol.time, sol.obj));
    catch
        sol = struct('stat',-1);
    end
    
    model = prevModel;
    tries = tries + 1;
    if sol.stat ~= 1
        continue
    end
    i = i + 1;
    raw = sol.full(1:p);
    if enforceRange
        raw = min(raw, maxval);
        raw = max(raw, minval);
    end
    sampling.samples(i,:) = raw;
    
    tcurrent = toc(tstart);
    [~,l,u,r] = maxgap([minval'; maxval'; sampling.samples(1:i,:)]);
    sampling.coverage(i) = 1 - mean(r(vars));
    sampling.minGap(i) = min(r(vars));
    sampling.medianGap(i) = median(r(vars));
    sampling.maxGap(i) = max(r(vars));
    sampling.elapsed(i) = tcurrent;
    
    if mod(i,reportInterval) == 0
        % report
        fprintf(fmt, i, n, sampling.coverage(i)*100, sampling.minGap(i), ...
            sampling.medianGap(i), sampling.maxGap(i), ...
            tcurrent, tcurrent/i*(n - i), tries - i);
    end
end

if i < n
    % hit maxTries; return only found samples
    sampling.samples = sampling.samples(1:i,:);
    sampling.coverage = sampling.coverage(1:i);
    sampling.minGap = sampling.minGap(1:i);
    sampling.medianGap = sampling.medianGap(1:i);
    sampling.maxGap = sampling.maxGap(1:i);
    sampling.elapsed = sampling.elapsed(1:i);
    
    if reportInterval > 0
        fprintf('\n*gapsplit terminated after %i attempts with %i solutions.\n',maxTries,i);
    end
end

end

function [width,lower,upper,rel] = maxgap(X)
% rows of X are samples
% columns of X are variables
n = size(X,1);
X = sort(X); % each col in ascending order
gaps = X(2:n,:) - X(1:n-1,:);
[width,idx] = max(gaps,[],1);
linIdx = idx + (0:length(idx)-1)*n;
lower = X(linIdx);
upper = X(linIdx+1);
rel = width ./ (X(n,:) - X(1,:));
end

function [hdr,fmt] = makeFormatStrs(n)
wN = length(int2str(n));
fracWidth = 2*wN + 1;
samplehdr = 'Samples';
samplefmt = ['%' int2str(wN) 'i/%' int2str(wN) 'i'];
if fracWidth < length(samplehdr)
    samplefmt = [repmat(' ',1,length(samplehdr)-fracWidth) samplefmt];
elseif fracWidth > length(samplehdr)
    samplehdr = [repmat(' ',1,fracWidth-length(samplehdr)) samplehdr];
end
hdr = ['\n' samplehdr '   Coverage   MinGap   Median   MaxGap     Elapsed     Remaining   Infeasible\n'];
fmt = [samplefmt '    %6.2f%%   %6.4f   %6.4f   %6.4f   %9.2f     %9.2f   %10i' '\n'];
end


