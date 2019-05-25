function [cobra] = tiger_to_cobra_milp(tiger)
% TIGER_TO_COBRA_MILP  Convert a TIGER model to COBRA MILP structure

% first let's check if this is already a COBRA MILP model
needed_names = {'A','b','c','lb','ub','osense','csense','vartype'};
if isfield(tiger,needed_names)
    cobra = tiger;
    return
end

% TIGER models have bounds and indicator constraints that are not in the
% A matrix. These are added right before solution by CMPI.PREPARE_MIP.
tiger = cmpi.prepare_mip(tiger);

cobra = struct();
cobra.A = tiger.A;
cobra.b = tiger.b;
cobra.c = tiger.obj;
cobra.lb = tiger.lb;
cobra.ub = tiger.ub;
cobra.osense = tiger.sense; % -1 = max, 1 = min

% TIGER uses <, =, > for constraint types; COBRA uses 'L','E', and 'G'.
cobra.csense = regexprep(tiger.ctypes(:)',{'<','=','>'},{'L','E','G'});

% TIGER uses 'c','b','i' for continuous, binary, and integer variables.
% COBRA uses the same letters but in UPPERCASE.
cobra.vartype = upper(tiger.vartypes);
