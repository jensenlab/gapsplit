
% ============== Create linear model ==============

%             1  2  3  4  5  6  7  8  9 10
cobra.S =  [ -1  0  0  0 -1  0  0  0  0  0;   % A
              0 -1  0 -1  0  0  0  0  0  0;   % C
              0  0 -1  0  0 -1  0  0  0  0;   % F
              0  0  0  0  1  1 -1  0  0  0;   % B
              0  0  0  1  0  0 -1 -1  0  0;   % D
              0  0  0  0  0  0  1  0 -1  0;   % E
              0  0  0  0  0  0  0  1  0 -1 ]; % G

cobra.lb = [ -1 -1 -1 -1 -1  0 -1  0  0 -1 ]';
cobra.ub = [  1  1  1  1  1  1  1  1  1  1 ]';

cobra.c  = [  0  0  0  0  0  0  0  0  1  0 ]';

[m,n] = size(cobra.S);

cobra.b = zeros(m,1);

cobra.rxns = arrayfun(@(x) sprintf('r%i',x),1:n,'Uniform',false);
cobra.mets = {'A','C','F','B','D','E','G'}';

cobra.genes = {'g4','g5a','g5b', ...
               'g6','g7a','g7b','g8a','g8b'}';

cobra.grRules = {'';
                 '';
                 '';
                 'g4';
                 'g5a & g5b';
                 'g6';
                 'g7a & g7b';
                 'g8a | g8b';
                 '';
                 ''};

cobra.rules = {'';
               '';
               '';
               'x(1)';
               'x(2) & x(3)';
               'x(4)';
               'x(5) & x(6)';
               'x(7) | x(8)';
               '';
               ''};

% ============== Create nonlinear model ==============

model = cobra;
x0 = zeros(size(model.lb));

options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);

logit = @(x) 1/(1 + exp(-x));
nonlincon = @(x) deal( x(8) - logit(-x(1)), model.S*x - model.b );

problem = struct();
problem.objective = @(x) 0;
problem.x0 = x0;
problem.lb = model.lb;
problem.ub = model.ub;
problem.nonlcon = nonlincon;
problem.solver = 'fmincon';
problem.options = options;

% ============== Sampling ==============

sampling = gapsplit_nonlinear(problem,100,'secondaryFrac',0.1);

%writetable(table([(1:100)' sampling.coverage]),'nonlinear_coverage.csv')
