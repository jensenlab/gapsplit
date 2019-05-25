function [sampling] = test_gapsplit(suite)

if nargin < 1
    suite = 'lp';
end

if ismember(suite,{'lp','all'})
    % =========== create a test COBRA model ===========

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

    cobra.rxns = array2names('r%i',1:n);
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

    sampling = gapsplit(cobra,100);
end

if ismember(suite,{'tiger','all'})
    %          r1 r2 r3 r4 r5
    cobra.S = [ 1 -1  0  0  0;  % A
                0  1  0 -1  0;  % B
                0  0  1 -1  0;  % C
                0  0  0  1 -1]; % D

    cobra.c = [ 0  0  0  0  1]';
    cobra.lb = -10*ones(size(cobra.c));
    cobra.ub =  10*ones(size(cobra.c));

    cobra.b = zeros(size(cobra.S,1),1);

    cobra.rxns = {'rxn1';'rxn2';'rxn3';'rxn4';'rxn5'};
    cobra.mets = {'A';'B';'C';'D'};

    cobra.genes = {'AB1';'AB2';'BCD1';'BCD2'};
    cobra.grRules = {'';'AB1 or AB2';'';'BCD1 and BCD2';''};
    cobra.rules = {'';'x(1) | x(2)';'';'x(3) & x(4)';''};

    tiger = cobra_to_tiger(cobra);
    sampling = gapsplit(tiger,100,'tiger',true);
end

