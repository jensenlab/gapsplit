function [cover,steps] = coverage(X,stepsize,mins,maxs)

n = size(X,1);

if nargin > 1
    X = [mins(:)'; X];
end
if nargin > 2
    X = [maxs(:)'; X];
end
bndpts = nargin - 2; % number of added points

steps = 1:stepsize:n;
cover = zeros(size(steps));
for i = 1:length(steps)
    [~,~,~,rel] = maxgap(X(1:(steps(i)+bndpts),:));
    cover(i) = 1 - mean(rel(isfinite(rel)));
end
%cover = cover(bndpts+(1:n),:);

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
