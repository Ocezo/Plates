% Special k-means - 12/02/14
% Jean-Marc Berthommé
%
% - 12/02/14:
%   . special kmeans that updates only the "Ku" first clusters. The other
%   ones stay fix.

% Step 0: Provide a data set X, a number of clusters K and eventually an
% initialization set Xi
function [C, Xm, it] = kmeans(X, Xi, Ku)
stop = false; % flag

% Step 1: No random init° here, Xi defines the "K" first clusters

% Step 2: Assign each object to the nearest cluster
d = euc_dist(double(Xi),double(X));
[~,C] = min(d,[],1); C = C'; % classes vector

it = 0;
while ~stop
    % Step 3: Compute the cluster means
    Xm = Xi; % Xmean
    for k=1:Ku, Xm(k,:) = mean(X(C==k,:), 1); end;
    
    % Step 4: Assign each object to the nearest cluster
    Cold = C; d = euc_dist(double(Xm),double(X));
    [~,C] = min(d,[],1); C = C';
    
    it = it + 1;
    if isequal(C, Cold), stop = true; end;
end

function d = euc_dist(X,Y)
nx = size(X,1); ny = size(Y,1);
d = sqrt(sum(X.^2,2)*ones(1,ny) + ones(nx,1)*sum(Y.^2,2)' - 2*(X*Y'));
