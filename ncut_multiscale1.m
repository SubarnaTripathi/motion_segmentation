%based on code of Timothy Cour
function [classes,X,lambda,Xr,W,C] = ncut_multiscale1(image_x,image_y, nsegs,pvec);

%% compute multiscale affinity matrix W and multiscale constraint matrix C
[W,C]=compute_w_c_multiscale1(image_x,image_y, pvec);

% fprintf('size of W =%d\n', size(W));
% pause

%% compute constrained normalized cuts eigenvectors
if ~isempty(C)
    [X,lambda,timing] = computencutconstraint_projection(W,C,nsegs);
else
    [X,lambda,timing] = computekfirsteigenvectors(W,nsegs);
end

%% compute discretisation
[p,q,r]=size(image_x);
indPixels = (1:p*q)';
X = reshape(X(indPixels,:),p,q,nsegs);
[classes,Xr] =discretisation(X);

