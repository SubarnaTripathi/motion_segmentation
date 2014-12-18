%based on code of Timothy Cour
%modification by Andrew Rabinovich amrabino@ucsd.edu
function dataW = computeParametersW1(image_x,image_y, pvec);
% sets parameters for computing multiscale image affinity matrix W
% dataW.edgeVariance: edge variance for texture cue (1-p)
% dataW.sigmaI: intensity variance for intensity cue (p)

% dataW.edgeVariance=0.5*(1-pvec)+eps;
% dataW.sigmaI=0.5*pvec+eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataW.xVariance = 0.5*pvec+eps; %%%%%% (flow x)
dataW.yVariance=0.5*pvec+eps; %%%%%%% (flow y)