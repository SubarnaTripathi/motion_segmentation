%based on code of Timothy Cour
function [W,C]=compute_w_c_multiscale(image,pvec);
% compute multiscale image affinity matrix W and multiscale constraint
% matrix C from input image

[p,q,r] = size(image);
dataW = computeparametersw(image,pvec);
layers = computeparameterslayers(p,q); %only depends on image size

[C,C12]=computemultiscaleconstraints(layers);

% compute each layers(i).location as subsamples of the finest layer
layers=computelocationfromconstraints(C12,layers);

%W=computemultiscalew(image,layers,dataW);
W=computemultiscalew(image,layers,dataW);
