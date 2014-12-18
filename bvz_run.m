function assignment = bvz_run(warps,I,r,lambda);
%  assignment = bvz_run(warps,I,r,lambda);

warps(isnan(warps))=1e7;
warps=warps.*warps;
mask=fspecial('gaussian',2*r+1,r);
mask=lambda*mask(:);
N=prod(size(I));
iseg=im2col2(I+2,r);
cseg=im2col2(reshape(1:N,size(I)),r);
cseg=cseg-1;
creategraph('test1.graph',warps,I,cseg,mask,r)
%! ./smooth test1.graph
system('smooth test1.graph');
% load output;
% assignment = output;
load my_output;
assignment=my_output;


function V=im2col2(M,winrad);
[N1,N2]=size(M);
A=zeros(N1+2*winrad,N2+2*winrad);
A(winrad+1:winrad+N1,winrad+1:winrad+N2)=M(:,:);
V=im2col(A,[2*winrad+1,2*winrad+1],'sliding');
