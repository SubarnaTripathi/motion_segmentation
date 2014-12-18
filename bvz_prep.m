function W=bvz_prep(I1,I2,Diff,nwarps)
% 

for count=1:(nwarps)
  eval(['Diff' int2str(count) '=Diff(:,:,' int2str(count) ');']);
end

for count=1:(nwarps)
  eval(['Diff' int2str(count) '(isnan(Diff' int2str(count) '))=exp(10);']);
  eval(['IW' int2str(count) '=abs(Diff' int2str(count) ');']);
end
%eval(['IW' int2str(nwarps+1) '=abs((I2-I1));']);

W=zeros([length(IW1(:)) nwarps]);
for count=1:(nwarps)
  eval(['W(:,' int2str(count) ')=IW' int2str(count) '(:);']);
end

W=W';
