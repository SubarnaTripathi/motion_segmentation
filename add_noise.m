%Andrew Rabinovich
function image=add_noise(image,sigma);
rand('state',sum(100*clock));
[x,y]=size(image);

rand('state',sum(100*clock));
r=randn(x,y)*sigma;
%keyboard;
image=image+r*sigma;
image(image>1)=1;
image(image<0)=0;
