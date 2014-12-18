%Estimate stability between a number of clusterings, based on stability defined in 
%Rabinovich et al. CVPR 06. 
%Andrew Rabinovich amrabino@ucsd.edu

function final_stability=comp_clust(labels,num_seg);

[num_points,num_iter]=size(labels);

%incorporates the k-means scores to weigh permulations
%pos=find(error==min(error));
%weight=zeros(length(error),1);
%weight=error(pos(1))./error;
weight=ones(num_iter,1);
for(i=1:num_points)
    for(j=1:num_seg)
        temp=find(labels(i,:)==j);
        count(j)=length(temp);
    end
    [m,p]=max(count);
    tmp=zeros(num_iter,1);
    tmp(find(labels(i,:)==p))=1;
    stability(i)=sum(tmp.*weight)/sum(weight);
end
final_stability=(sum(stability)-num_points/num_seg)/(num_points-num_points/num_seg);



%%%%%%%%%%%pairwise%%%%%%%%%%%%%%
% for(i=1:num_points)
%     for(j=1:num_seg)
%         temp=find(labels(i,:)==j);
%         count(j)=length(temp);
%     end
%     stability(i)=max(count)/num_iter;%this give a values between 0 and 1
%     %keyboard;
% end
% final_stability=(sum(stability)*num_seg-sum(stability))/(num_points*num_seg);