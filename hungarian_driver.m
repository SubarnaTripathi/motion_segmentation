%Andrew Rabinovich, amrabino@ucsd.edu
function labels=hungarian_driver(temp_labels,nbCluster,anchor)

%initially labels in a matrix whose columns represent labeling for 
%each iteration. With this code, an optimum assignment is going to 
%be given. It is required to choose an anchor, for example iteration
%1 could have all the rest aligned to it.
[dim, iters]=size(temp_labels);

%convert represendation to binary from decimal
% for(i=1:iters)
%     temp=zeros(dim,nbCluster);
%     for(j=1:nbCluster)
%         x=find(labels(:,i)==j);
%         temp(x,j)=1;
%     end
%     temp_labels{i}=temp;
% end

%this is the anchor. ideally, instead of trying a random anchor, need to try all the achnors, and then
% pick the max
label1=temp_labels{anchor};
[x,y]=size(label1);
if(y==1)
    temp=zeros(x,max(label1));
    for(q=1:x)
        temp(q,label1(q))=1;
    end
    label1=temp;
end
labels_new{1}=label1;

align_list=[2:iters];

for(i=1:iters-1)
    label2=temp_labels{align_list(i)};
    
    [x,y]=size(label2);
    if(y==1)
        temp=zeros(x,max(label2));
        for(q=1:x)
            temp(q,label2(q))=1;
        end
        label2=temp;
    end
    
    A=label1'*label2;
       
    [C,T]=hungarian(A);
    for(j=1:nbCluster)
        x=find(A(j,:)==max(A(j,:)));
        temp(:,j)=label2(:,x(1));
    end
    labels_new{align_list(i)}=temp;
end
%%%%%%%%%convert back to decimal from binary%%%%%%
labels=zeros(dim,iters);
for(i=1:iters)
    for(j=1:nbCluster)
        x=find(labels_new{i}(:,j)==1);
        labels(x,i)=j;
    end
end

