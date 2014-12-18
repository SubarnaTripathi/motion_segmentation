function ST_layersegmentation1
    addpath('mex');
    frame_type = 'jpg'; %'ppm'; %'jpg'

    if(nargin<3)
      r=2;
      lambda=0.285;
    end

    %bandwidth = 0.45; 
    
    %example = 'lovebird'; 
    example = 'stefan_cif';


    % set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
    alpha = 0.012;
    ratio = 0.75;
    minWidth = 20;
    nOuterFPIterations = 7;
    nInnerFPIterations = 1;
    nSORIterations = 30;

    para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

    frame_start = 3; %201
    frame_end = 5; %209

    movie_in = VideoWriter('input.avi');
    NS1 = VideoWriter('naive_assignment1.avi');
    NS2 = VideoWriter('naive_assignment2.avi');
    SS1 = VideoWriter('smooth_assignment1.avi');
    SS2 = VideoWriter('smooth_assignment2.avi');

    open(movie_in);
    open(NS1);
    open(NS2);
    open(SS1);
    open(SS2);
    
    extend = 0; %%%% augment the image data with optic flow (0 = only optic flow data)
    %%%% canonical ref frame for affine affinity
%     c_height = 64;
%     c_width = 64;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    for frame_pair = frame_start:frame_end
        % load the two frames
        if (strcmp(frame_type,'ppm'))
            img1 = sprintf('\\00%d.ppm', frame_pair);  %5
            img2 = sprintf('\\00%d.ppm', frame_pair+1); %6   
        elseif (strcmp(frame_type,'jpg'))
            img1 = sprintf('\\00%d.jpg', frame_pair);  %5
            img2 = sprintf('\\00%d.jpg', frame_pair+1); %6
        else
            disp('unknown image format... exiting');
            return;
        end

        im1 = im2double(imread([example img1])); %1.jpg
        im2 = im2double(imread([example img2]));

        I1 = rgb2gray(im1);
        I2 = rgb2gray(im2);
        
        [c_height, c_width, cl] = size(I1);

        % this is the core part of calling the mexed dll file for computing optical flow
        % it also returns the time that is needed for two-frame estimation
        tic;
        [vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
        toc
        % output gif
        clear volume;
        volume(:,:,:,1) = im1;
        volume(:,:,:,2) = im2;
        if exist('output','dir')~=7
            mkdir('output');
        end
        frame2gif(volume,fullfile('output_flow',[example '_input.gif']));
        volume(:,:,:,2) = warpI2;
        frame2gif(volume,fullfile('output_flow',[example '_warp.gif']));

        % visualize flow field
        clear flow;
        flow(:,:,1) = vx;
        flow(:,:,2) = vy;
        imflow = flowToColor(flow);
        clear r_flow;        
        %%%%% create backward flow 
        
        r_flow(:,:,1) = -vx;
        r_flow(:,:,2) = -vy;
        for h=1:c_height
            for w = 1: c_width
                h1 = round(h+flow(h,w,2));
                w1 = round(w+flow(h,w,1));
                if( h1 > c_height) h1 =c_height; end
                if( h1 <= 0) h1 = 1; end
                if( w1 > c_width) w1 = c_width; end
                if( w1 <= 0) w1 = 1; end
                r_flow(h1,w1,1) = -flow(h,w,1);
                r_flow(h1,w1,2) = -flow(h,w,2);
            end
        end       
        r_imflow = flowToColor(r_flow);
        
        figure(3);imshow(imflow);
        %pause
        figure(3);imshow(r_imflow);
        
        imwrite(imflow,fullfile('output_flow',[example '_flow.jpg']),'quality',100);
        imwrite(r_imflow,fullfile('output_flow',[example '_r_flow.jpg']),'quality',100);
        
        max_num_seg = 6;
        
        clear eImage;
        clear r_eImage;    
        
        if ( frame_pair == frame_start)
            if(1 == extend)
                eImage = imflow;
                eImage(:,:,3) = I1;  

                r_eImage = r_imflow;
                r_eImage(:,:,3) = I2;  

                %figure(5000), imshow(eImage), title('test');
                [segmentations, stabilities]=stable_segmentms1(eImage, max_num_seg,1);
                disp('draft segmentation on forward optic flow done');
                [r_segmentations, r_stabilities]=stable_segmentms1(r_eImage, max_num_seg,1);
                disp('draft segmentation on backward optic flow done');
            else
    %             eImage(:,:,1:2) = flow;
    %             eImage(:,:,3) = 0;
                eImage = imflow;

                r_eImage = r_imflow;
                [segmentations, stabilities]=stable_segmentms1(eImage, max_num_seg,0);
                disp('draft segmentation on forward optic flow done');
                [r_segmentations, r_stabilities]=stable_segmentms1(r_eImage, max_num_seg,0);
                disp('draft segmentation on backward optic flow done');
            end
            num_segments = size(segmentations,2);
            seg_size = zeros(num_segments);
        elseif ( frame_pair == frame_start+1)
            num_segments = size(r_segmentations,2); 
        else
            %%%% use segmentation from backward warping of the last frame-pair
            num_segments = num_segments;
        end
        

        
        for seg = 1:num_segments
            if ( frame_pair == frame_start)
                SI = cell2mat(segmentations(seg));
            elseif ( frame_pair == frame_start+1)
                SI = cell2mat(r_segmentations(seg));
            else
                %%% intitialize with segmentation from backward warping in
                %%% the last frame-pair
                SI = SI;
            end
            figure; imagesc(SI)
            
            if ( frame_pair <= frame_start+1)
                fprintf('\nstability of seg# %d = %f\n', seg+1, stabilities(seg));
            end
        
            [C, k1, k2] = unique(SI);
            num_centers = size(C,1);    
            clear seg_size;
            
            seg_size = zeros(num_segments);
            %%%% affine estimation over the detected regions
            for w = 1:num_centers               
                [affine_mat,seg_warp, curr_seg_size] = region_affine(im1,I1, I2, SI,w, flow, seg, num_segments,0, 1);
                Warps(w,:) = affine_mat(:);
                
                if (seg == num_segments)
                    seg_size(w) = curr_seg_size;
                    prim_warp_seg = genvarname(sprintf('prim_warp_seg%d',w));
                    eval([prim_warp_seg '= seg_warp;']);                
                    clear seg_warp;
                end
            end
%           Warps
                        
            if (seg == num_segments)
                %%%%%% create affine divergence matrix
                affine_affinity_mat = zeros(num_centers+1, num_centers+1); %% use additional dimension for better visualization
                %affine_affinity_mat = zeros(num_centers, num_centers);
                
                for m=1:num_centers
                    seg_size(m)
                    warp_seg_input = eval(sprintf('prim_warp_seg%d',m));
                    figure(100); imshow(warp_seg_input);
                    A = Warps(m,:);
                    A = reshape(A,3,3);
                    A(:,3)= [0 0 1];
                    tform = maketform('affine',A);
                    J = imtransform(warp_seg_input,tform);
                    %figure(300), imshow(J), title('apply corresponding affine motion to region');
                    
                    tinv  = inv(A);
                    tinv(:,3) = [0,0,1];
                    tform1 = maketform('affine',tinv);
                    J1 = imtransform(warp_seg_input,tform1); 
                    
                    [h w c] = size(J1);                    
                    sw = imresize(J1, [c_height c_width]);                    
                    %figure(400), imshow(sw), title('transform orig to canonical');
                    
                    for p=1:num_centers
                        if (p ~= m)
                            test_reg = eval(sprintf('prim_warp_seg%d',p));
                            figure(101); imshow(test_reg);
                            
                            A = Warps(p,:);
                            A = reshape(A,3,3);
                            A(:,3)= [0 0 1];
                            tform = maketform('affine',A);
                            J = imtransform(warp_seg_input,tform);
                            %figure(500), imshow(J), title('apply other affine motion to region');
                            
                            tinv  = inv(A);
                            tinv(:,3) = [0,0,1];
                            tform1 = maketform('affine',tinv);
                            %J1 = imtransform(warp_seg_input,tform1); 
                            J1 = imtransform(J,tform1);
                            
                            sw1 = imresize(J1, [c_height c_width]); 
                            
                            %figure(600), imshow(sw), title('transform other to canonical');
                            %Z = imsubtract(sw,sw1);
                            Z = abs(sw-sw1);
                            figure(700), imshow(imsubtract(sw,sw1)), title('take difference between two canonicals');
                            aff = sum(sum(Z(:,:,1)+Z(:,:,2)+Z(:,:,3)));
                            if (seg_size(m) ~= 0 )
                                affine_affinity_mat(m,p)= aff/(seg_size(m));
                            else                                                             
                                affine_affinity_mat(m,p) = aff;
                            end
                        else
                            affine_affinity_mat(m,p) = 0;
                        end
                        
                        %pause
                        
                    end
                end
                affine_affinity_mat(1:num_centers, 1:num_centers)
                %affine_affinity_gray = mat2gray(affine_affinity_mat);
                affine_affinity_mat(num_centers+1,:) = 0;
                affine_affinity_mat(:,num_centers+1) = 0;
                %figure(200); imagesc(affine_affinity_gray);
                figure(200); pcolor(affine_affinity_mat');
                
                
                
                %%% start post-processing for joining regions
                tr1 = affine_affinity_mat';
                %sym_aff_mat = mean2(affine_affinity_mat+ tr1);
                sym_aff_mat = max(affine_affinity_mat,tr1);
                %sym_aff_mat = triu(max(affine_affinity_mat,tr1),1);
                %sym_aff_mat(1:num_centers,1:num_centers);
                figure(300); pcolor(sym_aff_mat');

                for hier = 1:5
                    sprintf('h-level = %d\n', hier);
                    flag = zeros(1,num_centers);
                    k = num_centers;
                    k1 = 0;
                    for m=1:num_centers 
                        warp_seg_input = eval(sprintf('prim_warp_seg%d',m));
                        prim_warp_seg_h = genvarname(sprintf('prim_warp_seg%d_%d',m,hier));
                        eval([prim_warp_seg_h '= warp_seg_input;']);
                        
                        for p=m+1:num_centers
                            warp_seg_input1 = eval(sprintf('prim_warp_seg%d',p));
                            if ((flag(1,m) ~= 1)&& (flag(1,p) ~= 1) && (sym_aff_mat(m,p)< 0.1*hier) )
                                warp_seg = eval(prim_warp_seg_h) + warp_seg_input1;
                                %warp_seg_input = warp_seg_input +warp_seg_input1;
                                prim_warp_seg_h = genvarname(sprintf('prim_warp_seg%d_%d',m,hier));
                                eval([prim_warp_seg_h '= warp_seg;']); 
                                k = k-1;
                                flag(1,p) = 1;
                            end
                        end
                        
                        %%% estimate/refine affine parameters of
                        %%% the merged regions; and save in
                        %%% new_Warps
                        if ( max(flag(1,:)) > 0 )
                            [affine_mat,warp_seg, curr_seg_size] = region_affine(warp_seg,I1, I2, SI,w, flow, seg, num_segments,1,1);
                        end
                        
                        if (flag(1,m) == 1 )
                            continue;
                        end
                        k1 = k1+1;
                        new_Warps(k1,:) = affine_mat(:);
                        
                    end
                    h_num_seg = k;                    
                    sprintf('segments merged at h-level =%d\n',hier);                    
                    k1 = 0;
                    for i1 = 1:num_centers
                        figure(350);
                        if (flag(1,i1) == 1 )
                            continue;
                        end
                        k1 = k1+1;                        
                        subplot(1,h_num_seg,k1), subimage(eval(sprintf('prim_warp_seg%d_%d',i1,hier))), hold on
                    end
                    %pause
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    k1 = 0;
                    %%%%%%% debug
                    new_Warps(1,:) = Warps(1,:);
                    for i1 = 1:num_centers
                        figure(350);
                        if (flag(1,i1) == 1 )
                            continue;
                        end
                        k1 = k1+1;
                        
                        %% new set of warps for later MRF smoothing
                        new_Warps(k1,:) = Warps(i1,:);
                        
                        subplot(1,h_num_seg,k1), subimage(eval(sprintf('prim_warp_seg%d_%d',i1,hier))), hold on
                    end
                    new_num_centers = k1;
                    
                    %%%%% better show with semi-transparent segments in
                    %%%%% single image
                    iter = 0;
                    %test_im = im1;
                    for i = 1:num_centers
                        if (flag(1,i) == 1 )
                            continue;
                        end
                        test_im = eval(sprintf('prim_warp_seg%d_%d',i,hier));
                       
                        ind = find(test_im);
                        [A, B, C] = ind2sub(size(test_im), ind);
                        %size(test_im);
                        cl = mod(i,3)+1;                             
                        for j=1:size(A,1)       
                            x1 = B(j);
                            y1 = A(j); 
                            test_im(y1,x1,cl) = i*20;
                            %%%%%%
                            final_image(y1,x1,:) =test_im(y1,x1,:);
                            final_image(y1,x1,cl)= i*20;
                        end
                        iter = iter+1;
                    end
                    figure(2702),imshow(final_image),title('segments in single image');
                    %pause
                    
%                     figure(2702),
%                     subplot(1,num_centers,i), subimage(final_image), hold on
                end

%                 THRESHOLD = 0.28;
%                 [eigVectors,eigValues] = eig(sym_aff_mat(1:num_centers, 1:num_centers));
%                 % first eigenvector is the one which has the highest eigenvalue
%                 sz = size(eigVectors);
%                 firstEig = eigVectors(:,sz(1,2));
% 
%                 % plot the eigen vector corresponding to the largest eigen value
%                 [xx1,yy1,val1] = find(firstEig > THRESHOLD);
%                 [xx2,yy2,val2] = find(firstEig <= THRESHOLD);          
%                 figure,
%                 hold on;
%                 plot(data(xx1,1),data(xx1,2),'g*')
%                 plot(data(xx2,1),data(xx2,2),'b*')
%                 hold off;
%                 title('Clustering Results after thresholding the First Eigenvector');
%                 grid on;shg


                clear affine_affinity_mat;
            end
            
            
            
            %%%% show final segments
                       
            new_Warps = Warps;
            %num_centers = new_num_centers;
            
            starttime=cputime;
            %generate naive assignments (reconstruction error) 
%             [Wims,WimsA,WarpsInv]=warpimages(im1,im2,Warps);
            [Wims,WimsA,WarpsInv]=warpimages(im1,im2,new_Warps);
            [Assign1,Assign2,DiffMtx,DiffaMtx]=assignwarps(im1,im2,Wims,WimsA,example,frame_pair);
            
            W1=bvz_prep(I1,I2,DiffMtx,num_centers);
            output1=bvz_run(W1,I2,r,lambda);
            SmoothAssign1=reshape(output1,size(I2));
            W2=bvz_prep(I2,I2,DiffaMtx,num_centers);
            output2=bvz_run(W2,I1,r,lambda);
            SmoothAssign2=reshape(output2,size(I1));  

            %display smoothed assignments
            [SS1, SI] = displaywarps(im2,SmoothAssign1,num_centers,6, SS1);
            [SS2, SI] = displaywarps(im1,SmoothAssign2,num_centers,7, SS2);
 
%             % %Intersection masking for occlusion
%             %[M1,M2]=get_int_mask(SmoothAssign2,SmoothAssign1,Warps,WarpsInv,num_centers);
%             % 
%             % %display final assignments
%             % MF1 = displaywarps(im2,M1,num_centers,8, MF1);
%             % MF2 = displaywarps(im1,M2,num_centers,9, MF2);
%             %   
%             endtime=cputime;
%             fprintf('MRF smoothing finished in  %d seconds\n',endtime-starttime);
%             % 
%             disp('save plots and press any key to continue');
%             pause

        end
        fprintf('frame-pair #%d done\n', frame_pair);
        
        %pause
    end %loop frame-pair

    close(movie_in);
    close(NS1);
    close(NS2);
    close(SS1);
    close(SS2);
end


function[Wims,WimsA,WarpsInv]=warpimages(I1,I2,Warps)
    % warp the images
    [nwarps blah]=size(Warps);
    for count=1:nwarps
      Attemp=reshape(Warps(count,:),3,3);
      Attempinv=inv(Attemp);
      WarpsInv(count,:)=Attempinv(:)';
    end
    [I1rows,I1cols, clr]=size(I1);
    Wims=zeros(I1rows,I1cols,clr,nwarps);
    WimsA=zeros(I1rows,I1cols,clr,nwarps);
    tformtype='affine';
    for count=1:nwarps
        affine_mat = reshape(Warps(count,:),3,3);
        affine_mat1 = reshape(WarpsInv(count,:),3,3);            
        affine_mat(:,3) = [0,0,1];
        affine_mat1(:,3) = [0,0,1];            
        Wims(:,:,:,count)=apply_tform(I1,affine_mat',tformtype);
        WimsA(:,:,:,count)=apply_tform(I2,affine_mat1',tformtype); 
        
        %Recover the original image by transforming the distorted image.
        htrans1 = vision.GeometricTransformer(...
                   'TransformMatrixSource', 'Property', ...
                   'TransformMatrix',affine_mat,...
                   'OutputImagePositionSource', 'Property',...
                   'OutputImagePosition', [0 0 I1cols I1rows]);
        warped = step(htrans1,im2single(I1));
        release(htrans1);    

        htrans1 = vision.GeometricTransformer(...
               'TransformMatrixSource', 'Property', ...
               'TransformMatrix',affine_mat1,...
               'OutputImagePositionSource', 'Property',...
               'OutputImagePosition', [0 0 I1cols I1rows]);
        recovered = step(htrans1, im2single(I2));
        release(htrans1);

%         figure(24),imshow(I2-warped), title('I2 - warped');
%         figure(25),imshow(I1-recovered), title('I1 - recovered');

        %Compare recovered image to the original.
%         figure(20), imshow(I1)
%         title('base')
%         figure(21), imshow(recovered)
%         title('recovered')
% 
%         Compare warped image to the distorted.
%         figure(22), imshow(I2)
%         title('input')
%         figure(23), imshow(warped)
%         title('warped')
%         pause        
 
    end
end

%generate naive assignments (reconstruction error) 
% function[Assignment1,Assignment2,DiffMtx,DiffaMtx]=assignwarps(I1,I2,Wims,WimsA, Label)  
%     [rows,cols,clr,nwarps]=size(Wims);
%     
%     DiffMtx=zeros(rows,cols,nwarps);
%     DiffaMtx=zeros(rows,cols,nwarps);   
%     
%     color_error = 0;
%     for count=1:nwarps       
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       if (0 == color_error)
%           %%%%% warping error for Intensity only 
%           D1 = rgb2gray(I2);
%           D2 = rgb2gray(Wims(:,:,:,count));
%           D3 = rgb2gray(I1);
%           D4 = rgb2gray(WimsA(:,:,:,count));
%           diff = abs(D1-D2);
%           diffa = abs(D3-D4);
%       else
%           %%%% Warping error in RGB %%%%%%%%%%%%%%%%%%%%%
%           Diff = (I2-Wims(:,:,:,count));
%           Diffa =(I1-WimsA(:,:,:,count));       
%           diff(:,:) = Diff(:,:,1)+ Diff(:,:,2)+Diff(:,:,3);
%           diffa(:,:) = Diffa(:,:,1)+ Diffa(:,:,2)+Diffa(:,:,3);
%       end
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%       DiffMtx(:,:,count)= diff;
%       DiffaMtx(:,:,count)=diffa;      
%       ErrArray1(:,count)= diff(:);
%       ErrArray2(:,count)= diffa(:);
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     [Min,ArgMin]=min(ErrArray1,[],2);
%     asg1=reshape(ArgMin,[rows cols]);
%     [Min,ArgMin]=min(ErrArray2,[],2);
%     asg2=reshape(ArgMin,[rows cols]);
%     
%     Assignment1= asg1;
%     Assignment2= asg2; 
% end
 
%generate naive assignments (reconstruction error) 
function[Assignment1,Assignment2,DiffMtx,DiffaMtx]=assignwarps(I1,I2,Wims,WimsA,example,frame_pair) 
    %[rows, cols, nSegs] =size(BM1);    
    [rows,cols,clr,nwarps]=size(Wims);
    
    DiffMtx=zeros(rows,cols,nwarps);
    DiffaMtx=zeros(rows,cols,nwarps);   
    
    color_error = 1;
    seg_level = 0;
       
    for count=1:nwarps       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (0 == color_error)
          %%%%% warping error for Intensity only 
          D1 = rgb2gray(I2);
          D2 = rgb2gray(Wims(:,:,:,count));
          D3 = rgb2gray(I1);
          D4 = rgb2gray(WimsA(:,:,:,count));
          diff = abs(D1-D2);
          diffa = abs(D3-D4);
      else
          %%%% Warping error in RGB %%%%%%%%%%%%%%%%%%%%%
          Diff = (I2-Wims(:,:,:,count));
          Diffa =(I1-WimsA(:,:,:,count));       
          diff(:,:) = Diff(:,:,1)+ Diff(:,:,2)+Diff(:,:,3);
          diff(:,:) = diff(:,:)/3;
          diffa(:,:) = Diffa(:,:,1)+ Diffa(:,:,2)+Diffa(:,:,3);
          diffa(:,:) = diffa(:,:)/3;
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      DiffMtx(:,:,count)= diff;
      DiffaMtx(:,:,count)=diffa;      
      ErrArray1(:,count)= diff(:);
      ErrArray2(:,count)= diffa(:);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    end
    
    [Min,ArgMin]=min(ErrArray1,[],2);
    asg1=reshape(ArgMin,[rows cols]);
    [Min,ArgMin]=min(ErrArray2,[],2);
    asg2=reshape(ArgMin,[rows cols]);

    figure(400), imagesc(asg1)
           
    if (1 == seg_level)
%         diff(:,:) = diff(:,:)/3;
%         diffa(:,:) = diffa(:,:)/3;
        %%% start propagating from superpixel level; read the segmentations
        fprintf('reading superpixel segmentation\n');
        img_seg_file = sprintf('\\seg\\0_00%d.ppm', frame_pair); 
        seg_input = im2double(imread([example img_seg_file]));
        figure(600), imshow(seg_input);title('super-pixel')
        
        [C, m, n] = unique(reshape(seg_input, [], 3), 'rows');
        color_counts = accumarray(n, 1);
        nSegs = size(C,1);
        
        B = zeros(rows,cols);
        %%%%%%% assign at segment level %%%
        ind1 = zeros(1,nwarps);
        %ind2 = zeros(1,nwarps);
        %% first create 2D array from error vector
        err2D1 = reshape(ErrArray1, [rows cols nwarps]);
        err2D1 = err2D1/3;
        %err2D2 = reshape(ErrArray2, [rows cols nwarps]); 
        %err2D2 = err2D2/3;
                
        for m = 1:nSegs 
            %% create 1D (dim = nWarps) array of segment-level reconstruction error for each segment
            %% select ArgMin for best Label           
            [max_count, idx] = max(color_counts);
            color_counts(idx) = 0; 
            bw1 = (n == idx);
            bw1 = reshape(bw1, size(seg_input, 1), size(seg_input, 2));  
           
            %bw2 = double(BM2(:,:,m)); 
            elems = sum(sum(bw1));
            for w = 1:nwarps
                seg_err_mat = ((err2D1(:,:,w)).*bw1); %ind1(w) = sum(sum(abs(err2D1(:,:,w).*bw1)));                 
                S = seg_err_mat(:);
                T = S(~isnan(S));
                ind1(w) = sum(T)/elems;  %mean(nonzeros(T));
                                               
                %fprintf('error for warp# %d => sum = %f, max = %f, min =%f, mean val = %f, mean over abs = %f\n',  w, sum(T),max(T), min(T), sum(T)/elems, sum(abs(T))/elems);
                %figure(75), imshow(bw1)                
                %pause                  
                %%%%% set values for DiffMtx for each warp for each segment (not at pixel level)
                B(bw1 == 1) = ind1(w);
                DiffMtx(:,:,w) = B;
%               B(bw2 == 1) = ind2(w);
%               DiffaMtx(:,:,w) = B;
            end
            [Min, argmin1] = min(ind1);
            %[Min, argmin2] = min(ind2);                                  
            asg1(bw1==1) = argmin1;
            %asg2(bw2==1) = argmin2; 
        end
        figure(150), imshow(err2D1(:,:,1))
        figure(151), imshow(err2D1(:,:,2))
        
        if (1 == seg_level)        
            figure(250), imshow(DiffMtx(:,:,1));
            figure(255), imshow(DiffMtx(:,:,2));
            figure(500),imagesc(asg1);
        end
    end
    clear B;
    
    Assignment1 = asg1;
    Assignment2 = asg2;
 
 end
    

function [F,SI] = displaywarps(I,input,nwarps,fnum, F)
    % here nwarps does not include the static layer
    if(nnz(input)~=prod(size(input)))
      input=input+1;
    end
    for l = 1:3
        Assignment1(:,:,l)= input;
    end
    SI = I;
    for count=1:nwarps
      eval(['Mask' int2str(count) '=Assignment1==' int2str(count) ';']);
      figure(fnum)
      %eval(['subplot(1,nwarps,' int2str(count) '),imj((I).*Mask' int2str(count) ',[0,1]);']);
      eval(['subplot(1,nwarps,' int2str(count) '),imshow((I).*Mask' int2str(count) ',[0,1]);']);
      
%      eval([SI ' = ' SI ' + ' int2str(count) '.*Mask' int2str(count) ';']);
    end
%     frame = getframe;
%     writeVideo(F,frame);
end


function affine_mat = make_affine(input_points, base_points, base, input, gte )
    % Compute the transformation from the distorted to the original image.
    [tform_matrix inlierIdx] = step(gte, input_points, ...
        base_points);
    
    %Step 5: Compute Scale and Angle
    tform_matrix = cat(2,tform_matrix,[0 0 1]'); % pad the matrix
    tinv  = inv(tform_matrix);
    
    affine_mat = tinv; %tinv;
end

function [affine_mat, seg_warp, curr_seg_size] = region_affine(im1, I1, I2, SI, w, flow, seg, num_segments, def_reg, hh)
    if (def_reg == 0)
        seg_warp = zeros(size(im1));
        [A, B, V] = find(SI == w);          

        curr_seg_size = size(A,1);
    else
        %%%% if (def_reg == 1)
        seg_warp = im1;
        ind = find(im1);
        [A, B, C] = ind2sub(size(im1), ind);
        curr_seg_size = size(A,1);
    end
    
    for j=1:size(A,1)
        x1 = B(j);
        y1 = A(j);            
        index_x(j) = x1;
        index_y(j) = y1;      
        v_x(j) = flow(y1,x1,1);
        v_y(j) = flow(y1,x1,2); 

        if (def_reg == 0)
        %if (hh == 1)
            if (seg == num_segments)
                seg_warp(y1,x1,1:3) = im1(y1,x1,1:3);  
            end
        %end
        end
            
    end % Storing the extracted image in an array 

        base_points(:,1)= index_x;
        base_points(:,2)= index_y;

        input_points(:,1) = index_x + v_x;
        input_points(:,2) = index_y + v_y; 
        
        gte = vision.GeometricTransformEstimator; % defaults to RANSAC
        %gte.ExcludeOutliers = 0;
        gte.Transform = 'Affine';
        gte.NumRandomSamplingsMethod = 'Desired confidence';
        gte.MaximumRandomSamples = 1000;
        gte.DesiredConfidence = 99.8;
        %gte.Method = 'Least Median of Squares';

        %% make affine_estimate from correspondence
        affine_mat = make_affine(input_points, base_points,I1,I2, gte); 
        
        release(gte);  

        clear index_x;
        clear index_y;
        clear v_x;
        clear v_y;
        clear input_points;
        clear base_points;
end


 


