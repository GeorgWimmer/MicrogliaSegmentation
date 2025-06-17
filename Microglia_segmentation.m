function [Binary_Segmentation,Visualize_Segmentation] = Microglia_segmentation(im, nuclei_model,  ...
    denoise_factor, process_recover_filter)

%Examplar application
% [nuclei_model,  process_recover_filter] = load_models_and_filters();
% im_path='/home/pma/gwimmer/texture/mfiles/SCI/App/example_image.png';
% denoise_factor=15;
% [Binary_Segmentation,Visualize_Segmentation] = Microglia_segmentation(im_path, nuclei_model, ...
%  denoise_factor, process_recover_filter);

% [nuclei_model,  process_recover_filter] = load_models_and_filters();
% im_path='/scratch3/gwimmer/Bilder/AIBIA/spinal_cord_injury/data_gt_new/14dpi/SCI-Blue/SCI-blue_GM/SCI-blue_GM_2mmcaudal/E18068950GMDAPIchNfhrbIba1gpGFAP10-1120xstack_combined_noprocessing.png';
% denoise_factor=15;
%  [Binary_Segmentation,Visualize_Segmentation] = Microglia_segmentation(im_path, nuclei_model, ...
%   denoise_factor, process_recover_filter);




%%%%%%%%%%%%%%%%thresholds and other parameters


threshold_edge_brigthness_high=50;
threshold_edge_brigthness_low=40;
conn=4;   %4-connectivity for morphological operations (alternatively 8-connectivity)
%parameters for Nuclei detection
AverageCellDiameter=30;
Cell_Threshold=-2;
Flow_Error_Threshold=1;
min_core_size=150;
filter_threshold_high=11;
filter_threshold_low=7;
sigma_denoise=5;


if ischar(im)
    %load image
    im=imread(im);
else %is already an image and not the path to the image
    im=uint8(im);
end


im=uint8(im);
im_core=im(:,:,3);
im_micro=im(:,:,1);
im_micro_original=imadjust(im_micro);
im_micro=imadjust(im_micro);
im_micro=medfilt2(im_micro,[5,5]);
im_micro=medfilt2(im_micro,[5,5]); %added: second time median filtering
im_micro=imadjust(im_micro);
im_micro=uint8(im_micro);  %undenoised microglia image
%%%%%%%%%%%%%%%%%%%%

    
%CNN-based method for nuclei detection
core_labels = segmentCells2D(nuclei_model,im_core,ImageCellDiameter=AverageCellDiameter, CellThreshold=Cell_Threshold, FlowErrorThreshold = Flow_Error_Threshold);
for i=1:max(core_labels(:))
    if sum(sum(core_labels==i))<min_core_size
        core_labels(core_labels==i)=0;
    end
end



%image denoising
im_micro_undenoised=im_micro;

img=double(im_micro);
%img2=img-imfilter(img, denoise_filter, 'symmetric');
img2=img -imgaussfilt(img,sigma_denoise);
img2(img2<0)=0;
img2(img2>50)=50;
NI=(power(50-img2,2))./power(img+1,1.5);%noise image
  
im_micro=uint8(img-denoise_factor*NI); 
  
 
T = adaptthresh(im_micro,0);
T=imgaussfilt(T,20); %smootthing the local threshold
bin_micro=imbinarize(im_micro,T);



% reduce brightness of maybe too big microglia structures
bright_blobs=bwareaopen(bin_micro,5000);
im_micro_bright=uint8(img-(denoise_factor*3)*NI);
im_micro(bright_blobs)=im_micro_bright(bright_blobs);


%%%RECOVER PROCESSES:  recover processes that could have been deleted
%%%by the denoising. For this we employ the Maximum response filtering
%%%with directional Gaussian filters.
for i=1:size(process_recover_filter,3)
    fim(:,:,i)=imfilter(im_micro_undenoised, process_recover_filter(:,:,i));
end
fmax=max(fim,[],3);
fmax(fmax<0)=0;
recovered_areas=hysthresh(fmax,filter_threshold_high, filter_threshold_low) & hysthresh(im_micro_undenoised,threshold_edge_brigthness_high,threshold_edge_brigthness_low) ;

im_micro=uint8(double(im_micro)+ 100*recovered_areas);%just add any value (here we just took 100) so that these areas are marked as microglia
%Final preprocessed image

bin_micro=imbinarize(im_micro,T); %Binarized image




%find nuclei that are inside cell bodies as seed points for microglia
%segmentation

SE_imopen=strel('disk',8);
CellSoma = imopen(bin_micro,SE_imopen);%morphological opening
core_candidates = bwareaopen(CellSoma,250);
somalabels=bwlabel(core_candidates);
core_candidates_labels=core_labels;
core_candidates_labels(core_candidates==0)=0;


 valid_cores=zeros(size(core_labels));
 valid_core_labels=zeros(size(core_labels));

%%%find  probable  candidates  for microglias cell bodies (clearly
%%%distinguished from the background)
probable_candidates=zeros(size(core_candidates));
SEdilate=strel('disk',5);
somas=bwlabel(core_candidates);
for i=2:max(max(somas))
    soma_area=somas==i;
    mean_soma_brightness=mean2(im_micro_undenoised(soma_area));
    surrounding=imdilate(soma_area, SEdilate) & ~ soma_area;
    mean_surrounding_brightness=mean2(im_micro_undenoised(surrounding));
    if mean_soma_brightness > 2*mean_surrounding_brightness
       probable_candidates(soma_area)=true;
    end
end



for c=1:max(max(core_labels))
    if  sum(sum(core_candidates_labels==c)) > 0.8*sum(sum(core_labels==c))%at least 80 percent of the area of a core must fulfil the conditions to be not filtered out 
            valid_cores(core_labels==c)=255; 
            valid_core_labels(core_labels==c)=c; 
            im_micro(core_labels==c)=255;
    
 
    elseif  sum(sum(core_candidates_labels==c)) > 0.4*sum(sum(core_labels==c)) ...%at least 30 percent of the area of a core must fulfil the conditions to be not filtered out 
       &  sum(sum(core_candidates_labels==c))>min_core_size ...  %and a minimum size of 150 pixels
       &  sum(sum(core_candidates_labels==c)) <= 0.8*sum(sum(core_labels==c)) % and if the core candidtaes are not already added

        if any(probable_candidates(core_candidates_labels==c))  %if the soma candidate is a likely soma candidate
            valid_cores(core_candidates_labels==c)=255;
            valid_core_labels(core_candidates_labels==c)=c; 
        elseif sum(sum(core_labels==c & bin_micro))> 0.8* sum(sum(core_labels==c))  %if at least 80% of the core area are part of the binarized microglia image
            valid_cores(core_labels==c)=255; 
            valid_core_labels(core_labels==c)=c; 
            % probable_cores(core_labels==c)=255; 
            % probable_core_labels(core_labels==c)=c; 
            im_micro(core_labels==c)=255;
       end
     

    elseif sum(sum(core_labels==c & imclose(bwareaopen(bin_micro,250), strel('disk',8))))> 0.9* sum(sum(core_labels==c)) & ...
              sum(sum(core_labels==c & hysthresh(im_micro_undenoised,threshold_edge_brigthness_high,threshold_edge_brigthness_low))) > 0.9* sum(sum(core_labels==c)) 
                % probable_cores(core_labels==c)=255; 
                % probable_core_labels(core_labels==c)=c; 
                valid_cores(core_labels==c)=255; 
                valid_core_labels(core_labels==c)=c; 
                im_micro(core_labels==c)=255;
    end
  

end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binarize images and remove blobs that are not connected to a nucleus 


bwth2 = imbinarize(im_micro,T);




bblob=bwlabel(bwth2, conn);
h=(bblob>0 & valid_cores>0);
connectedBlobs=bblob(h>0);
blobclusters=unique(connectedBlobs);
blobadd=zeros(size(im_micro));
for i=1:length( blobclusters)
    blobadd(bblob==blobclusters(i))=1;
end
Binary_microglia= blobadd>0 ;
 

Binary_microglia=bwareaopen(Binary_microglia,200,conn); 

%%%%%%%%%%%%%%%%%%%%
%save the segmented images 

BW2 = bwperim(Binary_microglia);
BW2=imdilate(BW2,strel('disk',1));

BW3=zeros(size(Binary_microglia));
for i=1:max(valid_core_labels(:))
    BW3=BW3+bwperim(valid_core_labels==i);
end


BW3=imdilate(BW3,strel('disk',1));
%BW3_safe=imdilate(BW3_safe,strel('disk',1));
BW4_separated=BW2+BW3;

Visualize_Segmentation=uint8(im_micro_original);
Visualize_Segmentation(:,:,3)=im_core;
Visualize_Segmentation(:,:,2)=uint8(255*BW4_separated);%visualization of the microglia and nucleus borders

Binary_Segmentation=uint8(Binary_microglia)*255;
Binary_Segmentation(:,:,3)=uint8(valid_cores>0)*255;






end
