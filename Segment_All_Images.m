function [Binary_Segmentation,Visualize_Segmentation] = Segment_all_Images(Data_in, Data_out, varargin)
    
%   
%   Parameters are:
% 
%       'Data_in'             Directory in which are the images to be
%                           segmented
% 
%       'Data_out'            Directory where the segmented images are written out 
%   
%
%      'ImageFormat'   What image format is used ('Tiff' or 'PNG'). Either an Tiff image,
%                       where, we need to define which channels are
%                       containing microglia data (Iba1 channels) and which Nucleus data (DAPI channels).
%                        Or PNG images, which are already processed to have
%                        the microglia information (maximum intensity projection over the Iba1-channels) in the first and the
%                        Nucleus information (maximum intensity projection over the DAPI-channels)  in the third channel.
%                        The second image channel is empty (zeros)    
%
%     
%      'NucleusChannels'  What channels are Nucleus channels (only needed for tiff images). Needs a
%                           starting nummer and a stop number. E.g [1,7]
%                           means that the first seven channels of the Tiff image
%                           are Nucleus channels (Maximum intensity
%                           projection is carried out to reduce the information
%                           to one channel)
%       
%      'MicrogliaChannels'   What channels are microglia channels  (only needed for tiff images). Needs a
%                           starting nummer and a stop number.
%
%      'DenoiseFactor'   How strict is the denoising process. The higher the
%                          value, the more image parts are removed by the
%                          denoising process
%
% Outputs are the binary segmentation mask (microglia and their nuclei)
%
%%%%%Exemplar application
%% for tiff images
% Data_in='path_to_the_tiff_images'; %directory where the tiff images are that should be segmented
% Data_out='path_were_the_segmented_images_should_be_written_out';
% [Binary_Segmentation,Visualize_Segmentation] = Segment_All_Images(Data_in, Data_out)
%% or for png images  
%[Binary_Segmentation,Visualize_Segmentation] = Segment_All_Images(Data_in, Data_out, 'ImageFormat','PNG')






    p = inputParser;
    p.addRequired('Data_in',@isstr);
    p.addRequired('Data_out',@isstr);
    p.addParamValue('ImageFormat','Tiff', @(x)any(strcmpi(x,{'Tiff','PNG'})));
    p.addParamValue('NucleusChannels',[1,7], @(x) x(2)> x(1) && length(x)==2 );
    p.addParamValue('MicrogliaChannels',[8,14], @(x) x(2)> x(1) && length(x)==2 );
    p.addParamValue('DenoiseFactor',15,  @(x)x>=1 && x<= 50 );

    p.parse(Data_in, Data_out, varargin{:});
    ImageFormat = p.Results.ImageFormat;
    NucleusChannels = p.Results.NucleusChannels;
    MicrogliaChannels = p.Results.MicrogliaChannels;
    DenoiseFactor = p.Results.DenoiseFactor;

   


    switch ImageFormat
        case 'Tiff'
            ending='*.tif';
        case 'PNG'
            ending='*.png';
    end
  

    if ~ isdir(Data_out)
        mkdir(Data_out);
    end





    [nuclei_model, process_recover_filter] = load_models_and_filters();

    filelist=dir(fullfile(Data_in, ending));
    display(sprintf('%d images are to segment', length(filelist)));
    for im_counter=1:length(filelist)
        impath=fullfile(filelist(im_counter).folder, filelist(im_counter).name);
        display(sprintf('Image nr. %d is processed : %s', im_counter, filelist(im_counter).name));

        
        switch ImageFormat
            case 'Tiff'
                im = load_tiff_image(impath, MicrogliaChannels, NucleusChannels);
            case 'PNG'
                im=imread(impath);
        end

        [Binary_Segmentation,Visualize_Segmentation] = Microglia_segmentation(im, nuclei_model, ...
            DenoiseFactor, process_recover_filter);



        [pathstr,name,ext] = fileparts(impath);
        out_binary=fullfile(Data_out, [name,'_binary.png']);
        out_visualize_segmentation=fullfile(Data_out, [name,'_visualization.png']);
        imwrite(Binary_Segmentation, out_binary);
        imwrite(Visualize_Segmentation,out_visualize_segmentation);

    end




