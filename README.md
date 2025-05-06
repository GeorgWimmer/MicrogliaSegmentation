# MicrogliaSegmentation

Here is the code for microglia segmentation from the paper 'Microglia cell segmentation using a hand-crafted method capable of handling high noise levels in image data'.

Requirements: Matlab (at least version R2023b (for the necessary cellpose library)) 

First load the Matlab code and the two example images (PNG images) down in a directory.
Start Matlab, and add the path to the directory ("addpath(genpath('path_to_directory'))")

The code can be applied to 
    (1) Tiff images including Iba-1 channels and DAPI-channels. (If your images are in 'czi' format, e.g. use the code in https://github.com/cgohlke/czifile/blob/master/  (czi2tif.py) to transform it to tiff format. (Example images could not be provided because they have over 200 MB)
    (2) PNG images, where the first channel is a maximum intensity projection over the Iba1-layers and the third channel is the maximum intensity projection over the DAPI layers. (2 exemplar png images are provided, which were processed from  tiff images).

To apply the matlab code to all images in a directory (dir_in) and write out the segmentation results in the directory 'dir_out', apply the following code in matlab:

    (1) Tiff images:
          [Binary_Segmentation,Visualize_Segmentation] = Segment_all_Images(dir_in, dir_out);
    
    (2) PNG images (2 PNG images are included in the code to try the provided code on them):
          [Binary_Segmentation,Visualize_Segmentation] = Segment_all_Images(dir_in, dir_out, 'ImageFormat','PNG');
    

See the description in the function  function 'Segment_all_Images.m' for further information .    
