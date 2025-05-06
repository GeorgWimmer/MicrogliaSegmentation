function [im_combined] = load_tiff_image(source, microglia_channels, nucleus_channels)
%Load tiff_image and convert it to standard RGB image with microglia
%information in the first and nucleus information in the third channel. 

tiff_image=loadtiff(source);

max_value=max(tiff_image(:));
tiff_image=double(tiff_image) / double(max_value) *255;
tiff_image=uint8(tiff_image);

%max intensity projection of the microglia channels
im_microglia = tiff_image(:, :, microglia_channels(1):microglia_channels(2));
im_microglia = max(im_microglia, [], 3);

%max intensity projection of the nucleus channels
im_nucleus = tiff_image(:, :, nucleus_channels(1):nucleus_channels(2));
im_nucleus = max(im_nucleus,[],3);

im_combined=im_microglia;
im_combined(:,:,3)=im_nucleus;

end
