function [nuclei_model, Process_Recover_Filter] = load_models_and_filters();

%example:
%[nuclei_model,  process_recover_filter] = load_models_and_filters();


%Process_Recover_Filter
F=makeRFSfilters_sci;%creates oriented Gaussian zero-mean filters with different standard deviations.
Process_Recover_Filter=-F(:,:,31:40); % the 10 oriented Gaussian zero-mean filters of size 17x17 with sigma_x=2 and sigma_y=5

% %Denoise_Filter
% k=23;
% Denoise_Filter = fspecial('gaussian', k, (k-1)/6);

%Model for nuclei detection
 nuclei_model= cellpose(Model="nuclei");

end