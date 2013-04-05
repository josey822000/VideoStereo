% Load default options for the desired algorithm
algs = {'my'};
options = ojw_default_options(algs{1});

% In the case of stereo, download the desired sequence
% sequence = 'cones';
% options.input_type = 'stereo';
% download_stereo(sequence);
sequence = 'NewRoad';
download_img_sequence(sequence);
cd(sequence);
mkdir('GMM');
% Change any options here
options.dim_out = [0 0 450 375]; % Output image dimensions: [start_x-1 start_y-1 width height]
% cones 450 375
% tsukuba 384 288
% venus 434 383
options.imout = 1; % Index of projection matrix to use for output
% only account for two images?
options.nclosest = [1 2];
% options.nclosest = [1 7]; % Input images to use, in terms of distance of camera centres from output view
% cones [3 7]
% tsukuba [1 4]
% venus [3 7]
% Generate output matrices
Pout = zeros(0, 0, 1); % Use input matrix defined by options.im_out
% Pout = ojw_genview('steady', [1 63], (0:29)/29); % 30 frame steadicam sequence 
% Pout = ojw_genview('stereo', 1:63, 'l'); % Generate left views of stereo pairs

for a = 1:size(Pout, 3)
    % Call the rendering function
    [A out] = ibr_render(options, Pout(:,:,a));
	save('OutImg','A');
	save('outModel','out');
    % Save the output
    clf; sc(A);
end