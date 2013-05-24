% Load default options for the desired algorithm
% algs = {'InitialKey'};
% algs = {'InitialAll'};
algs = {'TubeFuse'};
options = ojw_default_options(algs{1});
close all;
% In the case of stereo, download the desired sequence
% sequence = 'cones';
% download_stereo(sequence);
sequence = 'Angkor';
options.sequence = sequence;
download_img_sequence(sequence);
if ~exist('data.mat')
cd(sequence);
end
mkdir('GMM');
% Change any options here
% cones 450 375
% tsukuba 384 288
% venus 434 383
options.imout = 1; % Index of projection matrix to use for output
% NewRoad test
options.nclosest = 1:105;
% options.nclosest = [15 1 30];
options.KeyFrame = [1 15 30 45 60 75 90 105];
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
% 	save('OutImg','A');
% 	save('outModel','out');
    % Save the output
    clf; sc(A);
end