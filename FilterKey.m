function [ output_args ] = FilterKey( vals, options, Key)
%INITIALALL Summary of this function goes here
%   Detailed explanation goes here
    
    % rec key frame 
    for i = 1:size(vals.P,3)
        Expan = false;
        Dproposals = Dproposals(:,:,1:7);
        tmpvals.visibility = false;
        %
        Expan = true;
		tmpvals.visibility = false;
        
        options.max_iters = 5;
        options.average_over = 5;
        options.imout = i;
        % here only two proposals used, fucking lack information, what is
        % the pros. to mark the black region? black region is occluded, the
        % info can be recoverd by another proposal
        close all;
		[D info.InitialEach ObjModel] = obj_fuse_proposals(tmpvals, D, options, Model, Expan);
		delete(['Frame_' num2str(i) '_info.mat'],['Final_' num2str(i) '_D.mat'],['Final_' num2str(i) '_Model.mat']);
		save(['Frame_' num2str(i) '_info.mat'],'info');
		save(['Final_' num2str(i) '_D.mat'],'D');
		save(['Final_' num2str(i) '_Model.mat'],'ObjModel');
			
	end
end

