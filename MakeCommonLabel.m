function [ output_args ] = MakeCommonLabel( keyF,Dir)
%INITIALALL Summary of this function goes here
%   Detailed explanation goes here
	step = 15;
	DirName = dir(Dir);
	Key = cell(numel(DirName),1);
	KeyFrameNum = numel(DirName);
	for i= 1:numel(DirName)
		path = DirName{i}.name
		tmp = load('Final_D');
		Key{i}.D = tmp.tmp;
		tmp = load('Final_Model');
		Key{i}.Model = tmp.tmp;
	end
	% load P
	
	
	
	% warp to now frame
	Dproposals = Dproposals(size(Key{1}.D,1),size(Key{1}.D,2),KeyFrameNum);
	Model = cell(KeyFrameNum,1);
	% hole stay -1 label and assign big cost(no information from this)
	Model.frame = k*step;
	[D Model] = Warp2CurrentView(Key,P,keyFrameNum,keyF,step);
	vals.visibility = false;
	[D info.InitialEach ObjModel] = obj_fuse_proposals(vals, Dproposals, options, info.segpln_Obj, Expan);
	delete(['Frame_' num2str(i) '_info.mat'],['Final_' num2str(i) '_D.mat'],['Final_' num2str(i) '_Model.mat']);
	save(['Frame_' num2str(i) '_info.mat'],'info');
	save(['Final_' num2str(i) '_D.mat'],'D');
	save(['Final_' num2str(i) '_Model.mat'],'ObjModel');
	% warp back to each View
	
end

