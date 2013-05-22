function [ DepthProposals] = DepthProposal_gen( info,d_min,d_max,pid)
%DEPTHPROPOSAL_GEN Summary of this function goes here
%   Detailed explanation goes here
segNum = max(max(info.segments(:,:,pid));
sz = size(info.disp);
DepthProposals = zeros(sz(1),sz(2),segNum);
for s = 1:segNum
	[X Y] = meshgrid(sz(2),sz(1));
	tmp = -(X * info.plane{pid}.N(s,1) + Y * info.plane{pid}.N(s,2) + info.plane{pid}.N(s,3));
	tmp(find(tmp<d_min || tmp>d_max)) = realmax;
	DepthProposals(:,:,s) = tmp;
end
end

