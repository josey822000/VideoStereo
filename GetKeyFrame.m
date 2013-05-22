function [ output_args ] = GetKeyFrame( Dir )
%GETKEYFRAME Summary of this function goes here
%   Detailed explanation goes here
    DirName = dir(Dir);
    for i= 1:numel(DirName)
        path = DirName{i}.name
        tmp = load('Final_D');
        D = tmp.tmp;
        tmp = load('Final_info');
        info = tmp.tmp;
    end
end

