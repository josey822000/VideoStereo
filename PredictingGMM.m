function [out out2] = PredictingGMM( samples,filename )
%PREDICTINGGMM Summary of this function goes here
%   Detailed explanation goes here
samples = reshape(permute(samples,[3 2 1]),3,[]);
[out out2] = PredictGMM(samples,filename);

end

