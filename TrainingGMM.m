function TrainingGMM(samples,filename)

samples = reshape(permute(samples,[3 2 1]),3,[]);
TrainGMM(samples,filename);
end