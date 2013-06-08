function MaxD = DeterminMaxD(Dir, sz)

dirpath = dir(Dir);
MaxD = 0;
MinD = 10000;
for i=3:size(dirpath,1)
    disp(dirpath(i).name);
	fp = fopen(fullfile(Dir,dirpath(i).name), 'rb');
	data = zeros(sz(2),sz(1),'single');
	data = fread(  fp, prod(sz), 'single');
	MaxD = max(data);
	MinD = min(data);
    data = reshape(data,sz(2),sz(1))';
    
	imshow(data/MaxD);
	pause;
end