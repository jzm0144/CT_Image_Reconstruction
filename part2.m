clc; clear; load project_1_new_data.mat;

ksamps = kspace;
numProjections = size(ksamps, 2);
readout_points = size(ksamps, 1);

ktraj = getRadialTraj(readout_points , numProjections);


figure(1); colormap('gray'); plot(reshape(ktraj, 1, []));

dcf = calcdcflut(ktraj)';

dcf = dcf/max(dcf);


kmax = max(abs(ktraj));
if (kmax > 0.5)
	disp('Warning:  Maximum k-space > 0.5.  ');
	disp('		Conventionally k-space files should not exceed 0.5');
	disp(' 		This may result in problems later.');
	disp(' ');
	q = input('Normalize k-space to 0.5 max magnitude? (y/n)','s');
	if (q(1)=='y')
		disp('Normalizing...');
		kraj = 0.48*kraj/kmax;
	end;
end;

[gdat] = gridkb(reshape(ktraj, 1, []), ksamps, dcf, 362, 3, 2);

im = fftshift(fft2(fftshift(gdat)));

ax = figure;
cmap = [0:255].'*[1 1 1] / 256;
colormap(cmap);
imagesc(abs(im)/5120); title('Gridding');
saveas(ax, 'Gridding_w3.png');