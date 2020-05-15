% This code show a demo run for a set of parameters 
% based on FFT TGV code

function fft_demo_gpu(a,b,c,d)
tic
% ADMM parameters
nite = gpuArray(100); % number of iterations
% balancing weights for Total Variation
alpha = gpuArray(a);  % 1st order
beta = gpuArray(b); % 2nd order


% Loading data   

ffcdata=load('GroundTruth');
GroTru=gpuArray(ffcdata.data(:,:,1,1));

%fname = 'parrotgray.png';
%I = im2double( imread( fname ) );
%I0 = I; % original as the reference

% Adding noise
rng(0); %reproducibility
noise_img = imnoise(GroTru,'gaussian',c,d);
%noise_img = GroTru+randn(size(GroTru))*0.3; % white Gaussian noise added to the image
	

%denoise_img = zeros( size(GroTru),'gpuArray' ); %Preallocation

% Calling TGV fucntion
%for c = 1:size(GroTru,3)
	denoise_img = arrayfun(@fft_tgv,noise_img, alpha, beta, nite );
%end


% PSNR
psnr_noisy = arrayfun(@psnr,noise_img,GroTru);
psnr_tgv = arrayfun(@psnr,denoise_img,GroTru);

% SSIM
[mssim, ssim_map] = arrayfun(@ssim,denoise_img,GroTru);

% MS SSIM
[ms_ssim, ms_ssim_map] = arrayfun(@multissim,denoise_img,GroTru);

% FSIM
fsim=arrayfun(@FeatureSIM,GroTru,denoise_img);

% ESSIM
essim=arrayfun(@ESSIM,GroTru,denoise_img);

%PSNR HVS
psnr_hvs=arrayfun(@psnrhvsm,GroTru,denoise_img);

% PSNR HMA
psnr_hma=arrayfun(@psnrhma,GroTru, denoise_img);

% SR SIM
srsim=arrayfun(@SR_SIM,GroTru, denoise_img);


% Displaying results

figure(1), imshow( [GroTru, noise_img, denoise_img] );
title( sprintf('From the left,  original,  noisy %.6fdB, TGV %.6fdB \n SSIM: %.8f FSIM: %.8f\n ESSIM: %.8f PSNR HVS: %.8f PSNR HMA: %.8f\n MS SSIM: %.8f SR SIM: %.8f', psnr_noisy, psnr_tgv, mssim, fsim,essim,psnr_hvs,psnr_hma, ms_ssim,srsim) );
toc
end