% This code show a demo run for a set of parameters 
% based on UC TGV code

function uc_demo

	maxIterations= 10000;	% number of iterations
    checkIterations= 100;   % check step
	alpha0 = 0.008298; 	% regularization parameter
	alpha1 = 0.000093; % regularization parameter
	lambda = 0.000861;		%fidelity parameter

    % Loading data   

    ffcdata=load('GroundTruth');
    GroTru=ffcdata.data(:,:,1,1);

	%f  = double(imread('parrotgray.png'))/255;   % Initial image
	
    % Additional noise
    
    rng(0); %reproducibility
	noise_img = GroTru+randn(size(GroTru))*0.1;
	
	% Call TGV function:
	denoise_img = uc_tgv(noise_img, lambda, alpha0, alpha1, maxIterations, checkIterations);
	
    %imwrite(noise_img,'noisy.png');
	%imwrite(denoise_img,'TGVdenoised.png');
    
    % PSNR:

    psnr_noisy = psnr(noise_img,GroTru);
    psnr_tgv = psnr(denoise_img,GroTru);
    
    % SSIM:
    mssim=ssim(denoise_img,GroTru);

% Displaying results

figure(1), imshow( [GroTru, noise_img, denoise_img] );
title( sprintf('From the left,  original,  noisy %.6fdB,  TGV %.6fdB, SSIM %.8f', psnr_noisy, psnr_tgv,mssim ) );

end
