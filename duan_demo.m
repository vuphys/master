function duan_demo

    alpha1=0.1;
    alpha0=0.3;
    theta1=5;
    theta2=5;
    iter=1000;

	%y  = double(imread('parrot.png'))/255;   % Initial image
	
    % load data
    ffcdata=load('GroundTruth');
    GroTru=ffcdata.data(:,:,1,1);
	
    % Adding noise
    rng(0); %reproducibility
	noise_img = GroTru+randn(size(GroTru))*0.1; % white Gaussian noise added to the image
	
    %call TGV
	denoise_img = duan_tgv(noise_img,alpha1,alpha0,theta1,theta2,iter);
	%imwrite(f,'noisy.png');
	%imwrite(x1,'TGVdenoised.png');
    
    %PSNR
    psnr_noisy = psnr(noise_img,GroTru);
    psnr_tgv = psnr(denoise_img,GroTru);

%SSIM
[mssim, ssim_map] = ssim(denoise_img, GroTru);


% Displaying results

figure(1), imshow( [GroTru, noise_img, denoise_img] )
title( sprintf('From the left,  original,  noisy %.6fdB,  TGV  %.6fdB and SSIM : %.8f', psnr_noisy, psnr_tgv, mssim ) );
toc
end