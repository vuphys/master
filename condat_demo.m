% Denoising/smoothing a given image y with the second order 
% total generalized variation (TGV), defined in 
% K. Bredies, K. Kunisch, and T. Pock, "Total generalized variation,"
% SIAM J. Imaging Sci., 3(3), 492-526, 2010.
%
% The iterative algorithm converges to the unique image x 
% (and the auxiliary vector field r) minimizing 
%
% ||x-y||_2^2/2 + lambda1.||r||_1,2 + lambda2.||J(Dx-r)||_1,Frobenius
%
% where D maps an image to its gradient field and J maps a vector 
% field to its Jacobian. For a large value of lambda2, the TGV 
% behaves like the TV. For a small value, it behaves like the 
% l1-Frobenius norm of the Hessian.
%		
% The over-relaxed Chambolle-Pock algorithm is described in
% L. Condat, "A primal-dual splitting method for convex optimization 
% involving Lipschitzian, proximable and linear composite terms", 
% J. Optimization Theory and Applications, vol. 158, no. 2, 
% pp. 460-479, 2013.
%
% Code written by Laurent Condat, CNRS research fellow in the
% Dept. of Images and Signals of GIPSA-lab, Univ. Grenoble Alpes, 
% Grenoble, France.
%
% Version 1.0, Oct. 12, 2016

% This code show a demo run for a set of parameters 
% based on Condat TGV code

function condat_demo
tic
	Nbiter= 5000;	% number of iterations
	lambda1 = 0.102021; 	% regularization parameter
	lambda2 = 4.512587;	% regularization parameter
  
	tau = 0.01;		% proximal parameter >0; influences the
                    % convergence speed

	%y  = double(imread('parrot.png'))/255;   % Initial image
	
    % load data
    ffcdata=load('GroundTruth');
    GroTru=ffcdata.data(:,:,1,1);
	
    % Adding noise
    rng(0); %reproducibility
	noise_img = GroTru+randn(size(GroTru))*0.1; % white Gaussian noise added to the image
	
    %call TGV
	denoise_img = condat_tgv(noise_img,lambda1,lambda2,tau,Nbiter);
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