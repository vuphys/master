% This code show a demo run for a set of parameters 
% based on He TGV code

function he_demo

%blur kernel definition and selection
K     =   fspecial('average',1); %for denoising

%K     =   fspecial('Gaussian',9,3); % for debluring

OrgSigma=0.1; %noise standard deviation of the noise image

%fprintf('The noise std of the observed image: %g.\n', OrgSigma); 


% Generating Ground Truth image
ffcdata=load('GroundTruth');
GroTru=ffcdata.data(:,:,1,1);

% Adding noise:
N     =   numel(GroTru);             
[m,n] =   size(GroTru);
blur_im  = imfilter(GroTru,K,'circular','conv');

BSNR=20*log10(norm(blur_im(:)-mean(blur_im(:)),'fro')/sqrt(N)/OrgSigma);
%fprintf('BSNR of the observed image: %g dB.\n', BSNR);

rng(0); %reproducibility
noise_img = GroTru + randn(m,n)*0.1;
%noise_img        = blur_im + OrgSigma*randn(m,n);       
%PSNR_F    =psnr(noise_img,GroTru);
%fprintf('PSNR of the observed image: %g dB.\n\n', PSNR_F);


tao0=0.006; tao1=0.03;% slightly tuning may cause more appealing result 
if size(K)==1
tao  = -BSNR*tao1+1.09; %for denoising
else
tao  = -BSNR*tao0+1.09; %for deblurring
end
c    =  tao*m*n*OrgSigma.^2; % upper bound for the constraint

% Input parameters of this TGV model:
Param.OrigIm     = GroTru;      Param.MaxIter    = 10000; 
Param.SolRE      = 1e-6;    Param.UpBound    = c;
Param.Beta       = 4.556444;       Param.Gamma      = 1;
Param.Tao        = 1;       Param.BSNR       = BSNR;
Param.alpha0=0.684097;             Param.alpha1=0.032299;

% Call TGV function:
output = he_tgv(noise_img, K, Param); %% main program

% Output value:
denoise_img     = output.Sol;           Reglambda  = output.Reglambda;
PSNR       = output.PSNR;     mse        = output.MSE;
IterTime   = output.IterTime; Fvalue     = output.Fvalue;
fprintf('Proposed APE_TGV_ADMM: Elapsed time is %g seconds.\n', IterTime(end));
fprintf('Proposed APE_TGV _ADMM: Total iterations is %g.\n', length(IterTime));
fprintf('Proposed APE_TGV_ADMM: Final regularization parameter is %g.\n', Reglambda(end));
fprintf('Proposed APE_TGV_ADMM: Final PSNR is %g dB.\n\n\n', PSNR(end));


% PSNR:
psnr_noisy = psnr(noise_img,GroTru);
psnr_tgv = psnr(denoise_img,GroTru);

% SSIM:
mssim= ssim(denoise_img,GroTru);

% Displaying results

figure(1), imshow( [GroTru, noise_img, denoise_img] );
title( sprintf('From the left,  original,  noisy %.6fdB,  TGV %.6fdB, ssim: %.8f', psnr_noisy, psnr_tgv, mssim ) );

end