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

% This code run and plot Monte Carlo simulation based on Condat code.
% There are five random parameters in this simulation which are
% lambda_1(alpha 1), lambda_2(alpha 0) for TGV and seed, mean, standard
% deviation for random Gaussian noise image

% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function  monte_sim_condat_all(M,N)
%% Preallocate:
H=zeros(N,9);   % preallocate storage for all value
T=zeros(N,1);   % preallocate storage for lambda_1
V=zeros(N,1);   % preallocate storage for lambda_2
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for MSE
Y=zeros(N,1);   % preallocate storage for Brisque
Z=zeros(N,1);   % preallocate storage for Niqe
R=zeros(N,1);   % preallocate storage for Piqe

denoise=cell(N,1); % create cell array for denoise images storage

result=struct();   % create struct for all results storage

%% TGV parameters:
Nbiter= 600;	% number of iterations
%lambda1 = 0.1; 	% regularization parameter
%lambda2 = 0.2;	% regularization parameter
tau = 0.01;		% proximal parameter >0; influences the
                %    convergence speed
%% Loading data   

ffcdata=load('GroundTruth');        %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);       %load a Ground Truth image

tic    
%% Adding noise
for k=1:M
    
    ans1=sprintf('Condat code create noise image %d/%d',k,M)
    
    rng(k,'twister')
    s=rand(1,1);
    t=s+k; %for reproducibility of the noise image
    rng(t,'twister')
    std=rand(1,1)*2; % create random standard deviation for Gaussian noise
                     % from [0,2]
    m=2.*rand(1,1)-1;   % create random mean for Gaussian noise
                        % from [-1,1]
    ans2=sprintf('seed: %.3f, mean: %.3f, standard deviation: %.3f',t,m,std)
    
    
    result.I(k).seed=s;     % store seed to struct
    result.I(k).mean=m;     % store mean to struct
    result.I(k).std=std;    % store standard deviation to struct
    

noise_img = imnoise(GroTru,'gaussian',m,std);   % adding Gaussian noise

result.I(k).n_img=noise_img; %store created noise image to struct
	
%% Monte Carlo simulation:

%count=0;
format long;

     parfor i=1:N % parallel computing
        

        % Random seed:
        rng(i,'twister') %for different seed in different stream
        lambda1=rand(1,1)*5;
        lambda2=rand(1,1)*5;
       
        % Store parameters
        T(i)=lambda1;
        V(i)=lambda2;
        
        %Call TGV:
        denoise_img = condat_tgv(noise_img,lambda1,lambda2,tau,Nbiter);
	
        %SSIM
        [mssim, ssim_map] = ssim(denoise_img, GroTru);
        [mulmssim,mulssim_map]=multissim(denoise_img,GroTru);
        
        %MSE
        err = immse(denoise_img, GroTru);
        
        %Brisque
        brisque_score = brisque(denoise_img);

        %PSNR
        psnr_tgv = psnr(denoise_img,GroTru);
        
        %Niqe
        niqe_score=niqe(denoise_img);
        
        %Piqe
        piqe_score=piqe(denoise_img);
        
        % Store MSSIM and PSNR:
        U(i)=mssim;
        P(i)=psnr_tgv;
        Q(i)=mulmssim;
        X(i)=err;
        Y(i)=brisque_score;
        Z(i)=niqe_score;
        R(i)=piqe_score;
        %count=count+1;

        % Displaying results
        ans4=sprintf('lambda1: %.3f, lambda2: %.3f and mssim: %.8f', lambda1, lambda2, mssim)
        
         denoise{i}=denoise_img; % save denoise image to cell array       
    end
    
    %% Save data:
    
    H(:,1)=T; %lambda1
    H(:,2)=V; %lambda2
    H(:,3)=U; %MSSIM
    H(:,4)=P; %PSNR
    H(:,5)=Q; %MS SSIM
    H(:,6)=X; %MSE
    H(:,7)=Y; %Brisque
    H(:,8)=Z; %Niqe
    H(:,9)=R; %Piqe
         
    result.I(k).factor=H;
    result.I(k).den_img=denoise; %save cell array to result struct

    ans5=sprintf('finish update H, continue to next noise level')
end
toc
    savefile=sprintf('DATA/condat_%d_%d_all.mat',M,N); % create file name
    save(savefile,'result'); % save struct result to file
    
    ans6=sprintf('finish Condate code')
    
end
