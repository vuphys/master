% This code run and plot Monte Carlo simulation based on UC code
% There are six random parameter in the simulation which are alpha 0,
% alpha 1, lambda for TGV, and seed, mean, standard deviation for random 
% Gaussian noise image

% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function monte_sim_uc_all(M,N)

%% Preallocate:
H=zeros(N,10);   % preallocate storage for all value
T=zeros(N,1);   % preallocate storage for alpha_1
V=zeros(N,1);   % preallocate storage for alpha_0
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for MSE
Y=zeros(N,1);   % preallocate storage for Brisque
Z=zeros(N,1);   % preallocate storage for Niqe
V1=zeros(N,1);  % preallocate storage for lambda
R=zeros(N,1);   % preallocate storage for Piqe

denoise=cell(N,1); % create cell array for denoise images storage

result=struct();    % create struct for all results storage
image=struct(); %create struct for all images storage


%% TGV parameters

	maxIterations= 10000;	% number of iterations
    checkIterations= 100;   % step of iteration check
	%alpha0 = 0.0009;       % regularization parameter
	%alpha1 = 0.0017;       % regularization parameter
	%lambda = 0.01;         % fidelity parameter

%% Loading data   

ffcdata=load('GroundTruth');        %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);       %load a Ground Truth image
	
%% Adding noise
tic
for k=1:M
    
    ans1=sprintf('UC code create noise image %d/%d',k,M)
    
    rng(k,'twister')
    s=round(rand(1,1),8);
    t=s+k; %for reproducibility of the noise image
    rng(t,'twister')
    std=round(rand(1,1)*2,8); % create random standard deviation for Gaussian noise
                     % from [0,2]
    m=round(2.*rand(1,1)-1,8);   % create random mean for Gaussian noise
                        % from [-1,1]
    ans2=sprintf('seed: %.3f, mean: %.3f, standard deviation: %.3f',t,m,std)
    
    
    result.I(k).seed=t;     % store seed to struct
    result.I(k).mean=m;     % store mean to struct
    result.I(k).std=std;    % store standard deviation to struct
    
    rng(t)
    noise_img = imnoise(GroTru,'gaussian',m,std);   % adding Gaussian noise
	
    image.I(k).n_img=noise_img; %store created noise image to struct
%% Monte Carlo simulation:

%count=0;
format long;
     parfor i=1:N % parallel computating
        
        % Random seed:
        rng(i,'twister') %for different seed in different stream
        alpha1=rand(1,1)/100;
        alpha0=rand(1,1)/100;
        
        lambda=0;
        while(lambda==0) %To make sure lambda is not zero
        lambda=rand(1,1)/10;
        end
        
        % Store parameters
        T(i)=alpha1;
        V(i)=alpha0;
        V1(i)=lambda;
        
        %Call TGV:
        denoise_img = uc_tgv(noise_img, lambda, alpha0, alpha1, maxIterations, checkIterations);	
        
        %SSIM
        [mssim, ssim_map] = ssim(denoise_img, GroTru);
        [mulmssim,mulssim_map]=multissim(denoise_img,GroTru);
        
        %MSE
        err = immse(denoise_img, GroTru);
        
        %Brisque
        score_brisque = brisque(denoise_img);

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
        Y(i)=score_brisque;
        Z(i)=niqe_score;
        R(i)=piqe_score;
        %count=count+1;

        % Dsiplaying results
        ans3=sprintf('alpha1: %.3f, alpha0: %.3f, lambda: %.3f and mssim: %.8f', alpha1, alpha0, lambda, mssim)
        
        
        denoise{i}=denoise_img; % save denoise image to cell array       
     end
  
 
    
    H(:,1)=T; %alpha1
    H(:,2)=V; %alpha0
    H(:,3)=U; %MSSIM
    H(:,4)=P; %PSNR
    H(:,5)=Q; %MS SSIM
    H(:,6)=X; %MSE
    H(:,7)=Y; %Brisque
    H(:,8)=Z; %Niqe
    H(:,9)=V1; %lambda
    H(:,10)=R; %Piqe

    result.I(k).factor=H;
    image.I(k).den_img=denoise; %save cell array to image struct

    ans4=sprintf('finish update H, continue to next noise level')
end    
toc
%% Save data:
    savefile=sprintf('DATA/uc_%d_%d_par.mat',M,N) % create file name
    save(savefile,'result'); % save struct result to file only store parameters
    
    saveimage=sprintf('DATA/uc_%d_%d_img.mat',M,N)
    save(saveimage,'image');
    
    ans5=sprintf('finish UC code')
end

    
        