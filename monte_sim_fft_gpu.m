% This code run and plot Monte Carlo simulation based on FFT code
% There are five random parameters in this simulation which is alpha(alpha
% 1), beta(alpha 0), seed, mean and standard deviation.
% M is the number of generated noise image
% N is the number of MC step

% Not complete yet

function monte_sim_fft_gpu(M,N)
%% Preallocate:
H=zeros(N,9,'gpuArray');   % preallocate storage for all value
T=zeros(N,1,'gpuArray');   % preallocate storage for alpha
V=zeros(N,1,'gpuArray');   % preallocate storage for beta
U=zeros(N,1,'gpuArray');   % preallocate storage for SSIM
P=zeros(N,1,'gpuArray');   % preallocate storage for PSNR
Q=zeros(N,1,'gpuArray');   % preallocate storage for MS SSIM
X=zeros(N,1,'gpuArray');   % preallocate storage for MSE
Y=zeros(N,1,'gpuArray');   % preallocate storage for Brisque
Z=zeros(N,1,'gpuArray');   % preallocate storage for Niqe
R=zeros(N,1,'gpuArray');   % preallocate storage for Piqe

result=struct();

%% ADMM parameters
nite = gpuArray(100); % number of iterations
% balancing weights for Total Variation
% alpha = 0.09;  % 1st order
% beta = 0.11; % 2nd order

% load an image
%% Loading data   

ffcdata=load('GroundTruth');        %load Ground Truth data
GroTru=gpuArray(ffcdata.data(:,:,1,1));       %load a Ground Truth image

%fname = 'parrotgray.png';
%I = im2double( imread( fname ) );
%I0 = I; % original as the reference

%% Adding noise
tic
for k=1:M
    
    sprintf('create image %d/%d',k,M)
    
    rng(k,'twister')
    s=rand(1,1);
    t=s+k; %for reproducibility of the noise image
    rng(t,'twister')
    std=gpuArray(rand(1,1)*2); % create random standard deviation for Gaussian noise
                     % from [0,2]
    m=gpuArray(2.*rand(1,1)-1);   % create random mean for Gaussian noise
                        % from [-1,1]
    sprintf('seed: %.3f, mean: %.3f, standard deviation: %.3f',t,m,std)
    
    
    result.I(k).seed=s;     % store seed to struct
    result.I(k).mean=m;     % store mean to struct
    result.I(k).std=std;    % store standard deviation to struct
    
%rng(0);                   %reproducibility 
%rng shuffle
x='gaussian';
noise_img = imnoise(GroTru,x,m,std);   % adding Gaussian noise
	

%% Monte Carlo simulation:

%count=0;
format long;


     parfor i=1:gpuArray(N)   % parallel computing
        sprintf('step %d/%d',i,N)
        % Random seed:
        rng(gather(i),'twister') %for different seed in different stream
        alpha=gpuArray(rand(1,1));
        beta=gpuArray(rand(1,1));
       
        % Store parameters
        T(i)=alpha;
        V(i)=beta;
        
        %Call TGV:
        %denoise_img=calltgv(GroTru,noise_img,alpha,beta,nite);
        denoise_img=arrayfun(@calltgv,GroTru,noise_img,alpha,beta,nite);
        %for c = 1:size(GroTru,3)
        %denoise_img(:,:,c) = fft_tgv( noise_img(:,:,c), alpha, beta, nite );
        %end
        
        %SSIM
        [mssim, ssim_map] = gpuArray(ssim(denoise_img, GroTru));
        [mulmssim,mulssim_map]=gpuArray(multissim(denoise_img,GroTru));
        
        %MSE
        err = gpuArray(immse(denoise_img, GroTru));
        
        %Brisque
        brisque_score = gpuArray(brisque(denoise_img));

        %PSNR
        psnr_tgv = gpuArray(psnr(denoise_img,GroTru));
        
        %Niqe
        niqe_score=gpuArray(niqe(denoise_img));
        
        %Piqe
        piqe_score=gpuArray(piqe(denoise_img));
        
        % Store MSSIM and PSNR:
        U(i)=mssim;
        P(i)=psnr_tgv;
        Q(i)=mulmssim;
        X(i)=err;
        Y(i)=brisque_score;
        Z(i)=niqe_score;
        R(i)=piqe_score;
        %count=count+1;

        % Dsiplaying results
        sprintf('alpha: %.3f, beta: %.3f and mssim: %.8f', alpha, beta, mssim)
        
     end
     
     
     
     %% Save data:
     
    H(:,1)=T; %alpha
    H(:,2)=V; %beta
    H(:,3)=U; %MSSIM
    H(:,4)=P; %PSNR
    H(:,5)=Q; %MS SSIM
    H(:,6)=X; %MSE
    H(:,7)=Y; %Brisque
    H(:,8)=Z; %niqe
    H(:,9)=R;  %piqe
    
    result.I(k).factor=H;
    
    sprintf('finish update H, continue to next noise level')
end
     toc
 
         
    savefile=sprintf('DATA/fft_%d_%d_all.mat',M,N); %create file name
    save(savefile,'result'); % save matrix H to file
    
    sprintf('finish code')
    
end

%% call TGV function for parallel computating
function output=calltgv(GroTru,noise_img,alpha,beta,nite)
output=zeros(size(noise_img)); % preallocate
for c = 1:size(GroTru,3)
        output(:,:,c) = fft_tgv( noise_img(:,:,c), alpha, beta, nite );
end
        
end