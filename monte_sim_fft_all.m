% This code run and plot Monte Carlo simulation based on FFT code
% There are five random parameters in this simulation which is alpha(alpha
% 1), beta(alpha 0), seed, mean and standard deviation.
% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function monte_sim_fft_all(M,N)
%% Preallocate:
H=zeros(N,9);   % preallocate storage for all value
T=zeros(N,1);   % preallocate storage for alpha
V=zeros(N,1);   % preallocate storage for beta
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for MSE
Y=zeros(N,1);   % preallocate storage for Brisque
Z=zeros(N,1);   % preallocate storage for Niqe
R=zeros(N,1);   % preallocate storage for Piqe

denoise=cell(N,1); % create cell array for denoise images storage

result=struct();  % create struct for all results storage
image=struct(); %create struct for all images storage

%% ADMM parameters
nite = 100; % number of iterations
% balancing weights for Total Variation
% alpha = 0.09;  % 1st order
% beta = 0.11; % 2nd order

% load an image
%% Loading data   

ffcdata=load('GroundTruth');        %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);       %load a Ground Truth image

tic
%% Adding noise

for k=1:M
    
    ans1=sprintf('FFT code create noise image %d/%d',k,M)
    
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

image.I(k).n_img=noise_img; %store created noise image to struct
	

%% Monte Carlo simulation:

%count=0;
format long;

     parfor i=1:N   % parallel computing
         
        ans3=sprintf('step %d/%d',i,N)
        
        % Random seed:
        rng(i,'twister') %for different seed in different stream
        alpha=rand(1,1);
        beta=rand(1,1);
       
        % Store parameters
        T(i)=alpha;
        V(i)=beta;
        
        %Call TGV:
        denoise_img=calltgv(GroTru,noise_img,alpha,beta,nite);        
        %for c = 1:size(GroTru,3)
        %denoise_img(:,:,c) = fft_tgv( noise_img(:,:,c), alpha, beta, nite );
        %end
        
        %SSIM
        [mssim, ssim_map] = ssim(denoise_img, GroTru);
        
        %Multi Scale SSIM
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

        % Dsiplaying results
        ans4=sprintf('alpha: %.3f, beta: %.3f and mssim: %.8f', alpha, beta, mssim)
   
        
        denoise{i}=denoise_img; % save denoise image to cell array       
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
    image.I(k).den_img=denoise; %save cell array to image struct

    ans5=sprintf('finish update H, continue to next noise level')
end

    toc
          
    savefile=sprintf('DATA/fft_%d_%d_all.mat',M,N) %create file name
    save(savefile,'result'); % save struct result to file
    
    saveimage=sprintf('DATA/fft_%d_%d_img.mat',M,N)
    save(saveimage,'image');
    
    ans6=sprintf('finish FFT code')
    
end



%% call TGV function for parallel computating
function output=calltgv(GroTru,noise_img,alpha,beta,nite)
output=zeros(size(noise_img)); % preallocate
for c = 1:size(GroTru,3)
        output(:,:,c) = fft_tgv( noise_img(:,:,c), alpha, beta, nite );
end
        
end