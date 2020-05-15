% This code run and plot Monte Carlo simulation based on He code
% There are several parameters in He TGV algorithm, this code run
% 6 main parameters which are alpha 0, alpha 1, beta for TGV and
% seed, mean, standard deviation for random Gaussian noise image
% For a full understanding of paramters please find more information in "he_tgv.m"

% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function monte_sim_he_all(M,N)

%blur kernel definition
K     =   fspecial('average',1); %for denoising

%K     =   fspecial('Gaussian',9,3); % for debluring
OrgSigma=0.1;
%fprintf('The noise std of the observed image: %g.\n', OrgSigma); 


%% Preallocate:
H=zeros(N,10);   % preallocate storage for all value
T=zeros(N,1);   % preallocate storage for alpha_1
V=zeros(N,1);   % preallocate storage for alpha_0
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for MSE
Y=zeros(N,1);   % preallocate storage for Brisque
Z=zeros(N,1);   % preallocate storage for niqe
V1=zeros(N,1);  % preallocate storage for beta
R=zeros(N,1);   %preallocate storage for piqe

denoise=cell(N,1); % create cell array for denoise images storage

result=struct();   % create struct for all results storage

%% Loading data   

ffcdata=load('GroundTruth');    %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);  %load a Ground Truth image
B     =   numel(GroTru);  %determine the number of elements in Ground Truth image matrix 
[m,n] =   size(GroTru);

%% Adding noise

blur_im  = imfilter(GroTru,K,'circular','conv');

BSNR=20*log10(norm(blur_im(:)-mean(blur_im(:)),'fro')/sqrt(B)/OrgSigma);
%fprintf('BSNR of the observed image: %g dB.\n', BSNR);

%rng(0); %reproducibility
%noise_img = GroTru+randn(m,n)*0.1; % white Gaussian noise added to the image
	
% Calculate parameter tao for this TGV algorithm
tao0=0.006; tao1=0.03;% slightly tuning may cause more appealing result 
if size(K)==1
tao  = -BSNR*tao1+1.09; %for denoising
else
tao  = -BSNR*tao0+1.09; %for deblurring
end
c    =  tao*m*n*OrgSigma.^2; % upper bound for the constraint

tic

for k=1:M
    
    ans1=sprintf('He code create noise image %d/%d',k,M)
    
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

 
        parfor i=1:N  % parallel computing
            

        % Input parameters:    
        Param = struct();
        Param.OrigIm     = GroTru;      
        Param.MaxIter    = 1000; 
        Param.SolRE      = 1e-6;    
        Param.UpBound    = c;
        Param.Beta       = 0;       
        Param.Gamma      = 1;
        Param.Tao        = 1;       
        Param.BSNR       = BSNR; 
        
        % Random seed:
        rng(i,'twister') %for different seed in different stream
        Param.alpha1=rand(1,1)*5;
        Param.alpha0=rand(1,1)*5;
        
        while(Param.Beta==0) % Make sure Beta is not zero
        Param.Beta=rand(1,1)*5;
        end
        
        % Store random parameters
        T(i)=Param.alpha1;
        V(i)=Param.alpha0;
        V1(i)=Param.Beta;
        
        % Call TGV:
        output = he_tgv(noise_img, K, Param); %% main program
        
        denoise_img     = output.Sol; 
        
        %Reglambda  = output.Reglambda;
        %PSNR       = output.PSNR;    
        %mse        = output.MSE;
        %IterTime   = output.IterTime; 
        %Fvalue     = output.Fvalue;
        
        %SSIM
        [mssim, ssim_map] = ssim(denoise_img,GroTru);
        
        %MS SSIM
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
        an4=sprintf('alpha1: %.3f, alpha0: %.3f, beta: %.3f and mssim: %.8f',Param.alpha1, Param.alpha0, Param.Beta, mssim)
        
        denoise{i}=denoise_img; % save denoise image to cell array       
        end
    H(:,1)=T; %alpha_1
    H(:,2)=V; %alpha_0
    H(:,3)=U; %MSSIM
    H(:,4)=P; %PSNR
    H(:,5)=Q; %MS SSIM
    H(:,6)=X; %MSE
    H(:,7)=Y; %Brisque
    H(:,8)=Z; %niqe
    H(:,9)=V1; %beta
    H(:,10)=R; %piqe        
      
    result.I(k).factor=H;
    result.I(k).den_img=denoise; %save cell array to result struct

    ans5=sprintf('finish update H, continue to next noise level')

end

toc
%% Save data:
               
    savefile=sprintf('DATA/he_%d_%d_all.mat',M,N); %create file name 
    save(savefile,'result');              % save struct result to file
    
    ans6=sprintf('finish He code')
end