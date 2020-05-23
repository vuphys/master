% This code run and plot Monte Carlo simulation based on He code
% There are several parameters in He TGV algorithm, this code run
% 6 main parameters which are alpha 0, alpha 1, beta for TGV and
% seed, mean, standard deviation for random Gaussian noise image
% For a full understanding of paramters please find more information in "he_tgv.m"

% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function monte_sim_he_all(M,N)

tic
%blur kernel definition
K     =   fspecial('average',1); %for denoising

%K     =   fspecial('Gaussian',9,3); % for debluring
OrgSigma=0.1;
%fprintf('The noise std of the observed image: %g.\n', OrgSigma); 


%% Preallocate:
H=zeros(N,10);  % preallocate storage for denoise image metric scores and parameters
T=zeros(N,1);   % preallocate storage for alpha_1
V=zeros(N,1);   % preallocate storage for alpha_0
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for PSNR HVS
Y=zeros(N,1);   % preallocate storage for PSNR HMA
Z=zeros(N,1);   % preallocate storage for VIF
V1=zeros(N,1);  % preallocate storage for beta
R=zeros(N,1);   %preallocate storage for FSIM
A=zeros(1,7);   % preallocate storage for noise image metric score


denoise=cell(N,1); % create cell array for denoise images storage

result=struct();   % create struct for all results storage
image=struct(); %create struct for all images storage

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



for k=1:M
    
    ans1=sprintf('He code create noise image %d/%d',k,M)
    
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

    % Store metric scores of noise image:

    A(1,1) = ssim(noise_img, GroTru); %SSIM
    A(1,2) = psnr(noise_img,GroTru); %PSNR
    A(1,3) = multissim(noise_img,GroTru); %Multi SSIM
    A(1,4) = psnrhvsm(noise_img, GroTru); %PSNR HVSM
    A(1,5) = psnrhma(noise_img,GroTru); %PSNR HMA
    A(1,6) = vif(GroTru,noise_img); %VIF
    A(1,7) = FeatureSIM(GroTru,noise_img); %FSIM
    
    result.I(k).noise_met=A;
	
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
        
        %PSNR HVS
        psnr_hvs = psnrhvsm(denoise_img, GroTru);
        
        %PSNR HMA
        psnr_hma = psnrhma(denoise_img,GroTru);

        %PSNR
        psnr_tgv = psnr(denoise_img,GroTru);
        
        %VIF
        vif_val=vif(GroTru,denoise_img);
        
        %FSIM
        fsim=FeatureSIM(GroTru,denoise_img);
        
        % Store MSSIM and PSNR:
        U(i)=mssim;
        P(i)=psnr_tgv;
        Q(i)=mulmssim;
        X(i)=psnr_hvs;
        Y(i)=psnr_hma;
        Z(i)=vif_val
        R(i)=fsim;
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
    H(:,6)=X; %PSNR HVS
    H(:,7)=Y; %PSNR HMA
    H(:,8)=Z; %VIF
    H(:,9)=V1; %beta
    H(:,10)=R;  %FSIM
      
    result.I(k).factor=H;
    image.I(k).den_img=denoise; %save cell array to image struct

    ans5=sprintf('finish update H, continue to next noise level')

end

 %% Create a data .txt file for later statistical process
allPar=zeros(M*N,14); % matrix store all random parameter and results
                     % of all image
h=0;
for i=1:M           % paste all random parameters of noise and TGV to allPar
    for j=1:N
        allPar(1+h:N+h,1)=i;
        allPar(1+h:N+h,2)=result.I(i).seed;
        allPar(1+h:N+h,3)=result.I(i).mean;
        allPar(1+h:N+h,4)=result.I(i).std;
        allPar(1+h:N+h,5)=result.I(i).factor(:,1); %alpha1
        allPar(1+h:N+h,6)=result.I(i).factor(:,2); %alpha0
        allPar(1+h:N+h,7)=result.I(i).factor(:,9); %beta
        allPar(1+h:N+h,8)=result.I(i).factor(:,3); %SSIM
        allPar(1+h:N+h,9)=result.I(i).factor(:,4); %PSNR
        allPar(1+h:N+h,10)=result.I(i).factor(:,5); %MS SSIM
        allPar(1+h:N+h,11)=result.I(i).factor(:,6); %PSNR HVS
        allPar(1+h:N+h,12)=result.I(i).factor(:,7); %PSNR HMA
        allPar(1+h:N+h,13)=result.I(i).factor(:,8); %VIF
        allPar(1+h:N+h,14)=result.I(i).factor(:,10); %FSIM
        h=i*N;
        if h==M*N
            break;
        end
    end
end
    % Give index name for data:
    tabPar=array2table(allPar,'VariableNames',{'img','seed'...
        ,'mean','std','alpha_1','alpha_0','beta','SSIM','PSNR',...
        'MS_SSIM','PSNR_HVS','PSNR_HMA','VIF','FSIM'});
   
    % Save data as .txt   
    txtfile=sprintf('DATA/datHe_%d_%d.txt',M,N)
    writetable(tabPar,txtfile,'Delimiter','tab');

%% Save data as matlab type      
               
    savefile=sprintf('DATA/he_%d_%d_par.mat',M,N) %create file name 
    save(savefile,'result');        % save struct result to file only store parameters
    
    saveimage=sprintf('DATA/he_%d_%d_img.mat',M,N)
    save(saveimage,'image','-v7.3');
    
    ans6=sprintf('finish He code')
    
 toc   
end