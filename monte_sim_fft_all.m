% This code run and plot Monte Carlo simulation based on FFT code
% There are five random parameters in this simulation which is alpha(alpha
% 1), beta(alpha 0), seed, mean and standard deviation.
% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function monte_sim_fft_all(M,N)

tic
%% Preallocate:
H=zeros(N,9);   % preallocate storage for all value
T=zeros(N,1);   % preallocate storage for alpha
V=zeros(N,1);   % preallocate storage for beta
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for PSNR HVS
Y=zeros(N,1);   % preallocate storage for PSNR HMA
Z=zeros(N,1);   % preallocate storage for VIF
R=zeros(N,1);   % preallocate storage for FSIM
A=zeros(1,7);   % preallocate storage for noise image metric score
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


%% Adding noise

for k=1:M
    
    ans1=sprintf('FFT code create noise image %d/%d',k,M)
    
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
        Z(i)=vif_val;
        R(i)=fsim;
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
    H(:,6)=X; %PSNR HVS
    H(:,7)=Y; %PSNR HMA
    H(:,8)=Z; %VIF
    H(:,9)=R;  %FSIM
    
    result.I(k).factor=H;
    image.I(k).den_img=denoise; %save cell array to image struct

    ans5=sprintf('finish update H, continue to next noise level')
end

   
    
    %% Create a data .txt file for later statistical process
allPar=zeros(M*N,13); % matrix store all random parameter and results
                     % of all image
h=0;
for i=1:M           % paste all random parameters of noise and TGV to allPar
    for j=1:N
        allPar(1+h:N+h,1)=i;
        allPar(1+h:N+h,2)=result.I(i).seed;
        allPar(1+h:N+h,3)=result.I(i).mean;
        allPar(1+h:N+h,4)=result.I(i).std;
        allPar(1+h:N+h,5)=result.I(i).factor(:,1); %alpha(alpha1)
        allPar(1+h:N+h,6)=result.I(i).factor(:,2); %beta(alpha0)
        allPar(1+h:N+h,7)=result.I(i).factor(:,3); %SSIM
        allPar(1+h:N+h,8)=result.I(i).factor(:,4); %PSNR
        allPar(1+h:N+h,9)=result.I(i).factor(:,5); %MS SSIM
        allPar(1+h:N+h,10)=result.I(i).factor(:,6); %PSNR HVS
        allPar(1+h:N+h,11)=result.I(i).factor(:,7); %PSNR HMA
        allPar(1+h:N+h,12)=result.I(i).factor(:,8); %VIF
        allPar(1+h:N+h,13)=result.I(i).factor(:,9); %FSIM
        h=i*N;
        if h==M*N
            break;
        end
    end
end

    % Give index name for data:
    tabPar=array2table(allPar,'VariableNames',{'img','seed'...
        ,'mean','std','alpha_1','alpha_0','SSIM','PSNR',...
        'MS_SSIM','PSNR_HVS','PSNR_HMA','VIF','FSIM'});
      
    % Save data as .txt   
    txtfile=sprintf('DATA/datFft_%d_%d.txt',M,N)
    writetable(tabPar,txtfile,'Delimiter','tab');

%% Save data as matlab type          
    savefile=sprintf('DATA/fft_%d_%d_par.mat',M,N) %create file name
    save(savefile,'result'); % save struct result to file only store parameters
    
    saveimage=sprintf('DATA/fft_%d_%d_img.mat',M,N)
    save(saveimage,'image','-v7.3');
    
    ans6=sprintf('finish FFT code')
toc    
end



%% call TGV function for parallel computing
function output=calltgv(GroTru,noise_img,alpha,beta,nite)
output=zeros(size(noise_img)); % preallocate
for c = 1:size(GroTru,3)
        output(:,:,c) = fft_tgv( noise_img(:,:,c), alpha, beta, nite );
end
        
end