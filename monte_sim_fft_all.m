% This code run and plot Monte Carlo simulation based on FFT code
% There are five random parameters in this simulation which is alpha(alpha
% 1), beta(alpha 0), seed, mean and standard deviation.
% Input parameters:
% M is the number of generated noise image
% N is the number of MC step

function monte_sim_fft_all(M,N)

tic
%% Preallocate:
H=zeros(N,9);   % preallocate storage for denoise image metric scores and parameters
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
%GroTru = GroTru./(max(abs(GroTru)));

%% Adding noise

for k=1:M
    
    ans1=sprintf('FFT code create noise image %d/%d',k,M)
    
    rng(k,'twister')
    s=round(unifrnd(0,1),8);
    t=s+k+1000; %for reproducibility of the noise image
    rng(t,'twister')
    std=round(unifrnd(0,1),8); % create a uniform distribution random standard deviation for Gaussian noise
                     % from [0,1]
    ans2=sprintf('seed: %.3f, standard deviation: %.3f',t,std)
    
    
    result.I(k).seed=t;     % store seed to struct
    result.I(k).std=std;    % store standard deviation to struct
    
    rng(t)
    noise_img = GroTru+randn(size(GroTru))*std;  % adding Gaussian noise

    image.I(k).n_img=noise_img; %store created noise image to struct

    % Store metric scores of noise image:

    A(1,1) = ssim(noise_img, GroTru); %SSIM
    A(1,2) = psnr(noise_img,GroTru); %PSNR
    A(1,3) = multissim(noise_img,GroTru); %Multi SSIM
    A(1,4) = psnrhvsm(noise_img, GroTru); %PSNR HVSM
    A(1,5) = psnrhma(noise_img,GroTru); %PSNR HMA
    A(1,6) = VIF_FR(GroTru,noise_img); %VIF
    A(1,7) = FeatureSIM(GroTru,noise_img); %FSIM
    
    result.I(k).noise_met=A;
    

%% Monte Carlo simulation:

%count=0;
format long;

 
     parfor i=1:N   % parallel computing
       
        ans3=sprintf('step %d/%d',i,N)
        rng(i)
        checkrand=randi([100,10000],1);                        
        %Random seed:
        rng(k+i*checkrand+100000,'twister') %for different seed in different stream
        alpha=unifrnd(0,1)*2;
        beta=unifrnd(0,1)*2;
                  
        %Call TGV:
        denoise_img=calltgv(GroTru,noise_img,alpha,beta,nite);        
               
        %SSIM
        [mssim, ~] = ssim(denoise_img, GroTru);
                                    
        %Multi Scale SSIM
        [mulmssim,~]=multissim(denoise_img,GroTru);
        
        %PSNR HVSM
        psnr_hvsm = psnrhvsm(denoise_img, GroTru);
        
        %PSNR HMA
        psnr_hma = psnrhma(denoise_img,GroTru);

        %PSNR
        psnr_tgv = psnr(denoise_img,GroTru);
        
        %VIF
        vif_val=VIF_FR(GroTru,denoise_img);
        
        %FSIM
        fsim=FeatureSIM(GroTru,denoise_img);
        
                       
         % Store parameters
        T(i)=alpha;
        V(i)=beta;
        U(i)=mssim;
        P(i)=psnr_tgv;
        Q(i)=mulmssim;
        X(i)=psnr_hvsm;
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
allPar=zeros(M*N,12); % matrix store all random parameter and results
                     % of all image
h=0;
for i=1:M           % paste all random parameters of noise and TGV to allPar
    for j=1:N
        allPar(1+h:N+h,1)=i;
        allPar(1+h:N+h,2)=result.I(i).seed;
        allPar(1+h:N+h,3)=result.I(i).std;
        allPar(1+h:N+h,4)=result.I(i).factor(:,1); %alpha(alpha1)
        allPar(1+h:N+h,5)=result.I(i).factor(:,2); %beta(alpha0)
        allPar(1+h:N+h,6)=result.I(i).factor(:,3); %SSIM
        allPar(1+h:N+h,7)=result.I(i).factor(:,4); %PSNR
        allPar(1+h:N+h,8)=result.I(i).factor(:,5); %MS SSIM
        allPar(1+h:N+h,9)=result.I(i).factor(:,6); %PSNR HVS
        allPar(1+h:N+h,10)=result.I(i).factor(:,7); %PSNR HMA
        allPar(1+h:N+h,11)=result.I(i).factor(:,8); %VIF
        allPar(1+h:N+h,12)=result.I(i).factor(:,9); %FSIM
        h=i*N;
        if h==M*N
            break;
        end
    end
end

    % Give index name for data:
    tabPar=array2table(allPar,'VariableNames',{'img','seed'...
        ,'std','alpha_1','alpha_0','SSIM','PSNR',...
        'MS_SSIM','PSNR_HVS','PSNR_HMA','VIF','FSIM'});
    
    % Create an identify name for data file name based on created time
    datestamp=datetime('now','Format','yyyy-MM-dd''T''HHmmss');
    
    % Save data as .txt   
    txtfile=sprintf('DATA/datFft_%d_%d_%s.txt',M,N,datestamp)
    writetable(tabPar,txtfile,'Delimiter','tab');

%% Save data as matlab type          
    savefile=sprintf('DATA/fft_%d_%d_par_%s.mat',M,N,datestamp) %create parameter file name
    save(savefile,'result'); % save struct result to file only store parameters
    
    saveimage=sprintf('DATA/fft_%d_%d_img_%s.mat',M,N,datestamp) %create image file name
    save(saveimage,'image','-v7.3'); %save all image in a seperate data
    
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