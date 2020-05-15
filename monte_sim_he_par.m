% This code run and plot Monte Carlo simulation based on He code
% There are several parameters in He TGV algorithm, this code run
% 3 main parameters which is alpha 0, alpha 1, and beta
% For a full understanding of paramters please find more information in "he_tgv.m"

function monte_sim_he_par(N)

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

%% Loading data   

ffcdata=load('GroundTruth');    %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);  %load a Ground Truth image
M     =   numel(GroTru);  %determine the number of elements in Ground Truth image matrix 
[m,n] =   size(GroTru);

%% Adding noise

blur_im  = imfilter(GroTru,K,'circular','conv');

BSNR=20*log10(norm(blur_im(:)-mean(blur_im(:)),'fro')/sqrt(M)/OrgSigma);
%fprintf('BSNR of the observed image: %g dB.\n', BSNR);

rng(0); %reproducibility
noise_img = GroTru+randn(m,n)*0.1; % white Gaussian noise added to the image
	
% Calculate parameter tao for this TGV algorithm
tao0=0.006; tao1=0.03;% slightly tuning may cause more appealing result 
if size(K)==1
tao  = -BSNR*tao1+1.09; %for denoising
else
tao  = -BSNR*tao0+1.09; %for deblurring
end
c    =  tao*m*n*OrgSigma.^2; % upper bound for the constraint

%% Monte Carlo simulation:

%count=0;
format long;

 tic
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
        sprintf('alpha1: %.3f, alpha0: %.3f, beta: %.3f and mssim: %.8f',Param.alpha1, Param.alpha0, Param.Beta, mssim)
        
        end
 toc    

%% Save data:
    
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
         
    savefile=sprintf('DATA/he_%d_v1.mat',N); %create file name 
    save(savefile,'H');              % save matrix H to file
    
 %% find maximum MSSIM:

[zmax1,loc1] = max(H(:,3));
a1=loc1;
xmax1=H(a1,1);
ymax1=H(a1,2);
beta1=H(a1,9);
psnr_value1=H(a1,4);
output1=sprintf('Step: %d, max MSSIM= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f',a1,zmax1, psnr_value1,xmax1,ymax1,beta1)

%% find maximum Multi MSSIM:

[zmax2,loc2] = max(H(:,5));
a2=loc2;
xmax2=H(a2,1);
ymax2=H(a2,2);
beta2=H(a2,9);
psnr_value2=H(a2,4);
output2=sprintf('Step: %d, max Multi MSSIM= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f',a2,zmax2, psnr_value2,xmax2,ymax2,beta2)

%% find minimum MSE:

[zmin,loc3] = min(H(:,6));
a3=loc3;
xmin=H(a3,1);
ymin=H(a3,2);
beta3=H(a3,9);
psnr_value3=H(a3,4);
output3=sprintf('Step: %d, min MSE= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f',a3,zmin, psnr_value3, xmin,ymin,beta3)

%% find minimum Brisque:

[zmin1,loc4] = min(H(:,7));
a4=loc4;
xmin1=H(a4,1);
ymin1=H(a4,2);
beta4=H(a4,9);
psnr_value4=H(a4,4);
mssim4=H(a4,3);
output4=sprintf('Step: %d, min brisque= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f, MSSIM: %.6f',a4,zmin1, psnr_value4, xmin1,ymin1,beta4,mssim4)

%% find minimum Niqe:

[zmin2,loc5] = min(H(:,8));
a5=loc5;
xmin2=H(a5,1);
ymin2=H(a5,2);
lamb5=H(a5,9);
psnr_value5=H(a5,8);
mssim5=H(a5,3);
output5=sprintf('Step: %d, min niqe= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, lambda: %.6f, MSSIM: %.6f',a5,zmin2, psnr_value5, xmin2,ymin2,lamb5, mssim5)

%% find minimum Piqe:

[zmin3,loc6] = min(H(:,10));
a6=loc6;
xmin3=H(a6,1);
ymin3=H(a6,2);
lamb6=H(a6,9);
psnr_value6=H(a6,8);
mssim6=H(a6,3);
output6=sprintf('Step: %d, min Piqe= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, lambda: %.6f, MSSIM: %.6f',a6,zmin3, psnr_value6, xmin3,ymin3,lamb6, mssim6)
    
 %% Plotting:
 
 % Plot alpha0-alpha1-SSIM
 figure(1)
 plot3(H(:,1),H(:,2),H(:,3),'.')
 title('He code \alpha_0 vs \alpha_1 vs SSIM')
 xlabel('\alpha_1 value')
 ylabel('\alpha_0 value')
 zlabel('Mean SSIM')
 
 %figure(2)
% plot3(H(:,1),H(:,9),H(:,3),'b.')
 
 %figure(3)
 %plot3(H(:,2),H(:,9),H(:,3),'b.')
 toc
end


