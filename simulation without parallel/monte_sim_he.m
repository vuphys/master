
function monte_sim_he(N)

%blur kernel definition
K     =   fspecial('average',1); %for denoising

%K     =   fspecial('Gaussian',9,3); % for debluring
OrgSigma=0.1;
%fprintf('The noise std of the observed image: %g.\n', OrgSigma); 


%% Preallocate:
H=zeros(N,9);   % preallocate storage for all value
T=zeros(N,1);   % preallocate storage for alpha_1
V=zeros(N,1);   % preallocate storage for alpha_0
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for MSE
Y=zeros(N,1);   % preallocate storage for Brisque
Z=zeros(N,1);   % preallocate storage for ESSIM
V1=zeros(N,1);  % preallocate storage for beta

%% Loading data   

ffcdata=load('GroundTruth');    %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);  %load a Ground Truth image
M     =   numel(GroTru);             
[m,n] =   size(GroTru);

%% adding noise

blur_im  = imfilter(GroTru,K,'circular','conv');

BSNR=20*log10(norm(blur_im(:)-mean(blur_im(:)),'fro')/sqrt(M)/OrgSigma);
%fprintf('BSNR of the observed image: %g dB.\n', BSNR);

%% Adding noise
rng(0);
noise_img = GroTru+randn(m,n)*0.1; % white Gaussian noise added to the image
	
%f        = blur_im + 0.1*randn(m,n);         %add noise
%PSNR_F    =psnr(f,y);
%fprintf('PSNR of the observed image: %g dB.\n\n', PSNR_F);


tao0=0.006; tao1=0.03;% slightly tuning may cause more appealing result 
if size(K)==1
tao  = -BSNR*tao1+1.09; %for denoising
else
tao  = -BSNR*tao0+1.09; %for deblurring
end
c    =  tao*m*n*OrgSigma.^2; % upper bound for the constraint

%% Monte Carlo simulation:

count=0;
format long;

 % Parameter of this TGV function:
        Param.OrigIm     = GroTru;      
        Param.MaxIter    = 1000; 
        Param.SolRE      = 1e-6;    
        Param.UpBound    = c;
        Param.Beta       = 1;       
        Param.Gamma      = 1;
        Param.Tao        = 1;       
        Param.BSNR       = BSNR;     
tic
        for i=1:N       
        % Random seed:
        rng shuffle
        Param.alpha1=rand(1,1)*5;
        rng shuffle
        Param.alpha0=rand(1,1)*5;
        while(Param.Beta==0)
        rng shuffle    
        Param.Beta=rand(1,1)*5;
        end
        % Store parameters
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
        [mulmssim,mulssim_map]=multissim(denoise_img,GroTru);
        
        %MSE
        err = immse(denoise_img, GroTru);
        
        %Brisque
        score = brisque(denoise_img);

        %PSNR
        psnr_tgv = psnr(denoise_img,GroTru);
        
        %ESSIM
        essim_value=ESSIM(GroTru,denoise_img);
        
        % Store MSSIM and PSNR:
        U(i)=mssim;
        P(i)=psnr_tgv;
        Q(i)=mulmssim;
        X(i)=err;
        Y(i)=score;
        Z(i)=essim_value;
        count=count+1;

        % Dsiplaying results
        sprintf('step: %d, alpha1: %.3f, alpha0: %.3f, beta: %.3f and mssim: %.8f', count, Param.alpha1, Param.alpha0, Param.Beta, mssim)
        
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
    H(:,8)=Z; %ESSIM
    H(:,9)=V1; %beta
         
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

%% find maximum ESSIM:

[zmax3,loc5] = max(H(:,8));
a5=loc5;
xmax5=H(a5,1);
ymax5=H(a5,2);
beta5=H(a5,9);
psnr_value5=H(a5,4);
output5=sprintf('Step: %d, max ESSIM= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta: %.6f',a5,zmax3, psnr_value5,xmax5,ymax5,beta5)

    
 %% Plotting:
 
 % Plot dot alpha0-alpha1-SSIM
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