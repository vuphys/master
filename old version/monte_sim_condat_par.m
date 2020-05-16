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
% There are two random parameters in this simulation which are
% lambda_1(alpha 1) and lambda_2(alpha 0)

% This is an old version without random noise


function monte_sim_condat_par(N)
%% Preallocate:
H=zeros(N,9);   % preallocate storage for all value
R=zeros(N,1);   % preallocate storage for lambda_1
V=zeros(N,1);   % preallocate storage for lambda_2
U=zeros(N,1);   % preallocate storage for SSIM
P=zeros(N,1);   % preallocate storage for PSNR
Q=zeros(N,1);   % preallocate storage for MS SSIM
X=zeros(N,1);   % preallocate storage for MSE
Y=zeros(N,1);   % preallocate storage for Brisque
Z=zeros(N,1);   % preallocate storage for Niqe
T=zeros(N,1);   % preallocate storage for Piqe

%% TGV parameters:
Nbiter= 600;	% number of iterations
%lambda1 = 0.1; 	% regularization parameter
%lambda2 = 0.2;	% regularization parameter
tau = 0.01;		% proximal parameter >0; influences the
                %    convergence speed
%% Loading data   

ffcdata=load('GroundTruth');        %load Ground Truth data
GroTru=ffcdata.data(:,:,1,1);       %load a Ground Truth image
    
%% adding noise

rng(0);                     %reproducibility for next time
noise_img = GroTru+randn(size(GroTru))*0.1;   % adding Gaussian noise


%% Monte Carlo simulation:

%count=0;
format long;
tic
     parfor i=1:N % parallel computing
        
        % Random seed:
        rng(i,'twister') %for different seed in different stream
        lambda1=rand(1,1)*5;
        lambda2=rand(1,1)*5;
       
        % Store parameters
        R(i)=lambda1;
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
        T(i)=piqe_score;
        %count=count+1;

        % Displaying results
        sprintf('lambda1: %.3f, lambda2: %.3f and mssim: %.8f', lambda1, lambda2, mssim)
        
    end
 toc   
    %% Save data:
    
    H(:,1)=R; %lambda1
    H(:,2)=V; %lambda2
    H(:,3)=U; %MSSIM
    H(:,4)=P; %PSNR
    H(:,5)=Q; %MS SSIM
    H(:,6)=X; %MSE
    H(:,7)=Y; %Brisque
    H(:,8)=Z; %Niqe
    H(:,9)=T; %Piqe
         
    savefile=sprintf('DATA/condat_%d_v1.mat',N); % create file name
    save(savefile,'H'); % save matrix H to file
    

    
    %% find maximum MSSIM:

[zmax1,loc1] = max(H(:,3));
a1=loc1;
xmax1=H(a1,1);
ymax1=H(a1,2);
psnr_value1=H(a1,4);
output1=sprintf('Step: %d, max MSSIM= %.8f, PSNR= %.8f \n lambda1= %.6f, lambda2= %.6f',a1,zmax1, psnr_value1,xmax1,ymax1)

%% find maximum Multi MSSIM:

[zmax2,loc2] = max(H(:,5));
a2=loc2;
xmax2=H(a2,1);
ymax2=H(a2,2);
psnr_value2=H(a2,4);
output2=sprintf('Step: %d, max Multi MSSIM= %.8f, PSNR= %.8f \n lambda1= %.6f, lambda2= %.6f',a2,zmax2, psnr_value2,xmax2,ymax2)

%% find minimum MSE:

[zmin,loc3] = min(H(:,6));
a3=loc3;
xmin=H(a3,1);
ymin=H(a3,2);
psnr_value3=H(a3,4);
output3=sprintf('Step: %d, min MSE= %.8f, PSNR= %.8f \n lambda1= %.6f, lambda2= %.6f',a3,zmin, psnr_value3, xmin,ymin)

%% find minimum Brisque:

[zmin1,loc4] = min(H(:,7));
a4=loc4;
xmin1=H(a4,1);
ymin1=H(a4,2);
psnr_value4=H(a4,4);
mssim4=H(a4,3);
output4=sprintf('Step: %d, min brisque= %.8f, PSNR= %.8f \n lambda1= %.6f, lambda2= %.6f, MSSIM: %.6f',a4,zmin1, psnr_value4, xmin1,ymin1, mssim4)

%% find minimum Niqe:

[zmin2,loc5] = min(H(:,8));
a5=loc5;
xmin5=H(a5,1);
ymin5=H(a5,2);
psnr_value5=H(a5,4);
output5=sprintf('Step: %d, min Niqe= %.8f, PSNR= %.8f \n lambda1= %.6f, lambda2= %.6f',a5,zmin2, psnr_value5,xmin5,ymin5)

%% find minimum Piqe:

[zmin3,loc6] = min(H(:,9));
a6=loc6;
xmin6=H(a6,1);
ymin6=H(a6,2);
psnr_value6=H(a6,4);
output6=sprintf('Step: %d, min Piqe= %.8f, PSNR= %.8f \n lambda1= %.6f, lambda2= %.6f',a5,zmin3, psnr_value6,xmin6,ymin6)
    
%% Plotting
    
 % Plot lambda1(alpha1)-lambda2(alpha0)-SSIM
 figure(2)
 plot3(H(:,1),H(:,2),H(:,3),'.')
 title('Condat code \alpha_0 vs \alpha_1 vs SSIM')
 xlabel('\alpha_1 value')
 ylabel('\alpha_0 value')
 zlabel('Mean SSIM')

end