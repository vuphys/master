function OutPut = he_tgv(f, K,Param)
%%
%major references:
%[1]C. He,etal. "A fast adaptive parameter estimation for total...
%variation image restoration," IEEE Trans. Image Process., vol. 23, no.
%12,pp. 4954-4967, Dec. 2014.

%[2]C. He, C. Hu,W. Zhang, B. Shi, and X. Hu, "An adaptive total generalized variation model with augmented
%lagrangian method for image denoising," Mathematical Problems in Engineering,
%vol. 2014, Article ID157893, 11 pages, 2014.

%[3]C. He, C. Hu,W. Zhang, B. Shi, and X. Hu, "Fast...
%total-variation image deconvolution with adaptive parameter estimation via
%split Bregman method," Mathematical Problems in Engineering,
%vol. 2014, Article ID617026, 9 pages, 2014.

%noise variance estimation: Y. Wen and R. H. Chan, "Parameter selection for Total-Variation-Based...
%image restoration using discrepancy principle," IEEE Trans. Image Process., vol. 21, no. 4, pp. 1770¨C1781, Apr. 2012.

%functions C = getC and [D,Dt] = defDDt are referring to FTVD-v4 at: http://www.caam.rice.edu/~optimization/L1/ftvd/ 
%%
% ADMM method is applied to constrained TGV-regularized image deconvolution problem.
% Suppose the mathematical model is given by
%
% $$f = K * u + n$$
% where f, u and n are the observed image with size $n\times m$, the original image 
% and the noise respectively. K is the blur matrix. 
%
% To restore f, We solve the constrained minimization problem
% 
% $$\min_u TGV(u)$$
% subject to $||K*u-f||_2 \leq c$ .
%% Input:
%         Param.OrigIm:     original image                   
%         Param.MaxIter:    default=500;                
%         Param.Sol:        initial solution, default=f 
%         Param.Disp:       display the result or not. default=1;
%         Param.UpBound:    upper bound for the norm of the residual $c$ 
%         Param.Reglambda:  regularization parameter  $\lambda$
%         Param.SolRE:      stop criterion for the relative difference.   default=1e-4.
%         Param.Beta        ADMM parameter.   default=1.
%         Param.Gamma       ADMM parameter.   default=1.618
%%  OutPut:
%         OutPut.Sol:       restored image u;
%         OutPut.PSNR:      the evolution of PSNR against iteration number
%         OutPut.IterTime:  the evolution of CPU time against iteration number

%%=========================================================================
%  Copyright(c), Sep. 2013, Dr.HE Chuan(hechuan8512@163.com)
%%========================================================================= 
u       = f;               MaxIter = 500;        
SolRE   = 1e-4;            beta    = 1;          gamma = 1.618;
gcv     = 0;               flag    = 1;          tao   = 1;
[m,n]   = size(f);
eta1  = zeros(m,n);      eta2  = eta1;        eta3  = eta1;
mu      = eta1;          xi1   = eta1;         xi2   = eta2;
P1      = eta1;          P2      = eta2; 
if nargin == 3
    if isfield(Param,'OrigIm'),     xtrue     = Param.OrigIm;     end
    if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
    if isfield(Param,'Beta'),       beta      = Param.Beta;       end                
    if isfield(Param,'Gamma'),      gamma     = Param.Gamma;      end
    if isfield(Param,'Tao'),        tao       = Param.Tao;        end
    if isfield(Param,'BSNR'),       BSNR      = Param.BSNR;       end
    if isfield(Param,'Disp'),       flag      = Param.Disp;       end
    if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
    if isfield(Param,'Reglambda'),  lambda    = Param.Reglambda;  end
    if isfield(Param,'UpBound'),    c         = Param.UpBound;    
                                    RegMethod = 'Discrepancy';    end
    if isfield(Param,'Sol'),        u         = Param.Sol;        end
    if isfield(Param,'alpha0'),    alpha0         = Param.alpha0;        end
    if isfield(Param,'alpha1'),     alpha1         = Param.alpha1;        end
end

beta1     = 10^(BSNR/10-1)*beta;
beta3     = 1* beta;
beta2     = 1* beta;

alfar0    = alpha0;
alfar1    = alpha1;
    

if ~isvar('c') && ~isvar('lambda') && ~gcv 
     sigma     = ImageStdDev(f);%noise variance estimation
     c         = tao*sigma.^2 * m * n;
     RegMethod = 'Discrepancy';
end

if isvar('lambda')
     RegMethod = 'RegFixed';
end


Param.Sol  = u;
OutPut     = Param;

X = zeros(m,n);
C = getC;
[D,Dt] = defDDt;

regpar = zeros(MaxIter,1);   PSNR=regpar;   iter_time = regpar; mse=regpar; Fvalue = zeros(MaxIter,1);
cont   = 1; k = 0; 
%finite diff
[D1U,D2U] = D(u);
[D1P1,D2P1]  =  D(P1);  [D1P2,D2P2]  =  D(P2);
%start_time = cputime;

while cont
    tic;
    k = k + 1;

     % ==================
     %   Z-subprolem Shrinkage Step
     % ==================
     Z1      =  D1P1 + eta1/beta3;
     Z2      =  D2P2 + eta2/beta3;
     Z3      =  (D2P1 + D1P2)/2 + eta3/beta3;
     W       =  Z1.^2 + Z2.^2 + 2*Z3.^2;
     W       =  sqrt(W);
     W(W==0) = 1;
     W       = max(W - alfar0/beta3, 0)./W;
     Z1      = Z1.*W;
     Z2      = Z2.*W;
     Z3      = Z3.*W;
     % ==================
     %   Y-subprolem Shrinkage Step
     % ==================
     Y1      =  D1U + xi1/beta2 - P1;
     Y2      =  D2U + xi2/beta2 - P2;
     V       =  Y1.^2 + Y2.^2;
     V       =  sqrt(V);
     V(V==0) = 1;
     V       = max(V - alfar1/beta2, 0)./V;
     Y1      = Y1.*V;
     Y2      = Y2.*V;

    % ==================
    %     U-subprolem
    % ==================
    if size(K)==1 %for denoising 
        KtF  = (beta1/beta2).*fft2(X-mu/beta1); 
        Unew = (KtF + fft2(Dt(Y1+P1-xi1/beta2,Y2+P2-xi2/beta2)))./(C.eigsDtD + (beta1/beta2));
    else  %for deblurring  
        KtF  = (beta1/beta2)*(conj(C.eigsK) .*fft2(X-mu/beta1)); 
       Unew = KtF + fft2(Dt(Y1+P1-xi1/beta2,Y2+P2-xi2/beta2));
        Unew = Unew./(C.eigsDtD + (beta1/beta2)*C.eigsKtK);
    end
    unew = real(ifft2(Unew));
    % ==================
    %     P-subprolem
    % ================== 
    [D1unew,D2unew] = D(unew);
    P1 = (beta3/beta2)*(Dt(Z1 - eta1/beta3, Z3 - eta3/beta3 - D1P2/2) + D1unew + xi1/beta2 - Y1);
    P1 = fft2(P1)./((beta3/beta2)* (C.eigsDtD1 + C.eigsDtD2/2) + ones(m,n));
    P1 = real(ifft2(P1));
    [D1P1,D2P1]  =  D(P1);
    P2 = (beta3/beta2)*(Dt(Z3 - eta3/beta3 - D2P1/2, Z2 - eta2/beta3) + D2unew + xi2/beta2 - Y2);
    P2 = fft2(P2)./((beta3/beta2)* (C.eigsDtD1/2 + C.eigsDtD2) + ones(m,n));
    P2 = real(ifft2(P2));   
    % ==================
    %   X-subprolem 
    % ==================
    switch RegMethod
        case 'Discrepancy'
            if  size(K)==1 %for denoising
                Kunew   =  real(ifft2(Unew));  
            else  % for deblurring
                Kunew   =  real(ifft2(C.eigsK .*Unew));
            end
            
            A       =  Kunew + mu/beta1;
            FA_norm =  norm(f-A,'fro');
            
        if  FA_norm^2<=c
            lambda  =  0;
            X       =  A;
        else
            lambda  =  beta1* FA_norm/sqrt(c)-beta1; 
            X       =  (lambda*f + beta1*A)/(lambda + beta1);
        end
            
            
        case 'RegFixed'
            Kunew   =  real(ifft2(C.eigsK .*Unew));
            A       =  Kunew + mu/beta1;
            X       = (lambda*f + beta1*A)/(lambda + beta1);
                 
    end
            
    %% update mu, eta and xi
    D1U          = D1unew; D2U          =  D2unew;
    [D1P2,D2P2]  =  D(P2);%[D1P1,D2P1]  =  D(P1);
     mu       = mu     - gamma*beta1*(X  - Kunew);
     eta1   = eta1 - gamma*beta3*(Z1 - D1P1);
     eta2   = eta2 - gamma*beta3*(Z2 - D2P2);
     eta3   = eta3 - gamma*beta3*(Z3 - (D2P1+D1P2)/2);
     xi1    = xi1  - gamma*beta2*(Y1 - D1unew +P1);
     xi2    = xi2  - gamma*beta2*(Y2 - D2unew +P2);       
    %% relative error and checking stopping rule 
    re   = norm(unew-u,'fro')/norm(u,'fro');
     if re<1e-4
       beta1=1.5*beta1;
     end
    u    = unew;
    fvalue = fval;
    cont = (k<MaxIter)&&(re>SolRE);   
    if isvar('xtrue')   
        PSNR(k) = psnr(u,xtrue);
        if mod(k,10)==0 && flag
            fprintf('He code %3d-th  psnr: %2.2f,   regpara: %1.2e  re: %1.2e\n', k,PSNR(k), lambda, re);
        end
    end
    tElapsed     = toc;
    if  k==1       
         iter_time(k)  =tElapsed; 
    else
        iter_time(k) = iter_time(k-1)+tElapsed;
    end
    Fvalue(k)    = fvalue;
    regpar(k)    = lambda;
    mse(k)       = norm(unew-xtrue,'fro')^2/m/n;
end
OutPut.Sol       = u;
OutPut.PSNR      = PSNR(1:k);
OutPut.IterTime  = iter_time(1:k);
OutPut.Reglambda = regpar(1:k);
OutPut.MSE       = mse(1:k);
OutPut.Fvalue    = Fvalue(1:k);
            
%% function used above 
function C = getC
        sizeF = size(f);
        C.eigsK = psf2otf(K,sizeF);
        C.eigsDtD1 =  abs(psf2otf([1,-1],sizeF)).^2;
        C.eigsDtD2 =  abs(psf2otf([1;-1],sizeF)).^2;
        C.eigsDtD  =  C.eigsDtD1 + C.eigsDtD2;
        C.eigsKtK = abs(C.eigsK).^2;
end



function Fvalue = fval
        Fvalue1 = sum(sum(sqrt((D1U-P1).^2 + (D2U-P2).^2)));
        Fvalue2 = sum(sum(sqrt((D1P1).^2 + (D2P2).^2 + 2*(.5*(D1P2 + D2P1)).^2)));
        KU_F    = real(ifft2(C.eigsK .* fft2(u))) - f;
        Fvalue  = Fvalue1 + Fvalue2 + lambda/2 * norm(KU_F,'fro')^2;
end


function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        % referring to FTVD code
        
        D = @(U) ForwardD(U);
        Dt = @(X,Y) Dive(X,Y);
end
        
function [Dux,Duy] = ForwardD(U)
         % Forward finite difference operator
         Dux = [diff(U,1,2), U(:,1) - U(:,end)];
         Duy = [diff(U,1,1); U(1,:) - U(end,:)];
end
        
function DtXY = Dive(X,Y) %Dt=-div
        % Transpose of the forward finite difference operator
        DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
        DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
end

function sigma = ImageStdDev(y)
h = [0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];       
h = h(end:-1:1); %% flip is used only to give exactly same result as previous version of the code (which was using filter2)
        
z = conv2(y,h,'same');
z=conv2(z,h','same');
sigma = median(abs(z(:)))/.6745;
end

function tf = isvar(name)
%function tf = isvar(name)
% determine if "name" is a variable in the caller's workspace

if nargin < 1
	help isvar
	error arg
end

tf = true(1);
evalin('caller', [name ';'], 'tf=logical(0);')
end

            
            
            
            
 
    
    

           
end