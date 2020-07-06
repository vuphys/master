%% Ask user to input image from computer:
sprintf('Please input your noised image')
[noise_img, InPath] = uigetfile('*.png', 'Please choose an image for denoising:');
if ~ischar(noise_img)
  disp('You did not choose an image. Program is closing');
  return;
end

noise_img  = fullfile(InPath, noise_img);
noise_img= im2double( imread( noise_img ) );
noise_img=noise_img(:,:,1);
%OutFile = fullfile(OutPath, OutFile);
figure(1)
imshow(noise_img)
title( sprintf('Your input noised image' ) );

%% Ask user to choose algorithm and metric for assessment
prompt1='Please choose algorithm (press 1,2 or 3)\n 1.Condat\n  2.FFT\n  3.UC\n Your choice: ';
f1=input(prompt1,'s');%save as string
prompt2='Please choose metric (press a value from 1 to 6)\n 1.SSIM\n 2.PSNR\n 3.MS-SSIM\n 4.PSNR HVSM\n 5.PSNR HMA\n 6.FSIM\n Your choice: ';
f2=input(prompt2,'s');

%% Load Ground Truth image
fname = 'test_image/GroTru.png';
GroTru = im2double( imread( fname ) );
GroTru=GroTru(:,:,1);


%% Automatical TGV filtering using fitting result:
if strcmp(f1,'1') %Condat algorithm
    if strcmp(f2,'1') %Condat SSIM
        x=ssim(noise_img,GroTru)
        alpha_1=-7.172*x^3+12.43*x^2-6.719*x+1.223
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=condat_tgv(noise_img,alpha_1,alpha_0,0.01,1000);
            run_time=toc
            figure(5)
            imshow(denoise_img)
        elseif strcmp(f3,'n') 
            tic
        denoise_img=condat_tgv(noise_img,alpha_1,100,0.01,1000);
        run_time=toc
        figure(5)
        imshow(denoise_img)
        end
        y=ssim(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
title( sprintf('From the left,  GT image,  noised image with SSIM= %.4f,  denoised image with SSIM=  %.4f',x,y ) );
    elseif strcmp(f2,'2') %Condat PSNR
         x=psnr(noise_img,GroTru)
        alpha_1=-4.865*10^(-5)*x^3+0.003814*x^2-0.0987*x+0.9125
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=condat_tgv(noise_img,alpha_1,alpha_0,0.01,1000);
            run_time=toc
            figure(5)
            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=condat_tgv(noise_img,alpha_1,100,0.01,1000);
        run_time=toc
        figure(5)
        imshow(denoise_img)
        end
        y=psnr(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR= %.4f,  denoised image with PSNR=  %.4f', x,y ) );

    elseif strcmp(f2,'3') %Condat MS SSIM
        x=multissim(noise_img,GroTru)
        alpha_1=-1.67*x^3+4.606*x^2-5.164*x+2.24
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=condat_tgv(noise_img,alpha_1,alpha_0,0.01,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=condat_tgv(noise_img,alpha_1,100,0.01,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=multissim(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with MS SSIM= %.4f,  denoised image with MS SSIM=  %.4f', x,y ) );

    elseif strcmp(f2,'4') %Condat PSNR HVSM
        x=psnrhvsm(noise_img,GroTru)
        alpha_1=-3.849*10^(-5)*x^3+0.008971*x^2-0.6954*x+18
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=condat_tgv(noise_img,alpha_1,alpha_0,0.01,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=condat_tgv(noise_img,alpha_1,100,0.01,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnrhvsm(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR HVSM= %.4f,  denoised image with PSNR HVSM=  %.4f', x,y ) );

    elseif strcmp(f2,'5') %Condat PSNR HMA
        x=psnrhma(noise_img,GroTru)
        alpha_1=-8.704*10^(-5)*x^3+0.01991*x^2-1.508*x+37.92
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=condat_tgv(noise_img,alpha_1,alpha_0,0.01,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=condat_tgv(noise_img,alpha_1,100,0.01,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnrhma(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR HMA= %.4f,  denoised image with PSNR HMA=  %.4f', x,y ) );
    else %Condat FSIM
        x=FeatureSIM(GroTru,noise_img)
        alpha_1=-208*x^3+608.1*x^2-595.7*x+195.7
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=condat_tgv(noise_img,alpha_1,alpha_0,0.01,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=condat_tgv(noise_img,alpha_1,100,0.01,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=FeatureSIM(GroTru,denoise_img);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with FSIM= %.4f,  denoised image with FSIM=  %.4f', x,y ) );
    end
elseif strcmp(f1,'2') %FFT algorithm
     if strcmp(f2,'1') %FFT SSIM
        x=ssim(noise_img,GroTru)
        alpha_1=-4.667*x^3+8.782*x^2-5.195*x+1.031
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=calltgv(GroTru,noise_img,alpha_1,alpha_0,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=calltgv(GroTru,noise_img,alpha_1,100,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=ssim(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
title( sprintf('From the left,  GT image,  noised image with SSIM= %.4f,  denoised image with SSIM=  %.4f', x,y ) );
    elseif strcmp(f2,'2') %FFT PSNR
         x=psnr(noise_img,GroTru)
        alpha_1=-2.213*10^(-5)*x^3+0.002388*x^2-0.07788*x+0.8303
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=calltgv(GroTru,noise_img,alpha_1,alpha_0,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=calltgv(GroTru,noise_img,alpha_1,100,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnr(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR= %.4f,  denoised image with PSNR=  %.4f', x,y ) );

    elseif strcmp(f2,'3') %FFT MS SSIM
        x=multissim(noise_img,GroTru)
        alpha_1=-0.9663*x^3+3.163*x^2-4.187*x+2.003
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=calltgv(GroTru,noise_img,alpha_1,alpha_0,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=calltgv(GroTru,noise_img,alpha_1,100,10000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=multissim(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with MS SSIM= %.4f,  denoised image with MS SSIM=  %.4f', x,y ) );

    elseif strcmp(f2,'4') %FFT PSNR HVSM
        x=psnrhvsm(noise_img,GroTru)
        alpha_1=-1.727*10^(-5)*x^3+0.004565*x^2-0.3952*x+11.28
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=calltgv(GroTru,noise_img,alpha_1,alpha_0,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=calltgv(GroTru,noise_img,alpha_1,100,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnrhvsm(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR HVSM= %.4f,  denoised image with PSNR HVSM=  %.4f', x,y ) );

    elseif strcmp(f2,'5') %PSnR HMA
        x=psnrhma(noise_img,GroTru)
        alpha_1=-2.235*10^(-5)*x^3+0.005898*x^2-0.5095*x+14.49
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=calltgv(GroTru,noise_img,alpha_1,alpha_0,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=calltgv(GroTru,noise_img,alpha_1,100,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        y=psnrhma(denoise_img,GroTru);
        end
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR HMA= %.4f,  denoised image with PSNR HMA=  %.4f', x,y ) );
     else %FFT FSIM
        x=FeatureSIM(GroTru,noise_img)
        alpha_1=-146.9*x^3+434.3*x^2-430.4*x+143.1
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=calltgv(GroTru,noise_img,alpha_1,alpha_0,1000);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=calltgv(GroTru,noise_img,alpha_1,100,1000);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=FeatureSIM(GroTru,denoise_img);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with FSIM= %.4f,  denoised image with FSIM=  %.4f', x,y ) );
     end
else % UC algorithm
       if strcmp(f2,'1') %UC SSIM
        x=ssim(noise_img,GroTru)
        alpha_1=-7.435*x^3+12.83*x^2-6.851*x+1.25
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=uc_tgv(noise_img, 1, alpha_0, alpha_1,10000, 100);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=uc_tgv(noise_img, 1, 100, alpha_1,10000, 100);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=ssim(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
title( sprintf('From the left,  GT image,  noised image with SSIM= %.4f,  denoised image with SSIM=  %.4f', x,y ) );
    elseif strcmp(f2,'2') %UC PSNR
         x=psnr(noise_img,GroTru)
        alpha_1=-4.875*10^(-5)*x^3+0.003874*x^2-0.09977*x+0.9169
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=uc_tgv(noise_img, 1, alpha_0, alpha_1,10000, 100);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=uc_tgv(noise_img, 1, 100, alpha_1,10000, 100);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnr(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR= %.4f,  denoised image with PSNR=  %.4f', x,y ) );

    elseif strcmp(f2,'3') %UC MS SSIM
        x=multissim(noise_img,GroTru)
        alpha_1=-0.2461*x^3+1.635*x^2-3.21*x+1.842
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=uc_tgv(noise_img, 1,alpha_0, alpha_1,10000, 100);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=uc_tgv(noise_img, 1, 100, alpha_1,10000, 100);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=multissim(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with MS SSIM= %.4f,  denoised image with MS SSIM=  %.4f', x,y ) );

    elseif strcmp(f2,'4') %UC PSNR HVSM
        x=psnrhvsm(noise_img,GroTru)
        alpha_1=-3.843*10^(-5)*x^3+0.009003*x^2-0.6998*x+18.14
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=uc_tgv(noise_img, 1, alpha_0, alpha_1,10000, 100)
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=uc_tgv(noise_img, 1, 100, alpha_1,10000, 100);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnrhvsm(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR HVSM= %.4f,  denoised image with PSNR HVSM=  %.4f', x,y ) );

    elseif strcmp(f2,'5') %UC PSNR HMA
        x=psnrhma(noise_img,GroTru)
        alpha_1=-0.0001285*x^3+0.02946*x^2-2.228*x+55.71
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=uc_tgv(noise_img, 1, alpha_0, alpha_1,10000, 100);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=uc_tgv(noise_img, 1, 100, alpha_1,10000, 100);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=psnrhma(denoise_img,GroTru);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with PSNR HMA= %.4f,  denoised image with PSNR HMA=  %.4f', x,y ) );
       else % UC PSNR FSIM
        x=FeatureSIM(GroTru,noise_img)
        alpha_1=-219.7*x^3+640.2*x^2-625.4*x+204.9
        prompt3='Default value of alpha_0 is 100 \n Do you want to change it (y/n)\n Your choice:';
        f3=input(prompt3,'s');
        if strcmp(f3,'y')
            prompt4='Please choose your value:'
            alpha_0=str2double(input(prompt4,'s'))
            tic
            denoise_img=uc_tgv(noise_img, 1, alpha_0, alpha_1,10000, 100);
            run_time=toc
                        figure(5)

            imshow(denoise_img)
        elseif strcmp(f3,'n')
            tic
        denoise_img=uc_tgv(noise_img, 1, 100, alpha_1,10000, 100);
        run_time=toc
                    figure(5)

        imshow(denoise_img)
        end
        y=FeatureSIM(GroTru,denoise_img);
        figure(100)
        imshow( [GroTru, noise_img, denoise_img] )
        title( sprintf('From the left,  GT image,  noised image with FSIM= %.4f,  denoised image with FSIM=  %.4f', x,y ) );
     end
end

%% call TGV function for FFT algorithm
function output=calltgv(GroTru,noise_img,alpha,beta,nite)
output=zeros(size(noise_img)); % preallocate
for c = 1:size(GroTru,3)
        output(:,:,c) = fft_tgv( noise_img(:,:,c), alpha, beta, nite );
end
        
end