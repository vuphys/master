%% use for he and uc code only

function findpeak_3para(H,n) 
fileID = fopen('log.txt','a+');
%% find maximum MSSIM:
newtext=('******************************************************************************************\n');
fprintf(fileID,newtext);
if (n==0)
    fprintf(fileID,'He code\n');
else 
    fprintf(fileID, 'UC code\n');
end

[zmax1,loc1] = max(H(:,3));
a1=loc1;
xmax1=H(a1,1);
ymax1=H(a1,2);
beta1=H(a1,9);
psnr_value1=H(a1,4);
output1=sprintf('Step: %d, max MSSIM= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f\n',a1,zmax1, psnr_value1,xmax1,ymax1,beta1)
fprintf(fileID,output1);
%% find maximum Multi MSSIM:

[zmax2,loc2] = max(H(:,5));
a2=loc2;
xmax2=H(a2,1);
ymax2=H(a2,2);
beta2=H(a2,9);
psnr_value2=H(a2,4);
output2=sprintf('Step: %d, max Multi MSSIM= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f\n',a2,zmax2, psnr_value2,xmax2,ymax2,beta2)
fprintf(fileID,output2);
%% find minimum MSE:

[zmin,loc3] = min(H(:,6));
a3=loc3;
xmin=H(a3,1);
ymin=H(a3,2);
beta3=H(a3,9);
psnr_value3=H(a3,4);
output3=sprintf('Step: %d, min MSE= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f\n',a3,zmin, psnr_value3, xmin,ymin,beta3)
fprintf(fileID,output3);
%% find minimum Brisque:

[zmin1,loc4] = min(H(:,7));
a4=loc4;
xmin1=H(a4,1);
ymin1=H(a4,2);
beta4=H(a4,9);
psnr_value4=H(a4,4);
mssim4=H(a4,3);
output4=sprintf('Step: %d, min brisque= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta= %.6f, MSSIM: %.6f\n',a4,zmin1, psnr_value4, xmin1,ymin1,beta4,mssim4)
fprintf(fileID,output4);
%% find minimum Niqe:

[zmin2,loc5] = min(H(:,8));
a5=loc5;
xmin2=H(a5,1);
ymin2=H(a5,2);
lamb5=H(a5,9);
psnr_value5=H(a5,8);
mssim5=H(a5,3);
output5=sprintf('Step: %d, min niqe= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta: %.6f, MSSIM: %.6f\n',a5,zmin2, psnr_value5, xmin2,ymin2,lamb5, mssim5)
fprintf(fileID,output5);
%% find minimum Piqe:

[zmin3,loc6] = min(H(:,10));
a6=loc6;
xmin3=H(a6,1);
ymin3=H(a6,2);
lamb6=H(a6,9);
psnr_value6=H(a6,8);
mssim6=H(a6,3);
output6=sprintf('Step: %d, min Piqe= %.8f, PSNR= %.8f \n alpha1= %.6f, alpha0= %.6f, beta: %.6f, MSSIM: %.6f\n',a6,zmin3, psnr_value6, xmin3,ymin3,lamb6, mssim6)
fprintf(fileID,output6);

fclose(fileID);
end