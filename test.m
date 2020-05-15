
function test
ffcdata=load('GroundTruth');
fname=ffcdata.data(:,:,1,1);
%for i=1:3
%t=rand(1,1);
rng(0);
%L=fname+randn(size(fname))*0.1;
K= imnoise(fname,'Gaussian',0,1);
%LL=psnr(L,fname);
KK=psnr(K,fname);
a=vifvec(fname,K);
sprintf('psnr: %.4f metric: %.8f',KK,a)
imshow([fname,K]);
%end
end