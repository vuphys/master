function plotting

s=load('DATA/fft_100_2000_par');
for i=1:100
    t=s.result.I(i).factor;
    plot(t(:,1),t(:,9),'.')
    hold on
end