function compare

s=load('DATA/fft_100_2000_par');
t=s.result.I(1).factor;
plot(t(:,1),t(:,3),'.')
hold on
plot(t(:,1),t(:,5),'.')
hold on
plot(t(:,1),t(:,9),'.')