function plotting

s=load('DATA/he_10_10000_par_2020-05-27T053933');
%x=zeros(10000,2);
%for i=1:2
 %   for j=1:1000
%x(j,i)=s.result.I(i).factor(j,4)-s.result.I(i).noise_met(2);
 %   end
%end
%x
for i=8:8
    t=s.result.I(i);
    plot(t.factor(:,1),t.factor(:,4),'.')
    hold on
    plot(t.factor(:,1),t.noise_met(2),'.')
end
