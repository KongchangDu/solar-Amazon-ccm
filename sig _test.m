load('solar_amazon1.mat')
a=[];
for i =1:200
r1=mean(amazon)+ std(amazon) .* randn(length(amazon),1);
eimf=eemd(r1,0.6,2000,50,6);
val=corr(sunspots_mean,eimf(:,4));
a=[a val];
end