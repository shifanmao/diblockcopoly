% Use K-S test to verify Lorentzian form
% of structure factor
clear;
% close all
addpath('../functions/')
addpath('../chainstats/')
addpath('../chainstats/')
addpath('../chainstats/eigcalc/')

N=1e4;
FA=0.5;
NK=1000;
NC=75;

RM=sqrt(r2(N));
% RM=1;
% k=logspace(-2,1,NK)/RM;
k=linspace(0.001,10,NK)/RM;
CHI=0;

% begin calculation
s=1./gamma2(N,FA,k,CHI)';
[~,ks,d2gam2]=spinodal(N,FA);

% find peak
kmax=find(s==max(s));
smax=1./gamma2(N,FA,ks,0);
alpha2 = (1/2)*RM^(-2)*d2gam2;

% find cdf
k0 = max([kmax-NC,1]);
kf = min([kmax+NC,length(k)]);
kc = k(k0:kf);
sc = s(k0:kf);
sfit = 1./(1/smax+alpha2.*(kc-ks).^2*RM^2);

% make a plot
figure;hold;
col='k';
set(gca,'fontsize',18)
area([kc(1),kc(end)],[.01,.01])
plot(k*RM,s/N,'-','linewidth',2,'color',col);
plot(ks*RM,smax/N,'o','color',col);
plot(kc*RM,sfit/N,'.-','linewidth',2,'color',col);
ylim([0,.1]);
xlabel('qR');ylabel('$\tilde{\psi}(q)\tilde{\psi}(-q)/N$','Interpreter','latex')

% plot cdf
s1 = sc/sum(sc)/(kc(2)-kc(1))/RM;
s2 = sfit/sum(sfit)/(kc(2)-kc(1))/RM;

ts1 = cumsum(s1)*(kc(2)-kc(1))*RM;
ts2 = cumsum(s2)*(kc(2)-kc(1))*RM;

figure;hold
plot(kc,ts1,kc,ts2);

% sample from cdf
nr = 1000;
kr1 = zeros(1,nr);
for ii = 1:nr
    r = rand();
    i0 = find(ts1>=rand());i0 = i0(1);
    kr1(ii) = kc(i0);
end

kr2 = zeros(1,nr);
for ii = 1:nr
    r = rand();
    i0 = find(ts2>=rand());i0 = i0(1);
    kr2(ii) = kc(i0);
end

% figure;hold
% plot(kc,s1);
% [f,x]=hist(kr1,20);%# create histogram from a normal distribution.
% bar(x,f/trapz(x,f));hold on
% 
% figure;hold
% plot(kc,s2);
% [f,x]=hist(kr2,20);%# create histogram from a normal distribution.
% bar(x,f/trapz(x,f));hold on

% Use KS test
kstest2(kr1,kr2)
% OUTPUT: if 1, rejects null hypothesis that kr1&kr2 sampled from same distr.
% i.e. if 1, not from same distr.
%      if 0, from same distr.