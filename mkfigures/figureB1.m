cd ../

NV=logspace(-1,3,51);  % number of statistical steps of total chain
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
[chis,ks,d2gam2]=spinodal(NV,0.5);

figure;hold;set(gca,'fontsize',20);
plot(NV,ks,'k-','linewidth',2);
x = logspace(-1,-0.5,20);y = power(x,-1);plot(x,5*y,'k--','linewidth',1.5)
x = logspace(2,3,20);y = power(x,-1/2);  plot(x,3.5*y,'k--','linewidth',1.5)
set(gca,'xscale','log');set(gca,'yscale','log');
xlabel('N');ylabel('2l_Pq^*');box on

figure;hold;set(gca,'fontsize',20);
plot(NV,d2gam2/2,'k-','linewidth',2);
x = logspace(-1,-0.5,20);y = power(x,1);plot(x,0.25*y,'k--','linewidth',1.5)
x = logspace(2,3,20);y = power(x,0);  plot(x,0.9*y,'k--','linewidth',1.5)
set(gca,'xscale','log');set(gca,'yscale','log');
xlabel('N');
ylabel('$c=\frac{1}{2} \left.\frac{\partial^2\Gamma_{2}}{\partial(2l_Pq)^{2}} \right|_{q=q^{*}}$','interpreter','latex');
box on

%%% figure B2 %%%
figure;hold;set(gca,'fontsize',20)
N = 100;
col = 'k';

NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
plot(FAV,-gam3*N,'-','linewidth',2,'color',col);
xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)');box on

N = 10;
col = 'b';

NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
plot(FAV,-gam3*N,'-','linewidth',2,'color',col);
xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)');box on

%%% figure B3 %%%
figure;hold;set(gca,'fontsize',20)
N = 100;
col = 'k';

NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
plot(FAV,gam4*N,'-','linewidth',2,'color',col);
xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*,\theta=0)');box on

N = 10;
col = 'b';

NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
plot(FAV,gam4*N,'-','linewidth',2,'color',col);
xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*,\theta=0)');box on

cd mkfigures/