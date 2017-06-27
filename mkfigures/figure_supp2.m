clear;

g = load('../data/fig2data');

NV = [100, 10];
cols = ['k', 'b'];

figure;
% subplot(1,2,1);
hold;set(gca,'fontsize',15)
for ii = 1:2
    N = NV(ii);
    col = cols(ii);
    
    IND1 = find(g(:,1)==N);
    plot(g(IND1, 2), -g(IND1, 3),'-','linewidth',2, 'color', col);
end

xlabel('f_A');ylabel('-N\Gamma_3(q^*)');box on
axis([.2, .5, 0, 200])




figure;
set(gca,'fontsize',18)
for ii = 1:2
    subplot(1,2,ii);
    hold;

    N = NV(ii);
    col = cols(ii);
    
    IND1 = find(g(:,1)==N);
    plot(g(IND1, 2), g(IND1, 4),'-','linewidth',2, 'color', col);
%     plot(g(IND1, 2), g(IND1, 5),'-.','linewidth',2, 'color', col);
    plot(g(IND1, 2), g(IND1, 6),':','linewidth',2, 'color', col);
    plot(g(IND1, 2), g(IND1, 7),'--','linewidth',2, 'color', col);
    
    xlabel('f_A');ylabel('N\Gamma_4(q^*)');box on
    axis([.3, .5, 100, 550])
    legend('\Gamma_{4}(0,0)', '\Gamma_{4}(0,2)', '\Gamma_{4}(1,2)')
    title(strcat('N=', num2str(N)))
end

% legend('\Gamma_{4}(0,0)', '\Gamma_{4}(0,1)', '\Gamma_{4}(0,2)', '\Gamma_{4}(1,2)')