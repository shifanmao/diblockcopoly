% function plotvectors(X, Y, V)
% 
% X1 = zeros(length(X), length(Y));
% X2 = zeros(length(X), length(Y));
% V1 = X1;V2 = X2;
% for ii = 1:length(X)
%     for jj = 1:length(Y)
%         X1(ii,jj) = X(ii);
%     	X2(ii,jj) = Y(jj);
%         V1(ii,jj) = V(ii, jj, 1);
%         V2(ii,jj) = V(ii, jj, 2);
%     end
% end
% 
% % figure;hold
% hq = quiver(X1,X2,V1,V2,'color','k', 'linewidth', 1.5,'AutoScale','off');
% hq = quiver(X1,X2,-V1,-V2,'color','k', 'linewidth', 1.5,'AutoScale','off');
% 

function plotvectors(X, Y, V)


X1 = zeros(length(X), length(Y));
X2 = zeros(length(X), length(Y));
V1 = X1;V2 = X2;

length(X)

for jj = 1:length(Y)
    if mod(jj,3)==0
        Xnew = linspace(Y(jj)/sqrt(3), ...
            2/sqrt(3)-Y(jj)/sqrt(3), length(X)) / Y(jj);
        for ii = 1:length(X)
            if mod(ii,1)==0
                X1(ii,jj) = Xnew(ii)*Y(jj);
                X2(ii,jj) = Y(jj);

                V1(ii,jj) = V(ii, jj, 1) - V(ii, jj, 2);
                V2(ii,jj) = sqrt(3)*(-V(ii, jj, 1) - V(ii, jj, 2));
            end
        end
    end
end

% figure;hold
% hq = quiver(X1,X2,V1,V2,'color','k', 'linewidth', 1.5,'AutoScale','off');
% hq = quiver(X1,X2,-V1,-V2,'color','k', 'linewidth', 1.5,'AutoScale','off');

Z = zeros(length(X), length(Y))+1e5;
W = zeros(length(X), length(Y));

hq = quiver3(X1,X2,Z,V1,V2,W,'color','k', 'linewidth', 1.5,'AutoScale','off');