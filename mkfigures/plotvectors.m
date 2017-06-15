function hq = plotvectors(X, Y, V, arrowlength)
% X: fA
% Y: phiP

X1 = zeros(length(Y), length(Y));
X2 = zeros(length(Y), length(Y));
V1 = X1;V2 = X2;

% triangular mesh
lenFAVV = zeros(length(Y), 1);
for jj = 1:length(Y)
    lenFAVV(jj) = jj;
end

for jj = 1:length(Y)
    if mod(jj,5)==1 || jj==length(Y)
        Xnew = linspace((1-Y(jj))/sqrt(3), ...
            2/sqrt(3)-(1-Y(jj))/sqrt(3), lenFAVV(jj)) / Y(jj);
        for ii = 1:5:lenFAVV(jj)
            X1(ii,jj) = Xnew(ii)*Y(jj);
            X2(ii,jj) = 1-Y(jj);

            V1(ii,jj) = V(ii, jj, 1) - V(ii, jj, 2);
            V2(ii,jj) = sqrt(3)*(-V(ii, jj, 1) - V(ii, jj, 2));

            V1(ii,jj) = V1(ii,jj)*arrowlength;
            V2(ii,jj) = V2(ii,jj)*arrowlength;
            
            if ~(jj==length(Y) && ii==1) || (jj==length(Y) && ii==lenFAVV(jj))
                plot3(X1(ii,jj), X2(ii,jj), 1e5, 'k.', 'markersize', 15)
            end
        end
    end
end

Z = zeros(length(Y), length(Y))+1e5;
W = zeros(length(Y), length(Y));

hq = quiver3(X1,X2,Z,V1,V2,W,'color','k', 'linewidth', 1.5,'AutoScale','off');