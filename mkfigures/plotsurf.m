function [X1, X2, V1] = plotsurf(X, Y, V)


%% first surface
X1 = NaN(length(X), length(Y));
X2 = NaN(length(X), length(Y));
V1 = NaN(length(X), length(Y));

% triangular mesh
lenFAVV = zeros(length(Y), 1);
for jj = 1:length(Y)
    lenFAVV(jj) = jj;
end

for jj = 1:length(Y)
    Xnew = linspace((1-Y(jj))/sqrt(3), ...
                     2/sqrt(3)-(1-Y(jj))/sqrt(3), lenFAVV(jj));
    for ii = 1:lenFAVV(jj)
        X1(ii,jj) = Xnew(ii);
    	X2(ii,jj) = 1-Y(jj);
        V1(ii,jj) = V(ii, jj);
    end
end

surface(X1, X2, V1, 'FaceColor', 'interp','EdgeColor','none');


%% plot again
X1 = NaN(length(X), length(Y));
X2 = NaN(length(X), length(Y));
V1 = NaN(length(X), length(Y));

% triangular mesh
lenFAVV = zeros(length(Y), 1);
for jj = 1:length(Y)
    lenFAVV(jj) = jj;
end

for jj = 1:length(Y)
    Xnew = linspace((1-Y(jj))/sqrt(3), ...
                     2/sqrt(3)-(1-Y(jj))/sqrt(3), lenFAVV(jj));
    Xnew = fliplr(Xnew);
    for ii = 1:lenFAVV(jj)
        X1(ii,jj) = Xnew(ii);
    	X2(ii,jj) = 1-Y(jj);
        V1(ii,jj) = V(ii, jj);
    end
end

surface(X1, X2, V1, 'FaceColor', 'interp','EdgeColor','none');
view([0,90])
% contour(X1, X2, V1, 'ShowText', 'on');

% plot3(X1, X2, V1, 'ko');