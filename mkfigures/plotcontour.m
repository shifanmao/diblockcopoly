function plotcontour(X, Y, V)

X1 = zeros(length(X), length(Y));
X2 = zeros(length(X), length(Y));
V1 = zeros(length(X), length(Y));


for jj = 1:length(Y)
    Xnew = linspace(Y(jj)/sqrt(3), 2/sqrt(3)-Y(jj)/sqrt(3), length(X));
    for ii = 1:length(X)
        
        X1(ii,jj) = Xnew(ii);
    	X2(ii,jj) = Y(jj);
        V1(ii,jj) = V(ii, jj);
    end
end

% contourf(X1, X2, V1, 10, 'ShowText', 'on');
contourf(X1, X2, V1, 'ShowText', 'on');