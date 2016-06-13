function [firstD,secondD]=WignerD(q1,q2,q3,ORDL)
    % Find unit vectors
    Q1_n=q1/norm(q1,2);
    Q2_n=q2/norm(q2,2);
    Q3_n=q3/norm(q3,2);
    
    % Find Euler angles
    aligned=10^-13; % angle in radians betwene to vectors before I assume they are the same
    if norm(q2) < aligned
        % if q2=0 we may as well point it perpendicular to the other two
        cosB1=0;
        cosB2=0;
        alpha1=dot(Q1_n,Q3_n);
    elseif abs(Q2_n-Q1_n) < aligned
        cosB1=1;
        alpha1=0;
        cosB2=dot(Q2_n,Q3_n);
    elseif abs(Q2_n-Q3_n) < aligned
        cosB1=dot(Q2_n,Q1_n);
        alpha1=0;
        cosB2=1;
    else
        cosB1=dot(Q2_n,Q1_n);
        cosB2=dot(Q2_n,Q3_n);
        v1=cross(Q2_n,Q1_n)/norm(cross(Q2_n,Q1_n));
        v2=cross(Q2_n,Q3_n)/norm(cross(Q2_n,Q3_n));
        alpha1=acos(dot(v1,v2));
    end 
    % Calculate Wigner D matrices
    firstD=WignerD_lm0(ORDL,alpha1,cosB1);
    secondD=WignerD_lm0(ORDL,0,cosB2);
end
function out=WignerD_lm0(ORDL,alpha,cosB)
% Wigner D Matrix when second L index is zero
% Returns a ORDL x ORDL matrix with indicies (lam+1,M+1)
out=zeros(ORDL,ORDL)*NaN;

for L=0:(ORDL-1)
    M=(0:L);
    out(L+1,1:(L+1))=legendre(L,cosB)'.*...
                     sqrt(factorial(L-M)./factorial(L+M))...
                     .*exp(1i*M*alpha);
end
end