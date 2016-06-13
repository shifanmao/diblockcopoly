function coeffinder
% generates coefficients in free energy expansion
% F = B*Gam2*phi^2+C*Gam3*phi^3+...
%       (D1*Gam4(T=0)+D2*Gam4(T=pi/3)+D3*Gam4(T=pi/2)+D4*Gam4(1,2))*phi^4
clear;

% LAM PHASE
Q1=[1,0,0]';
Q2=-Q1;
Q=[Q1,Q2];
[B,C,D1,D2,D3,D4]=findcoef(Q);
disp(['LAM PHASE // B=',num2str(B),'/ C=',num2str(C),'/ D1=',num2str(D1),'/ D2=',num2str(D2),'/ D3=',num2str(D3)])

% CYL PHASE (tetragonally)
clear;

Q1=[1,0,0]';
Q2=[0,1,0]';
Q3=-Q1;
Q4=-Q2;

Q=[Q1,Q2,Q3,Q4];
[B,C,D1,D2,D3,D4]=findcoef(Q);
disp(['CYL(TETRA) PHASE // B=',num2str(B),'/ C=',num2str(C),...
    '/ D1=',num2str(D1),'/ D2=',num2str(D2),'/ D3=',num2str(D3),'/ D4=',num2str(D4)])

% CYL PHASE (hexagonal)
clear;

Q1=[1,0,0]';
Q2=[-1/2,sqrt(3)/2,0]';
Q3=[-1/2,-sqrt(3)/2,0]';
Q4=-Q1;
Q5=-Q2;
Q6=-Q3;

Q=[Q1,Q2,Q3,Q4,Q5,Q6];
[B,C,D1,D2,D3,D4]=findcoef(Q);
disp(['CYL(HEXA) PHASE // B=',num2str(B),'/ C=',num2str(C),...
    '/ D1=',num2str(D1),'/ D2=',num2str(D2),'/ D3=',num2str(D3),'/ D4=',num2str(D4)])

% BCC PHASE
clear;

Q1=[1,1,0]';
Q2=[-1,1,0]';
Q3=[0,1,1]';
Q4=[0,1,-1]';
Q5=[1,0,1]';
Q6=[1,0,-1]';
Q1=Q1./norm(Q1);
Q2=Q2./norm(Q2);
Q3=Q3./norm(Q3);
Q4=Q4./norm(Q4);
Q5=Q5./norm(Q5);
Q6=Q6./norm(Q6);

Q7=-Q1;
Q8=-Q2;
Q9=-Q3;
Q10=-Q4;
Q11=-Q5;
Q12=-Q6;

Q=[Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12];
[B,C,D1,D2,D3,D4]=findcoef(Q);
disp(['BCC PHASE // B=',num2str(B),'/ C=',num2str(C),...
    '/ D1=',num2str(D1),'/ D2=',num2str(D2),'/ D3=',num2str(D3),'/ D4=',num2str(D4)])
end

function [B,C,D1,D2,D3,D4]=findcoef(Q)
    [~,NQ]=size(Q);
    B=0;
    C=0;     %theta=pi/3
    D1=0;   %theta=0
    D2=0;   %theta=pi/3
    D3=0;   %theta=pi/2
    D4=0;   %gam4(1,2)

    for j1=1:NQ
        for j2=1:NQ
            K1=Q(:,j1);
            K2=Q(:,j2);

            if norm(K1+K2)<1e-10
                B=B+1;
            end
        end
    end

    for j1=1:NQ
        for j2=1:NQ
            for j3=1:NQ
                K1=Q(:,j1);
                K2=Q(:,j2);
                K3=Q(:,j3);

                if norm(K1+K2+K3)<1e-10
                    if abs(abs(dot(K1,K2))-cos(pi/3))<1e-10
                        C=C+1;
                    else
                        disp(['something wrong',num2str(abs(dot(K1,K2)))])
                    end
                end
            end
        end
    end

    for j1=1:NQ
        for j2=1:NQ
            for j3=1:NQ
                for j4=1:NQ
                    K1=Q(:,j1);
                    K2=Q(:,j2);
                    K3=Q(:,j3);
                    K4=Q(:,j4);

                    if norm(K1+K2+K3+K4)<1e-10
                        if dot(cross(K1,K2),K3)==0  % on same plane
                            if abs(dot(K1,K2)-(-1))<1e-10   % K1=-K2
                                if abs(dot(K1,K3)-(-1))<1e-10   %K1=-K3
                                    if abs(dot(K1,K4)-1)<1e-10    
                                        D1=D1+1;
                                    elseif abs(abs(dot(K1,K3))-cos(pi/3))<1e-10
                                        D2=D2+1;
                                    elseif abs(abs(dot(K1,K3))-0)<1e-10
                                        D3=D3+1;
                                    else
                                        disp(['here something wrong',num2str(dot(K1,K3))])
                                    end
                                elseif abs(dot(K1,K3)-1)<1e-10    
                                    D1=D1+1;
                                elseif abs(abs(dot(K1,K3))-cos(pi/3))<1e-10
                                    D2=D2+1;
                                elseif abs(abs(dot(K1,K3))-0)<1e-10
                                    D3=D3+1;
                                else
                                    disp(['here something wrong',num2str(dot(K1,K3))])
                                end
                            elseif abs(dot(K1,K2)-1)<1e-10    
                                D1=D1+1;
                            elseif abs(abs(dot(K1,K2))-cos(pi/3))<1e-10
                                D2=D2+1;
                            elseif abs(abs(dot(K1,K2))-0)<1e-10
                                D3=D3+1;
                            else
                                disp(['something wrong',num2str(abs(dot(K1,K2)))])
                            end
                        else
                            p=[power(norm(K1+K2)/norm(K1),2),...
                            power(norm(K1+K3)/norm(K1),2),...
                            power(norm(K1+K4)/norm(K1),2)];
                            if ((length(find(abs(p-1)<1e-2))==2) &&...
                                    length(find(abs(p-2)<1e-2))==1 )
                                D4=D4+1;
                                [K1,K2,K3,K4];
                            else
                                disp(['something wrong  ',num2str(p)])
                            end
                        end
                    end
                end
            end
        end
    end
end