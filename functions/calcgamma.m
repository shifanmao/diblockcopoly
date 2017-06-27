function [gam3,gam4]=calcgamma(N,FAV,NQ)
%% calcgamma.m :: This code calculates the coefficients in free energy expansion of block
% copolymers. Coefficients are evaluated by finding the cubic and quartic
% order vertex functions at dominant peak given in quadratic order fluctuations
% Usage: [gam3,gam4]=calcgamma(N,FAV,NQ)
% Parameters:
%   N, number of statistical steps of total chain
%   FAV, range of A-type monomer fractions
%   NQ, number of wavevector sets in calculating GAM4
% Return:
%   gam3, cubic order vertex constant
%   gam4, quartic order vertex function

% results to return
gam3=zeros(length(FAV),1);
gam4=zeros(length(FAV),NQ);

filename='data/gamdata';
if exist(filename,'file')
    data=dlmread(filename);
    for ii=1:length(FAV)
        FA = FAV(ii);
        FB = 1-FA;
        ind1 = find(abs(data(:,2)-FA)<1e-4 & abs(data(:,1)-N)<1e-2);
        ind2 = find(abs(data(:,2)-FB)<1e-4 & abs(data(:,1)-N)<1e-2);
        if ~isempty(ind1)
            fprintf('Step 2: Loading vertices at FA=%.2f, N=%.2e\n',FA,N)
            gam3(ii)=data(ind1,3)/N;
            gam4(ii,1:NQ)=data(ind1,4:3+NQ)/N;
        elseif ~isempty(ind2)
            fprintf('Step 2: Loading vertices at FA=%.2f, N=%.2e\n',FA,N)
            gam3(ii)=-data(ind2,3)/N;
            gam4(ii,1:NQ)=data(ind2,4:3+NQ)/N;
        else
            [gam3(ii),gam4(ii,1:NQ)]=gamma(N,FA,NQ);
        end
    end
else
    for ii=1:length(FAV)
        FA=FAV(ii);
        [gam3(ii),gam4(ii,1:NQ)]=gamma(N,FA,NQ);
    end
end
end

function [gam3,gam4]=gamma(N,FA,NQ)
    fprintf('Step 2: Calculating vertices at FA=%.2f, N=%.2e\n',FA,N)
    
    % wavevectors for Gamma4 calculations
    Q1=zeros(3,NQ);
    Q2=zeros(3,NQ);
    Q3=zeros(3,NQ);
    Q4=zeros(3,NQ);

    if NQ==1
        % case1 :: theta=pi
        Q1(1:3,1)=[1,0,0];
        Q2(1:3,1)=rotz(pi)*Q1(1:3,1);
        Q3(1:3,1)=-Q2(1:3,1);
        Q4(1:3,1)=-Q1(1:3,1);

    elseif NQ==4
        % case1 :: theta=pi
        Q1(1:3,1)=[1,0,0];
        Q2(1:3,1)=rotz(pi)*Q1(1:3,1);
        Q3(1:3,1)=-Q2(1:3,1);
        Q4(1:3,1)=-Q1(1:3,1);

        % case2 :: theta=pi/3
        Q1(1:3,2)=[1,0,0];
        Q2(1:3,2)=rotz(pi/3)*Q1(1:3,2);
        Q3(1:3,2)=-Q1(1:3,2);
        Q4(1:3,2)=-Q2(1:3,2);

        % case3 :: theta=pi/2
        Q1(1:3,3)=[1,0,0];
        Q2(1:3,3)=rotz(pi/2)*Q1(1:3,3);
        Q3(1:3,3)=-Q1(1:3,3);
        Q4(1:3,3)=-Q2(1:3,3);

        % case4 :: gam4(1,2)
        Q1(1:3,4)=[-1,0,1];
        Q2(1:3,4)=[-1,0,-1];
        Q3(1:3,4)=[1,1,0];
        Q4(1:3,4)=[1,-1,0];
    else
        Qrot=linspace(0,pi/2,NQ);
        for IQ=1:NQ
            Qvec=Qrot(IQ);
            Q1(1:3,IQ)=[1,0,0];
            Q2(1:3,IQ)=rotz(Qvec)*Q1(1:3,IQ);
            Q3(1:3,IQ)=-Q1(1:3,IQ);
            Q4(1:3,IQ)=-Q2(1:3,IQ);
        end
    end

    % calculate spinodal and critical wavelength
    [~,ks]=spinodal(N,FA);

    % calculate free energy coefficients
    gam3=gamma3(N,FA,ks);
    gam4=zeros(1,NQ);

    for IQ=1:NQ
        K1=Q1(:,IQ);
        K2=Q2(:,IQ);
        K3=Q3(:,IQ);
        K4=Q4(:,IQ);
        gam4(IQ)=gamma4(N,FA,ks,K1,K2,K3,K4);
    end
end