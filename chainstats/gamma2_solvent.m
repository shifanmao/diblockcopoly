function [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP)

%% Inputs
PHIS = 1-PHIP;   % Fraction of solvent

%% Start Calculation
EIG = zeros(length(KV),2);
EIGV = zeros(length(KV),4);

for ii = 1:length(KV)
  k = KV(ii);  % wavevector
  
  %%% CALCULATE S FUNCTIONS %%%
  s2inv=s2inverse(N,FA,k);
  SAAINV = s2inv(1,1);
  SABINV = s2inv(1,2);
  SBBINV = s2inv(2,2);
  
  %%% CALCULATE GAM FUNCTIONS %%%
  GAMAA = -2*CHI(1,1)+N*SAAINV./PHIP+1/PHIS;
  GAMBB = -2*CHI(2,2)+N*SBBINV./PHIP+1/PHIS;
  GAMAB = (CHI(1,2)-CHI(1,1)-CHI(2,2))+N*SABINV/PHIP+1/PHIS;
  GAMBA = GAMAB;
  GAM = [GAMAA, GAMAB; GAMBA, GAMBB];

  %%% CALCULATE EIGENMODES %%%
  [EIV,EI] = eig(GAM);
  EI = [EI(1,1),EI(2,2)];
  
  sign = 1;
  [~,ind] = sort(EI,'ascend');ind1 = ind(1);ind2 = ind(2);
  EIG(ii, 1:2) = [EI(ind1); EI(ind2)];
  EIGV(ii, 1:4) = [EIV(:, ind1)*sign; EIV(:, ind2)*sign];
  
%   if EIV(1,1)*EIV(2,1) < 0
%       if EIV(1,1) < 0
%           sign = -1;
%       end
%       ind1 = 1;ind2 = 2;
%   else
%       if EIV(2,1) < 0
%           sign = -1;
%       end
%       ind1 = 2;ind2 = 1;
%   end
%   EIG(ii, 1:2) = [EI(ind1); EI(ind2)];
%   EIGV(ii, 1:4) = [EIV(:, ind1)*sign; EIV(:, ind2)*sign];
end

end