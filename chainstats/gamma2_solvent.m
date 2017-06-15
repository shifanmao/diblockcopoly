function [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP)

%% Inputs
PHIS = 1-PHIP;   % Fraction of solvent

%% Start Calculation
EIG = zeros(length(KV),2);
EIGV = zeros(length(KV),4);

for ii = 1:length(KV)
  k = KV(ii);  % wavevector
  
  if PHIP == 1
      if FA > 0 && FA < 1
          %%% CALCULATE S FUNCTIONS %%%
          s2inv = s2inverse(N,FA,k);
          SAAINV = s2inv(1,1);
          SABINV = s2inv(1,2);
          SBBINV = s2inv(2,2);

          %%% CALCULATE GAM FUNCTIONS %%%
          GAMAA = -2*CHI(1,1)+N*SAAINV./PHIP;
          GAMBB = -2*CHI(2,2)+N*SBBINV./PHIP;
          GAMAB = (CHI(1,2)-CHI(1,1)-CHI(2,2))+N*SABINV/PHIP;
          GAMBA = GAMAB;
          GAM = [GAMAA, GAMAB; GAMBA, GAMBB];

          EIG1 = GAM(1,1) + GAM(2,2) - GAM(1,2) - GAM(2,1);

          EIGV(ii, 1:4) = [1; -1; 1; 1]/sqrt(2);
          EIG(ii, 1:2) = [EIG1, 1e10];
      else
          EIGV(ii, 1:4) = [NaN; NaN; NaN; NaN];
          EIG(ii, 1:2) = [NaN, NaN];
      end
      
  elseif PHIP == 0
      EIGV(ii, 1:4) = [NaN; NaN; NaN; NaN];
      EIG(ii, 1:2) = [NaN, NaN];

  else
      
      if FA > 0 && FA < 1
          %%% CALCULATE S FUNCTIONS %%%
          s2inv = s2inverse(N,FA,k);
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

      elseif FA == 1
          %%% CALCULATE S FUNCTIONS %%%
          s2 = s2wlc(N,FA,k);
          SAAINV = 1/s2(1,1);

          %%% CALCULATE GAM FUNCTIONS %%%
          GAMAA = -2*CHI(1,1)+N*SAAINV./PHIP+1/PHIS;

          %%% CALCULATE EIGENMODES %%%
          EIGV(ii, 1:4) = [-1; 0; 0; -1];
          EIG(ii, 1:2) = [GAMAA, 1e10];

      elseif FA == 0
          %%% CALCULATE S FUNCTIONS %%%
          s2 = s2wlc(N,FA,k);
          SBBINV = 1/s2(2,2);

          %%% CALCULATE GAM FUNCTIONS %%%
          GAMBB = -2*CHI(2,2)+N*SBBINV./PHIP+1/PHIS;

          %%% CALCULATE EIGENMODES %%%
          EIGV(ii, 1:4) = [0; -1; -1; 0];
          EIG(ii, 1:2) = [GAMBB, 1e10];
      end
  end

end

end