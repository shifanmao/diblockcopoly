function valeq=S4_case3_int(FA,R1,R12,R3)
% on same monomer
%     if min([R1,R12,R3])<-700
    if min(real([R1,R12,R3]))<-500
        valeq=0;
        return
    end
    MIN=(10^-4);    
    FB=1-FA;
    vec=[R1,R12,R3];
    nzeros=sum(abs(vec)<MIN);
    if nzeros==3
        valeq=FA^2*FB^2*(R3+E12+FA*E1-FA*E3+3)/12;
    elseif nzeros==2
        if max(abs([R1,R3]))<MIN   % 0~R1~R3
            valeq=exp(1).^((-1).*FA.*R12).*R12.^(-4).*((-1)+exp(1).^(FA.*R12)+(-1).* ...
                       FA.*R12).*(exp(1).^R12+exp(1).^(FA.*R12).*((-1)+((-1)+FA).*R12));
        elseif max(abs([R1,R12]))<MIN  % 0~R1~R12
            valeq=(-1/2).*FA.^2.*R3.^(-2).*(1+(-1).*exp(1).^(R3+(-1).*FA.*R3)+R3+( ...
                        -1).*FA.*R3);
        elseif max(abs([R12,R3]))<MIN % 0~R12~R3
            valeq=(1/2).*((-1)+FA).^2.*R1.^(-2).*((-1)+exp(1).^(FA.*R1)+(-1).*FA.*R1);
        else
            error('not an option')
        end
    elseif nzeros==1
        if abs(R1-R12)<MIN && abs(R3)<MIN
            valeq=(-FB*FA*expl(2,R1*FA)*R1^2+...
                 (FA*expl(3,R1)+(2*FB-1)*expl(3,R1*FA))*R1+...
                 expl(4,FB*R1)-expl(4,R1)+expl(4,R1*FA) )/(R1^4);
        elseif abs(R12-R3)<MIN && abs(R1)<MIN
            valeq=(-FA*FB*expl(2,R3*FB)*R3^2+...
                 (FB*expl(3,R3)+(2*FA-1)*expl(3,R3*FB))*R3+...
                 expl(4,FA*R3)-expl(4,R3)+expl(4,R3*FB) )/(R3^4);
        elseif abs(R1-R3)<MIN && abs(R12)<MIN
            valeq=( R3*(-FB*expl(3,FA*R3)-FA*expl(3,R3*FB))+...
                expl(4,R3)-expl(4,FA*R3)-expl(4,R3*FB) )/(R3^4);
        elseif abs(R1)<MIN
            valeq=( R3*expl(2,R12*FB)-R12*expl(2,R3*FB) )*expl(2,FA*R12)/(R3*R12^3*(R12-R3));
        elseif abs(R12)<MIN   
            valeq=expl(2,FA*R1)*expl(2,R3*FB)/(R1^2*R3^2);
        elseif abs(R3)<MIN
            valeq=( R1*expl(2,R12*FA)-R12*expl(2,R1*FA) )*expl(2,FB*R12)/(R1*R12^3*(R12-R1));
        else
            error('not an option')
        end
    else
        if (abs(R1-R12)<MIN && abs(R12-R3)<MIN)
            valeq=exp(1).^((-1).*FA.*R1).*R1.^(-4).*(exp(1).^(FA.*R1)+exp(1).^R1.*(( ...
                      -1)+R1+(-1).*FA.*R1)).*(1+exp(1).^(FA.*R1).*((-1)+FA.*R1));
        elseif abs(R1-R12)<MIN
            valeq=exp(1).^((-1).*FA.*(R1+R3)).*R1.^(-3).*(1+exp(1).^(FA.*R1).*((-1)+ ...
                      FA.*R1)).*(R1+(-1).*R3).^(-1).*R3.^(-1).*((-1).*exp(1).^(FA.*R1+ ...
                      R3).*R1+exp(1).^(FA.*(R1+R3)).*(R1+(-1).*R3)+exp(1).^(R1+FA.*R3).* ...
                      R3);
        elseif abs(R1-R3)<MIN
            valeq=exp(1).^((-1).*FA.*(R1+R12)).*R1.^(-2).*(R1+(-1).*R12).^(-2).* ...
                  R12.^(-2).*((-1).*exp(1).^(FA.*R1+R12).*R1+exp(1).^(FA.*(R1+R12)) ...
                  .*(R1+(-1).*R12)+exp(1).^(R1+FA.*R12).*R12).*(R1+(-1).*exp(1).^( ...
                  FA.*R12).*R1+((-1)+exp(1).^(FA.*R1)).*R12);
        elseif abs(R12-R3)<MIN
            valeq=exp(1).^((-1).*FA.*R12).*R1.^(-1).*R12.^(-3).*((-1).*R1+R12).^(-1) ...
                  .*(((-1)+exp(1).^(FA.*R12)).*R1+R12+(-1).*exp(1).^(FA.*R1).*R12).* ...
                  (exp(1).^(FA.*R12)+exp(1).^R12.*((-1)+R12+(-1).*FA.*R12));
        else
            valeq=exp(1).^((-1).*FA.*(R12+R3)).*R1.^(-1).*(R1+(-1).*R12).^(-1).* ...
              R12.^(-2).*(((-1)+exp(1).^(FA.*R12)).*R1+R12+(-1).*exp(1).^(FA.* ...
              R1).*R12).*(R12+(-1).*R3).^(-1).*R3.^(-1).*(exp(1).^(FA.*R12+R3).* ...
              R12+(-1).*exp(1).^(R12+FA.*R3).*R3+exp(1).^(FA.*(R12+R3)).*((-1).* ...
              R12+R3));
        end
    end
    
    if isnan(valeq)
        disp([' E1=',num2str(R1),' // E12=',num2str(R12),' // E3=',num2str(R3)])
        error(' Error in case 3 :: valeq NaN ')
    end
end