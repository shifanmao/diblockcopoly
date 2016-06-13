function valeq=S4_case2_int(FA,R1,R12,R3)
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
        valeq=FA^2*FB*(3*R12+2*FA*R1-FA*R12+6)/12;
    elseif nzeros==2
        if max(abs([R1,R3]))<MIN   % 0~R1~R3
            valeq=(1/2).*((-1)+FA).*R12.^(-3).*(2+(-2).*exp(1).^(FA.*R12)+FA.*R12.*( ...
                        2+FA.*R12));
        elseif max(abs([R1,R12]))<MIN  % 0~R1~R12
            valeq=(1/2).*exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)) ...
                        .*R3.^(-4).*((-2)+2.*exp(1).^(FA.*R3)+(-1).*FA.*R3.*(2+FA.*R3));
        elseif max(abs([R12,R3]))<MIN % 0~R12~R3
            valeq=(1/2).*((-1)+FA).*R1.^(-3).*(2+(-2).*exp(1).^(FA.*R1)+FA.*R1.*(2+ ...
                        FA.*R1));
        else
            error('not an option')
        end
    elseif nzeros==1
        if abs(R1-R12)<MIN && abs(R3)<MIN
            valeq=(1-FA)*(-2*expl(3,FA*R1)+FA*R1*expl(2,FA*R1))/(R1^3);
        elseif abs(R12-R3)<MIN && abs(R1)<MIN
            valeq=(2*(expl(1,R3)*expl(3,-FA*R3)+2*coshl(4,FA*R3))...
                +FA*R3*(expl(1,R3)*expl(2,-FA*R3)-2*sinhl(3,FA*R3)))/(R3^4);
        elseif abs(R1-R3)<MIN && abs(R12)<MIN
            valeq=(2*(expl(1,R3)*expl(3,-FA*R3)+2*coshl(4,FA*R3))...
                +FA*R3*(expl(1,R3)*expl(2,-FA*R3)-2*sinhl(3,FA*R3)))/(R3^4);
        elseif abs(R1)<MIN
            valeq=(exp(FB*R3)-1)*(R12^2*expl(3,FA*R3)-R3^2*expl(3,FA*R12))/(R3^3*R12^2*(R3-R12));
        elseif abs(R12)<MIN   
            valeq=(exp(FB*R3)-1)*(R1^2*expl(3,FA*R3)-R3^2*expl(3,FA*R1))/(R3^3*R1^2*(R3-R1));
        elseif abs(R3)<MIN
            valeq=FB*(R1^2*expl(3,FA*R12)-R12^2*expl(3,FA*R1))/(R1^2*R12^2*(R12-R1));
        else
            error('not an option')
        end
    else
        if (abs(R1-R12)<MIN && abs(R12-R3)<MIN)
           valeq=(1/2).*exp(1).^((-1).*FA.*R1).*(exp(1).^R1+(-1).*exp(1).^(FA.*R1)) ...
                      .*R1.^(-4).*((-2)+exp(1).^(FA.*R1).*(2+FA.*R1.*((-2)+FA.*R1)));
        elseif abs(R1-R12)<MIN
           valeq=exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)).*R1.^( ...
                      -2).*(R1+(-1).*R3).^(-2).*R3.^(-2).*((-1).*R1.^2+exp(1).^(FA.*R3) ...
                      .*R1.^2+2.*R1.*R3+(-1).*R3.^2+exp(1).^(FA.*R1).*R3.*(R1.*((-2)+ ...
                      FA.*(R1+(-1).*R3))+R3));
        elseif abs(R1-R3)<MIN
           valeq=exp(1).^((-1).*FA.*R1).*(exp(1).^R1+(-1).*exp(1).^(FA.*R1)).*R1.^( ...
                      -3).*(R1+(-1).*R12).^(-2).*R12.^(-1).*(((-1)+exp(1).^(FA.*R1)).* ...
                      R12.^2+R1.^2.*((-1)+exp(1).^(FA.*R12)+exp(1).^(FA.*R1).*FA.*R12)+ ...
                      R1.*R12.*(2+(-1).*exp(1).^(FA.*R1).*(2+FA.*R12)));
        elseif abs(R12-R3)<MIN
           valeq=exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)).*R1.^( ...
                      -1).*(R1+(-1).*R3).^(-2).*R3.^(-3).*((-1).*R1.^2+2.*R1.*R3+(-1).* ...
                      R3.^2+exp(1).^(FA.*R1).*R3.^2+exp(1).^(FA.*R3).*R1.*(R1+(-1).*FA.* ...
                      R1.*R3+R3.*((-2)+FA.*R3)));
        else
           valeq=exp(1).^((-1).*FA.*R3).*(exp(1).^R3+(-1).*exp(1).^(FA.*R3)).*R1.^( ...
                      -1).*(R1+(-1).*R12).^(-1).*R12.^(-1).*(R1+(-1).*R3).^(-1).*(R12+( ...
                      -1).*R3).^(-1).*R3.^(-2).*(((-1)+exp(1).^(FA.*R1)).*R12.*(R12+(-1) ...
                      .*R3).*R3+R1.^2.*(((-1)+exp(1).^(FA.*R3)).*R12+R3+(-1).*exp(1).^( ...
                      FA.*R12).*R3)+R1.*((1+(-1).*exp(1).^(FA.*R3)).*R12.^2+((-1)+exp(1) ...
                      .^(FA.*R12)).*R3.^2));
        end 
    end
    
    if isnan(valeq)
        disp(num2str(min([R1,R12,R3])))
        disp([' E1=',num2str(R1),' // E12=',num2str(R12),' // E3=',num2str(R3)])
        error(' Error in case 2 :: valeq NaN ')
    end
end