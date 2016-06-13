function valeq=S4_case1_int(NM,E1,E12,E3)
MIN=(10^-3)/NM;
vec=[E1,E12,E3];
nzeros=sum(abs(vec)<MIN);
if nzeros==3
    valeq=NM^4*(NM*(E1+E12+E3)+5)/120;
elseif nzeros==2
    [~,index]=max(abs(vec));
    others=vec; others(index)=[];
    valeq=chicken(others(1),others(2),vec(index),NM);
elseif nzeros==1
    [~,index]=min(abs(vec));
    others=vec; others(index)=[];
    if abs(others(1)-others(2))<MIN
        val=0.5*(others(1)+others(2));
        valeq=(-6*expl(3,NM*val)+2*NM*val*expl(2,NM*val))/(2*val^4);
    else
        valeq=(expl(4,NM*others(2))/(others(2)^3)...
            -expl(4,NM*others(1))/(others(1)^3))...
            /(others(2)-others(1));
    end
else
    dif=abs([E12-E3,E3-E1,E1-E12]);
    if max(dif)<MIN
        val=mean([E1,E12,E3]);
        valeq=(1/2)*val^(-4)*((-2)*(3+val*NM)+exp(val*NM)*(6+val*NM*(( ...
            -4)+val*NM)));
    elseif min(dif)>MIN
        valeq=E1.^(-2).*(E1+(-1).*E12).^(-1).*E12.^(-2).*(E1+(-1).*E3).^(-1).*( ...
          E12+(-1).*E3).^(-1).*E3.^(-2).*(((-1)+exp(1).^(E1.*NM)).*E12.^2.*( ...
          E12+(-1).*E3).*E3.^2+E1.*E12.^2.*E3.^2.*((-1).*E12+E3).*NM+E1.^3.* ...
          ((-1).*((-1)+exp(1).^(E12.*NM)).*E3.^2+E12.*E3.^2.*NM+E12.^2.*(( ...
          -1)+exp(1).^(E3.*NM)+(-1).*E3.*NM))+E1.^2.*(((-1)+exp(1).^(E12.* ...
          NM)).*E3.^3+(-1).*E12.*E3.^3.*NM+E12.^3.*(1+(-1).*exp(1).^(E3.*NM) ...
          +E3.*NM)));
    else
        [~,index]=min(dif);
        val=vec(index);
        others=vec; others(index)=[];
        oval=0.5*sum(others);
        valeq=((-3*oval*expl(4,NM*oval)+NM*oval^2*expl(3,NM*oval))*val^2+...
          (2*expl(4,NM*oval)-NM*oval*expl(3,NM*oval))*val^3+...
          oval^3*expl(4,NM*val))/(oval^3*val^2*(oval-val)^2);       
        
    end
    
end
end