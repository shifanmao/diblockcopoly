function [val]=legendrep(P,ORDL)

val=zeros(length(P),ORDL);

val(:,1)=ones(length(P),1);
if ORDL>=2
   val(:,2)=P;
end
   
for N=3:ORDL
    L=N-2;
    val(:,L+1+1)=((2*L+1)*P.*val(:,L+1)-L*val(:,L-1+1))/(L+1);
end

val=transpose(val);
end