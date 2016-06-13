function s2=s2rr(N,FA,K)

s2=zeros(2,2);

KJ=abs(K);
for I1=1:N
    for I2=1:N
      sep = abs(I1-I2);
      pgc = zeros(2,2);
      if sep==0
          pgc(1,1)=1*pa1a2(1,1,FA,N,I1,I2);
          pgc(1,2)=1*pa1a2(1,2,FA,N,I1,I2);
          pgc(2,1)=1*pa1a2(2,1,FA,N,I1,I2);
          pgc(2,2)=1*pa1a2(2,2,FA,N,I1,I2);
      else
          pgc(1,1)=sin(KJ*sep)/(KJ*sep)*pa1a2(1,1,FA,N,I1,I2);
          pgc(1,2)=sin(KJ*sep)/(KJ*sep)*pa1a2(1,2,FA,N,I1,I2);
          pgc(2,1)=sin(KJ*sep)/(KJ*sep)*pa1a2(2,1,FA,N,I1,I2);
          pgc(2,2)=sin(KJ*sep)/(KJ*sep)*pa1a2(2,2,FA,N,I1,I2);
      end
      s2(1,1)=s2(1,1)+pgc(1,1);
      s2(1,2)=s2(1,2)+pgc(1,2);
      s2(2,1)=s2(2,1)+pgc(2,1);
      s2(2,2)=s2(2,2)+pgc(2,2);
    end
end

end
 
function [PA1A2]=pa1a2(A1,A2,FA,N,I1,I2)
    IND1=(I1-1)/N;
    IND2=(I2-1)/N;
    PA1A2=0;
    
    if ( ((IND1<FA)&&(A1==1) || (IND1>=FA)&&(A1==2)) &&...
         ((IND2<FA)&&(A2==1) || (IND2>=FA)&&(A2==2)) )
        PA1A2=1;
    end
end
