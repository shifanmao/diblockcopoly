function out = expl(m,x)
%This is a function used to acuratly evalute the tylor series expansion
% for exp(x) with the first m terms removed.
%Created 6/22/15
%I haven't tested this function for complex x.
    fact=[1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800];
    if floor(m)~=m
        error('First arguement must be an interger');
    end
    if isscalar(x)
        if abs(x) < 0.1
            if m==0
                out=exp(x);
                return
            end
            out=0;
            for j = m:11
                out=out+x^j/fact(j+1); %the j+1 is becuase matlab starts at 1
            end
        else
            out=exp(x);
            if m==0
                return
            end
            for j=0:m-1
                out=out-x^j/fact(j+1); %the j+1 is becuase matlab starts at 1
            end
        end
    else 
        out=zeros(size(x));
        for j = m:11
            out=out+(abs(x)<0.1).*x.^j/factorial(j);
        end
        out=out+(abs(x)>=0.1).*exp(x);
        if m==0
            return
        end
        for j=0:m-1
            out=out-(abs(x)>=0.1).*x.^j/factorial(j);
        end

    end
end
