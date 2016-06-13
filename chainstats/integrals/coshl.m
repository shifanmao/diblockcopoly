function out=coshl(m,x)
% Calculates cosh less the first m talor series values.
% Note that cosh is even so coshl(1,x)=coshl(2,x)
    if floor(m)~=m
        error('First arguement must be an interger');
    end
    out=0.5*(expl(m,x)+expl(m,-x));
end