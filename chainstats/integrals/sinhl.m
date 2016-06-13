function out=sinhl(m,x)
% Calculates cosh less the first m talor series values.
    if floor(m)~=m
        error('First arguement must be an interger');
    end
    out=0.5*(expl(m,x)-expl(m,-x));
end