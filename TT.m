function re=TT(h1,h2,parameters)
%transpose(T)
re=(h1==0).*(h2==0)*parameters.T0'+...
    (h1==-1).*(h2==0)*parameters.T1'+...
    (h1==0).*(h2==-1)*parameters.Tm1';
end