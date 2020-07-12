function [ Tn Tm Th ] = tau(V,T)
% input: V in mV    T in deg C

global Vr

[ An Am Ah ] = alpha(V,T);

[ Bn Bm Bh ] = beta(V,T);


Tn = 1./(An + Bn);

Tm = 1./ (Am + Bm);

Th = 1./ (Ah + Bh);

end



