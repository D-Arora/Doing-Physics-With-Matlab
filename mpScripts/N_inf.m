function [ n_inf m_inf h_inf ] = N_inf(V,T)
% input: membrane potential in mV    T in deg C
global Vr

[ An Am Ah ] = alpha(V,T);

[ Bn Bm Bh ] = beta(V,T);


n_inf = An ./(An + Bn);

m_inf = Am ./(Am + Bm);

h_inf = Ah ./ (Ah + Bh);

end



