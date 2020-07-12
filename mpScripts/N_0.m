function [ n_0 m_0 h_0 ] = N_0(T);
% V in mV and T in degC

global Vr

V_0 = Vr;

[ An Am Ah ] = alpha(V_0,T);

[ Bn Bm Bh ] = beta(V_0,T);


n_0 = An ./(An + Bn);

m_0 = Am ./ (Am + Bm);

h_0 = Ah ./ (Ah + Bh);

end



