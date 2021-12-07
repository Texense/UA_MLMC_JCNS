function [DSE,DSI] = ISI(spike, NE,NI, dt)
DSE = []; DSI = [];
for ii = 1:NE
    Sp_tii = find(spike(ii,:) == 1)*dt;
    dif_ii = Sp_tii(2:length(Sp_tii)) - Sp_tii(1:length(Sp_tii)-1);
    DSE = [DSE,dif_ii];
    clear Sp_tii dif_ii
end

for jj = (NE+1):(NE+NI)
    Sp_tjj = find(spike(jj,:) == 1)*dt;
    dif_jj = Sp_tjj(2:length(Sp_tjj)) - Sp_tjj(1:length(Sp_tjj)-1);
    DSI = [DSI,dif_jj];
    clear Sp_tjj dif_jj
end

end