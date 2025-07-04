% for the time being we only consider inter-layer delays
% further expansion should consider all to all delays between all pairs of
% neurons

function delay = initDelays(nN,nL)
delay = 10*ones(nL,nL);% for the time being we only consider inter-layer delays
for i=1:nL
    for j=1:nL
        if (i==j)
            delay(j,i)=1;
        else
             delay(j,i)=8;
        end
    end
end
delay = delay/1000;
% ind = 1:nL*nL;
% ind2sub
