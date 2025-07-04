function timeParams = initTimeParams(nN,nL,Ax,fx,dbsOnset,simTime,dbsType)
timeParams.dt =.1;%time steps in milli-seconds
timeParams.tStart = 0.0*1000;%start time in milli-secs
timeParams.tEnd = simTime*1000;%end time in milli-secs
timeParams.fdbs = fx/1000;%dbs frequency in hertz/milli-second
timeParams.dbsAmplitude = Ax;% parameter A in dbs definition(dbs amplitude)
timeParams.dbsStartTime = dbsOnset*1000; %dbs starting point
timeParams.fr = 20/1000;%poisson noise mean frequency in hertz/milli-second
timeParams.tVec = timeParams.tStart:timeParams.dt:timeParams.tEnd;%time vector
timeParams.lambda = timeParams.fr*timeParams.dt;% for poission distribution
% make dbs pulses
if strcmp(dbsType,'cDBS')
timeParams.dbs = genPatternCDBS(timeParams);
end
if strcmp(dbsType,'DBS A')
timeParams.dbs = genPatternA(timeParams,1000/130*20,130/1000,95/1000);
end
if strcmp(dbsType,'DBS B')
timeParams.dbs = genPatternB(timeParams,1000/130*20,130/1000,130/1000,37);
end


timeParams.dbs(length(timeParams.tVec)) = 0;% make last one zero (nothing special)
% multiply dbs pulses by amplitude

timeParams.dbs = timeParams.dbs.*timeParams.dbsAmplitude;
%timeParams.dbs = genPatternA(timeParams,1000/130*20,130/1000,95/1000);
%timeParams.dbs = genPatternB(timeParams,1000/130*20,130/1000,130/1000,37);
delay = zeros(nL,nL);% for the time being we only consider inter-layer delays
% make delays for intra-layer and inter-layer delays
for i=1:nL
    for j=1:nL
        if (i==j)
            delay(j,i) = 1;%intra-layer delay in milliseconds
        else
            delay(j,i) = 2;%inter-layer delay in milliseconds
        end
    end
end
timeParams.delay = round(delay/timeParams.dt) + 1;%convert to time step units add one for converting to index;
%these values are
%rounded for
%simplicity and can
%be interpolated for
%greater accuracy


timeParams.maxDelay = max(max(timeParams.delay));% retreive maximum delay step to change
% the size of delayed post synaptic current buffer


figure(1000)
plot(timeParams.dbs);

%timeParams.poissonSpikes = poissrnd(lambda,size(timeParams.tVec));%poission spike events
end
function dbs = genPatternCDBS(timeParams)

for i=1:length(timeParams.tVec)-1
    period = 1.0/timeParams.fdbs;
    t = timeParams.tVec(i);
    tNext = timeParams.tVec(i+1);
    if (t>timeParams.dbsStartTime)
        dbs(i) = floor(tNext/period)- floor(t/period);
    end
end
end

function dbs = genPatternA(timeParams,Tphase1,freq1,freq2)
for i=1:length(timeParams.tVec)-1
    t = timeParams.tVec(i);
    if(t<Tphase1+timeParams.dbsStartTime)
        period = 1.0/freq1;
    else
        period = 1.0/freq2;
    end

    tNext = timeParams.tVec(i+1);
    if (t>timeParams.dbsStartTime)
        dbs(i) = floor(tNext/period)- floor(t/period);
    end
end

end

function dbs = genPatternB(timeParams,Tphase1,freq1,freq2,Tpause)
counter = 0;
pause = 0;
pauseStart = 0;
for i=1:length(timeParams.tVec)-1
    t = timeParams.tVec(i);
    if(t<Tphase1+timeParams.dbsStartTime)
        period = 1.0/freq1;
        phase = 0;
    else
        period = 1.0/freq2;
        phase = 1;
    end

    tNext = timeParams.tVec(i+1);

    if (t>timeParams.dbsStartTime)

        dbs(i) = floor(tNext/period)- floor(t/period);
        if(phase==1)
            if ((dbs(i)>0)&&(~pause))
                counter = counter+1;
            end
            if(counter==4)
                pauseStart = t;
                pause=1;
                counter=0;
                continue;
            end
            if(pause)
                dbs(i)=0;
            end

            if(t-pauseStart>Tpause)&&pause
                pause = 0;
            end

        end
    end
end

end
