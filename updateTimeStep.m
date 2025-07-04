function [neurons,reset,mm,mmi]= updateTimeStep(neurons,omega,timeParams,Idbs)
global indIpost
% calculate the uv time dependent variables
sN = size(neurons.v);
uTemp = neurons.u;
vTemp = neurons.v;
neurons.u = neurons.u +  timeParams.dt*(neurons.a.*(neurons.b.*neurons.v - neurons.u));

n = size(omega,1)*size(omega,2)*size(omega,3)*size(omega,4);
%a = reshape(omega,1,size(omega,1)*size(omega,2)*size(omega,3)*size(omega,4));
%b = reshape(neurons.IpostSynaptic,1,size(neurons.IpostSynaptic,1)*size(neurons.IpostSynaptic,2)*size(neurons.IpostSynaptic,3));

%--Replace line 16 with 15 if you want to run with Mex on a windows machine
% A = multMex(omega,neurons.IpostSynaptic,n,indIpost);
A = reshape(omega,1,size(omega,1)*size(omega,2)*size(omega,3)*size(omega,4)).*neurons.IpostSynaptic(indIpost);
%A = times(reshape(omega,1,size(omega,1)*size(omega,2)*size(omega,3)*size(omega,4)),neurons.IpostSynaptic(indIpost));%%fixMe for delays applied here
%A = parMult(reshape(omega,1,size(omega,1)*size(omega,2)*size(omega,3)*size(omega,4)),neurons.IpostSynaptic(indIpost));
A = reshape(A,size(omega));
Asum = sum(A,1);
Asum1 = sum(Asum,2);
Asum1 = reshape(Asum1,size(Asum1,3),size(Asum1,4));

mm = max(Asum1,[],'all');
mmi = min(Asum1,[],'all');
%A = poissrnd(timeParams.lambda,sN);
poissonSpikes = poissrnd(timeParams.lambda,sN).*neurons.mu;%poission spike events must be multiplied by mu
A = normrnd(neurons.WhiteNoiseMean,neurons.WhiteNoiseSigma,size(neurons.v));
whiteNoise = neurons.WhiteNoiseIntensity*A;%?(t)



neurons.v = neurons.v + timeParams.dt*(0.04.*neurons.v.*neurons.v + 5*neurons.v ...
    - uTemp + 140 + neurons.inputCurrent ...
    + Asum1 ) + whiteNoise +  ...
    Idbs*neurons.IdbsCondition + poissonSpikes;

if(0)
    if(abs(sum(sum(Asum1)))>0)
        s1= sum(sum(neurons.v));
        s2 = sum(sum(0.04.*neurons.v.*neurons.v))*timeParams.dt;
        s3 = sum(sum(5*neurons.v))*timeParams.dt;
        s4 = sum(sum(uTemp))*timeParams.dt;
        s5 = sum(sum(neurons.inputCurrent))*timeParams.dt;
        s6 = sum(sum(Asum1))*timeParams.dt;
        mm=s6;
        s7 = sum(sum(whiteNoise));
        s8 = sum(sum(poissonSpikes));
        s9 = sum(sum(Idbs*neurons.IdbsCondition));
        
        stop=1;
    end
end
Pnoise = normrnd(neurons.PthreshMean,neurons.PthreshSigma,size(neurons.v));%?(t)

reset = vTemp >= neurons.vp + Pnoise;%resetting condition (spike events)
neurons.u(reset) = neurons.u(reset) + neurons.d(reset);
neurons.v(reset) = neurons.c(reset);

%in case of spike events calculate the post synaptic current of all neurons
indSpike(:,:) = reset;
neurons = computePostSynapticCurrent_TM(neurons,timeParams,indSpike);
%neurons.r = neurons.r + timeParams.dt*(-neurons.r./neurons.tauf) + neurons.U.*(1-neurons.r).*indSpike;
%neurons.x = neurons.x + timeParams.dt*((1-neurons.x)./neurons.taud) - neurons.r.*neurons.x.*indSpike;
%neurons.I = neurons.I + timeParams.dt*(-neurons.I./neurons.taus) + neurons.A.*neurons.r.*neurons.x.*indSpike;

%neurons.I(isnan(neurons.I))=0;

% aa = mean(poissonSpikes,'all')/20;
%pushing the last post synaptic current into the buffer
IpostTemp = neurons.IpostSynaptic;
%type = neurons.type;
%isExcitory   = type==1|type==2|type==5;
%isInhibitory = type==3|type==4|type==6;
%isInactive   = (~isInhibitory)&(~isExcitory);
%neurons.IpostSynaptic(:,:,1) = -(isInhibitory).*neurons.I + (isExcitory).*neurons.I;
neurons.IpostSynaptic(:,:,1)  = neurons.I;
neurons.IpostSynaptic(:,:,2:end) = IpostTemp(:,:,1:end-1);
% reset inactive neurons
 inactive = neurons.type==0;
neurons.u(inactive) = 0;
neurons.v(inactive) = 0;
neurons.I(inactive) = 0;


% reset = reset(:,5:6);
