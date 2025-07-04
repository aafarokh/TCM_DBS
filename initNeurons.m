function neurons =  initNeurons(nN,nL,maxDelay, dbsLayer, dbsFid)
%initial values of u and v
% neurons.u = zeros(nN,nL);

%v=-65*ones(Ne+Ni,1); % Initial values of v
%u=b.*v; % Initial values

% rng(27364283) %control the random numbers generation


neurons.vp   = ones(nN,nL)*30; %spike potential threshold
% type of neurons from 1 to 6 represents RS,IB,FS,LTS,TC,TR table S2 in supplementary material
neurons.type = zeros(nN,nL,'int32');% zero stands for inactive neurons
% S layer
neurons.type(1:end/2,1)         = 1;%1;%RS 50 percent
neurons.type(end/2 + 1:end,1)   = 2;%2;%IB 50 percent
% M layer
neurons.type(:,2)               = 1;%1;%RS 100 percent
% D layer
neurons.type(1:.7*end,3)        = 1;%1;%RS 70 percent
neurons.type(.7*end + 1:end,3)  = 2;%2;%IS 30 percent
% INs layer
neurons.type(1    :end/2,4)     = 3;%3;%FS 50 percent
neurons.type(end/2 + 1:end,4)   = 4;%4;%LTS 50 percent
% rel layer
neurons.type(:,5) = 5;%TC 100 percent
% ret layer
neurons.type(:,6) = 6;%TR 100 percent

% inactive neurons
neurons.type(41:end,6) = 0; %layer six only contains 40 active neurons

type = neurons.type;

%% use types to initilize intrinsic neuron parameters
% a values
neurons.a = zeros(nN,nL);
rr = .0*rand(nN,nL);%% amirali should set this parameter
neurons.a(type==1) = 0.02   + rr(type==1);
neurons.a(type==2) = 0.02   + rr(type==2);
neurons.a(type==3) = 0.1    + rr(type==3);
neurons.a(type==4) = 0.02   + rr(type==4);
neurons.a(type==5) = 0.02   + rr(type==5);
neurons.a(type==6) = 0.02   + rr(type==6);
% b values
neurons.b = zeros(nN,nL);%we should add randomness later
neurons.b(type==1) = 0.2;
neurons.b(type==2) = 0.2;
neurons.b(type==3) = 0.2;
neurons.b(type==4) = 0.25;
neurons.b(type==5) = 0.25;
neurons.b(type==6) = 0.25;
% c values
neurons.c = zeros(nN,nL);%we should add randomness later
neurons.c(type==1) = -65;
neurons.c(type==2) = -55;
neurons.c(type==3) = -65;
neurons.c(type==4) = -65;
neurons.c(type==5) = -65;
neurons.c(type==6) = -65;
%initial values of v
neurons.v = neurons.c;
%initial values of u
neurons.u = neurons.b.*neurons.v;
% d values
neurons.d = zeros(nN,nL);
neurons.d(type==1) = 8;
neurons.d(type==2) = 4;
neurons.d(type==3) = 2;
neurons.d(type==4) = 2;
neurons.d(type==5) = 0.05;
neurons.d(type==6) = 2.05;
% input current values
neurons.inputCurrent = zeros(nN,nL);%Iij is just above the firing threshold
neurons.inputCurrent(type==1) = 2.5;
neurons.inputCurrent(type==2) = 2.5;
neurons.inputCurrent(type==3) = 3.2;
neurons.inputCurrent(type==4) = 0.0;
neurons.inputCurrent(type==5) = 1.0-1.0;
neurons.inputCurrent(type==6) = 1.0-0.5;
% r,x,I initial values
neurons.r = zeros(nN,nL);
neurons.x = ones(nN,nL);%fix for correct initial value
%neurons.x = .5*ones(nN,nL);
neurons.I = zeros(nN,nL);% post Synaptic Currents

% synaptic coefficients
neurons.IpostSynaptic = zeros(nN,nL,maxDelay);% post Synaptic Currents buffer.third size is set to max delay for buffering previous values
neurons.taufExcitory   = [670;17;326];
neurons.taufInhibitory = [376;21;62];
neurons.taudExcitory   = [138;671;329];
neurons.taudInhibitory = [45;706;144];
neurons.UExcitory      = [0.09;0.5;0.29];
neurons.UInhibitory    = [0.016;0.25;0.32]; 
% neurons.taus = 2.5;%should add different values for inhibitory and excitatory neurons separately
neurons.A = 1.0;

neurons.FDPdistExcitory   = [0.2,0.63,0.17];%fdp distributions for excitatory
neurons.FDPdistInhibitory = [.08,0.75,0.17];%fdp distributions for inhibitory


isExcitory   = type==1|type==2|type==5;
isInhibitory = type==3|type==4|type==6;
rand1 = rand(size(neurons.v));
rand2 = rand(size(neurons.v));
rand3 = rand(size(neurons.v));
rand4 = rand(size(neurons.v));

neurons.a(isInhibitory) = neurons.a(isInhibitory) .*(1.0  +  rand1(isInhibitory)*0.08/0.02);
neurons.b(isInhibitory) = neurons.b(isInhibitory) .*(1.0  -  rand2(isInhibitory)*0.05/.25);
neurons.c(isExcitory)   = neurons.c(isExcitory)   .*(1.0  +  rand3(isExcitory).^2*15/-65);
neurons.d(isExcitory)   = neurons.d(isExcitory)   .*(1.0  +  rand4(isExcitory).^2*(-6/8));


%% tauf and taud and U
neurons.fdpType = zeros(nN,nL,'int32');
randDist = rand(nN,nL);

fepercent = 0.2; depercent = 0.63;
fipercent = 0.08; dipercent = 0.75;

neurons.fdpType(isExcitory)   = (randDist(isExcitory)<=fepercent)*1 + ...
                                (randDist(isExcitory)>fepercent&randDist(isExcitory)<=depercent +fepercent)*2 + ...
                                (randDist(isExcitory)>depercent + fepercent  &randDist(isExcitory)<=1.0)*3;
neurons.fdpType(isInhibitory) = (randDist(isInhibitory)<=fipercent)*1 + ...
                                (randDist(isInhibitory)>fipercent&randDist(isInhibitory)<=dipercent +fipercent)*2 + ...
                                (randDist(isInhibitory)>dipercent + fipercent  &randDist(isInhibitory)<=1.0)*3;
                          
neurons.tauf = zeros(nN,nL);
neurons.tauf(isExcitory)   = neurons.taufExcitory(neurons.fdpType(isExcitory)); 
neurons.tauf(isInhibitory) = neurons.taufInhibitory(neurons.fdpType(isInhibitory)); 

neurons.taud = zeros(nN,nL);
neurons.taud(isExcitory)   = neurons.taudExcitory(neurons.fdpType(isExcitory)); 
neurons.taud(isInhibitory) = neurons.taudInhibitory(neurons.fdpType(isInhibitory)); 

neurons.U = zeros(nN,nL);
neurons.U(isExcitory)   = neurons.UExcitory(neurons.fdpType(isExcitory)); 
neurons.U(isInhibitory) = neurons.UInhibitory(neurons.fdpType(isInhibitory)); 

neurons.taus = zeros(nN,nL);
neurons.taus(isExcitory) = 2;
neurons.taus(isInhibitory) = 8;
%% noise parameters
neurons.mu = 0.0*ones(nN,nL);% amplitude of poission random spikes


neurons.PthreshMean = .0;%ζ(t) is a Gaussian white noise added to individual neurons threshold potentials.
                          %The mean of this threshold noise was set to half
                          %of the additive white gaussian noise.


neurons.PthreshSigma = .1;%ζ(t) is a Gaussian white noise added to individual neurons threshold potentials.
                            %The mean of this threshold noise was set to half
                            %of the additive white gaussian noise.

neurons.WhiteNoiseMean = .0;  %white gaussian noise, ξ(t), to the membrane current 
                                %of each individual neuron adjusts the mean
                                %firing rates compatible with experimental recordings
                                
neurons.WhiteNoiseSigma = .5; %white gaussian noise, ξ(t), to the membrane current 
                                %of each individual neuron adjusts the mean
                                %firing rates compatible with experimental recordings
neurons.WhiteNoiseIntensity = 1.0; %white gaussian noise, ξ(t), to the membrane current 
                                  %of each individual neuron adjusts the mean
                                  %firing rates compatible with experimental recordings

neurons.IdbsCondition = zeros(size(neurons.v),'logical');
neurons.IdbsCondition(1:dbsFid,dbsLayer) = true;%% layer (:,j) is triggered with (i%,:) fidelity
%  neurons.IdbsCondition(1:50,[3,4]) = true;%% together with i% of CIs
