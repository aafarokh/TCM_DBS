%% parameters initialization

clearvars; 
rng(22823303,'twister') %control the random numbers generation

% close all;
%mex CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'multMex.c'
nL = 6;%number of layers
nN = 100;%number of neurons per layer
n_tot = 540; %total number of the neurons in the network
omega = initOmega(nN,nL);% inter-neuron post synaptic current correlation
%omega = initOmega_thmOnly(nN,nL);% inter-neuron post synaptic current correlation

simTime = 12; % simulation time in seconds
dbsOnset = 6; % DBS onset in seconds
fx = 130;     % fx is the dbs frequency, set it to zero for no stimulation
Ax = 125;      % Ax is dbs amplitude         
dbsLayer = [3]; %the layer(s) that receives DBS directly (add to the array the desired layer number)
dbsFid = 50;   %number of neurons in a layer to receive direct DBS injection 
dbsType = 'cDBS'; %conventional DBS, i.e., a continoius pulse train with a constant frequency 
% dbsType = 'DBS A'; %a pulse train with an initial phase of 20 pulses at 130 Hz and then switch to ~95 Hz.
% dbsType = 'DBS B'; %a pulse train with an initial phase of 20 pulses at 130 Hz and then pause for 37 ms then start 4 pulses at 130 Hz and so forth.

timeParams = initTimeParams(nN,nL,Ax,fx,dbsOnset,simTime,dbsType);
neurons = initNeurons(nN,nL,ceil(timeParams.maxDelay), dbsLayer, dbsFid);% neurons structure including all neurons time-dependent and time-independent variables

%% start the simulation and store the data

% initilize delay indices
[iPrime,jPrime,i,j]= ind2sub(size(omega),1:size(omega,1)*size(omega,2)*size(omega,3)*size(omega,4));

% set omega link to zero for inactive neurons
inactive     = neurons.type == 0;
inactiveIND1 = sub2ind(size(inactive),iPrime,jPrime);
inactiveIND2 = sub2ind(size(inactive),i,j);
in1 = inactive(inactiveIND1);
in2 = inactive(inactiveIND2);
omega(in1)=0;
omega(in2)=0;

indDelay = sub2ind(size(timeParams.delay),jPrime,j);
dd = timeParams.delay(indDelay);
global indIpost
indIpost = int32(sub2ind(size(neurons.IpostSynaptic),iPrime,jPrime,dd));%post synaptic buffer index
clear iPrime jPrime j indDelay dd in1 in2 inactiveIND1 inactiveIND2

%% start simulation
rasterData = zeros(length(timeParams.tVec),nN,nL,'logical');%for storing raster plot
tic
for i=1:length(timeParams.tVec)

    if(mod(i,1000)==0)
        clc
        if fx~=0
    disp(['Stimulation at ',num2str(fx),' Hz'])
        else 
    disp('No stimulation')
        end
    disp('progress percentage = ');
    disp(round(i/length(timeParams.tVec)*100));
    end
    
    
    [neurons1,rasterData(i,:,:),mm(i),mmi(i)] = updateTimeStep(neurons,omega,timeParams,timeParams.dbs(i));
    % save APs for layer 1
    vRS1(:,i) = neurons.v(1:50,1);
    vIB1(:,i) = neurons.v(51:100,1);
    
    % save APs for layer 2
    vRS2(:,i) = neurons.v(:,2);
 
    % save APs for layer 3
    vRS3(:,i) = neurons.v(1:70,3);
    vIB3(:,i) = neurons.v(71:100,3);
    
    % save APs for layer 4
    vFS4(:,i) = neurons.v(1:50,4);
    vLTS4(:,i) = neurons.v(51:100,4);
    
    % save APs for layer 5
    vTC5(:,i) = neurons.v(:,5);

    % save APs for layer 6
    vTR6(:,i) = neurons.v(:,6);
    
    % plot post synaptic current of layer six
    I1(i) = sum(neurons.IpostSynaptic(:,1,1));
    I2(i) = sum(neurons.IpostSynaptic(:,2,1));
    I3(i) = sum(neurons.IpostSynaptic(:,3,1));
    I4(i) = sum(neurons.IpostSynaptic(:,4,1));
    I5(i) = sum(neurons.IpostSynaptic(:,5,1));
    I6(i) = sum(neurons.IpostSynaptic(:,6,1));
    I14(i) = sum(neurons.IpostSynaptic(:,1,1))+sum(neurons.IpostSynaptic(:,4,1));
    I24(i) = sum(neurons.IpostSynaptic(:,2,1))+sum(neurons.IpostSynaptic(:,4,1));
    I34(i) = sum(neurons.IpostSynaptic(:,3,1))+sum(neurons.IpostSynaptic(:,4,1));
    I56(i) = sum(neurons.IpostSynaptic(:,5,1))+sum(neurons.IpostSynaptic(:,6,1));
    ECoG(i) = 0.6 * I1(i) + 0.3 * I2(i) + 0.1 * I3(i) - I4(i); %deifne the ECoG signal as a weighted sum of cortical postsynaptic currents
    neurons = neurons1;

end
toc;

%% Output Figures

% Output_Fig
% thm_fig

% Network_rasterplot;

%%
% figure(10)
% title('Rt example APs'); hold on
% plot(timeParams.tVec,vTR6(1,:),'k'); hold on
% plot(timeParams.tVec,vTR6(10,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
% plot(timeParams.tVec,vTR6(20,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
% plot(timeParams.tVec,vTR6(40,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
% zoom xon
% xlabel('Time (ms)')
% save('TCM_DBS_Sim_Results')

beep
