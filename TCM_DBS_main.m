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
% dbsType = 'cDBS'; %conventional DBS, i.e., a continoius pulse train with a constant frequency 
dbsType = 'DBS A'; %a pulse train with an initial phase of 20 pulses at 130 Hz and then switch to ~95 Hz.
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
%Later we should merge output fig and analysis together as
%AnalysisWoutputFig

% Output_Fig
% thm_fig

% Network_rasterp00lot;


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

%% ANALYSIS
% TCM_DBS_Analysis

%% Morgera's synchrony index analysis
if fx~=0
    pdFinit = dbsOnset*1000/timeParams.dt;
    pdInit = pdFinit - 5*1000/timeParams.dt; %select the last 5 second of the PD state
    dbsInit = dbsOnset*1000/timeParams.dt+1*1000/timeParams.dt; %select a portion during DBS 1 second after the onset
    dbsFinit = dbsInit + 5*1000/timeParams.dt; %for 1 second
    
    %choose all the cortical pyramidal cells for the estimation
    v_cortex_pd = [vRS1(1:2,pdInit:pdFinit);vIB1(1:2,pdInit:pdFinit);...
                    vRS2(1:2,pdInit:pdFinit);...
                    vRS3(1:2,pdInit:pdFinit);vIB3(1:2,pdInit:pdFinit);...
                    % vFS4(:,pdInit:pdFinit); vLTS4(:,pdInit:pdFinit)...
                    ];
    v_cortex_dbs = [vRS1(1:2,dbsInit:dbsFinit);vIB1(1:2,dbsInit:dbsFinit);...
                    vRS2(1:2,dbsInit:dbsFinit);...
                    vRS3(1:2,dbsInit:dbsFinit);vIB3(1:2,dbsInit:dbsFinit);...
                    % vFS4(:,dbsInit:dbsFinit); vLTS4(:,dbsInit:dbsFinit)
                    ];

    M_pd = Morgera_index(v_cortex_pd)
    M_dbs = Morgera_index(v_cortex_dbs)
else
    pdFinit = simTime*1000/timeParams.dt;    %choose the whole simulation
    pdInit = pdFinit - (simTime-1)*1000/timeParams.dt; %exclude the first second of the simulation

    %choose all the cortical pyramidal cells for the estimation
    v_cortex_pd = [vRS1(:,pdInit:pdFinit);vRS2(:,pdInit:pdFinit);vRS3(:,pdInit:pdFinit);vIB1(:,pdInit:pdFinit);vIB3(:,pdInit:pdFinit);...
                    %vFS4(:,pdInit:pdFinit); vLTS4(:,pdInit:pdFinit)];
                    ];
    M = Morgera_index(v_cortex_pd)
end

%% Kuramoto's order parameter
% v_test = [vRS1;vIB1;vRS2];
% rho = Kuramoto(v_test, timeParams.dt);
% figure; plot(rho)

%% PSDs
fs = 1000/timeParams.dt; frange = 1:200;

if fx~=0
    [pECoG1,fECoG1] = pwelch(ECoG(pdInit:pdFinit),fs,fs/2,frange,fs); %psd of the PD state ECoG
    [pECoG2,fECoG2] = pwelch(ECoG(dbsInit:dbsFinit),fs,fs/2,frange,fs); %psd of the PD+DBS state ECoG
%     pECoG2=pECoG2./sum(pECoG2); %normalize
%     pECoG1=pECoG1./sum(pECoG1);
    
    BetaPower_pd = trapz(pECoG1(13:30));
    BetaPower_pd_dbs = trapz(pECoG2(13:30));
    DeltaP = BetaPower_pd - BetaPower_pd_dbs
    
    figure;
    plot(fECoG1, pECoG1, 'LineWidth', 2); hold on
    plot(fECoG2, pECoG2, 'LineWidth', 2); zoom xon
%     xlim([5 200])
    xlabel('Frequency (Hz)')
    ylabel('PSD (a.u.)')
    legend('DBS off',['DBS on (',num2str(fx),' Hz)'])
else
    [pECoG1,fECoG1] = pwelch(ECoG(pdInit:pdFinit),fs,fs/2,frange,fs); %psd of the PD state ECoG
%     pECoG1=pECoG1./sum(pECoG1); %normalize

    BetaPower_pd = trapz(pECoG1(13:30));

%     figure;
%     plot(fECoG1, pECoG1, 'LineWidth', 2); hold on
% %     xlim([5 200])
%     xlabel('Frequency (Hz)')
%     ylabel('PSD (a.u.)')
%     legend('DBS off')
end
% set(gca,'FontSize',12,'FontWeight','bold')

beep