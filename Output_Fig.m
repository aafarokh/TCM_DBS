%>this is the output figure .m file for TCM_DBS

% All rasters column 1
spikes_S1=zeros(nN/2,50000); spikes_M=zeros(nN,50000);
spikes_S2=zeros(nN/2,50000);
spikes_D1=zeros(.7*nN,50000); spikes_IC1=zeros(nN/2,50000);
spikes_D2=zeros(.3*nN,50000); spikes_IC2=zeros(nN/2,50000);
spikes_TC=zeros(nN,50000); spikes_TR=zeros((nN/2)-(nN/10),50000);

fh = figure;
ax1 = subplot(6,5,[1,6,11,16,21,26]); title('Network Raster Plot'); hold on

nRT = (nN/2)-(nN/10);
for j=1:nRT
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vTR6(j,i)>neurons.vp(1,6)
            spikes_TR(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[j-1 j],'LineWidth',1,'Color','c')
            %     drawnow
        end
    end
end

for j=1:nN
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vTC5(j,i)>neurons.vp(1,5)
            spikes_TC(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[nRT+j-1 j+nRT],'LineWidth',1,'Color','m')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vFS4(j,i)>neurons.vp(1,4)
            spikes_IC1(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[nN+nRT+j-1 j+nN/2+nRT],'LineWidth',1,'Color','k')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vLTS4(j,i)>neurons.vp(1,4)
            spikes_IC2(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[nN+nN/2+nRT+j-1 j+nN+nRT],'LineWidth',1,'Color','g')
        end
    end
end

for j=1:30
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vIB3(j,i)>neurons.vp(1,3)
            spikes_D2(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[2*nN+nRT+j-1 j+nN+nN/2+nRT],'LineWidth',1,'Color','y')
        end
    end
end

for j=1:70
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vRS3(j,i)>neurons.vp(1,3)
            spikes_D1(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[2*nN+nN/2+nRT+j-1 j+2*nN+nRT],'LineWidth',1,'Color','y')
        end
    end
end

for j=1:nN
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vRS2(j,i)>neurons.vp(1,2)
            spikes_M(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[3*nN+nRT+j-1 j+3*nN+nRT],'LineWidth',1,'Color','y')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vIB1(j,i)>neurons.vp(1,1)
            spikes_S2(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[4*nN+nRT+j-1 j+3*nN+nN/2+nRT],'LineWidth',1,'Color','y')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vRS1(j,i)>neurons.vp(1,1)
            spikes_S1(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[4*nN+nN/2+nRT+j-1 j+4*nN+nRT],'LineWidth',1,'Color','y')
        end
    end
end

ylabel('Neuron #');

%% PSTHs column 2
ax2 = subplot(6,5,2); title('PSTH of S'); hold on
binlength = 100;
for bin=2:500
    psth_TC(bin-1) = sum(spikes_TC(1:100,(bin-1)*binlength:bin*binlength),'all');
end
% plot(psth)
bar(psth_TC)
% xlim([300 500])
% xlabel('Time (sec)')
ylabel('# of spikes')

ax3 = subplot(6,5,7); title('PSTH of M'); hold on
binlength = 100;
for bin=2:500
    psth_TC(bin-1) = sum(spikes_TC(1:100,(bin-1)*binlength:bin*binlength),'all');
end
% plot(psth)
bar(psth_TC)
% xlim([300 500])
% xlabel('Time (sec)')
ylabel('# of spikes')
% % 

ax4 = subplot(6,5,12); title('PSTH of D'); hold on
binlength = 100;
for bin=2:500
    psth_TC(bin-1) = sum(spikes_TC(1:100,(bin-1)*binlength:bin*binlength),'all');
end
% plot(psth)
bar(psth_TC)
% xlim([300 500])
% xlabel('Time (sec)')
ylabel('# of spikes')

ax5 = subplot(6,5,17); title('PSTH of CI'); hold on
binlength = 100;
for bin=2:500
    psth_TC(bin-1) = sum(spikes_TC(1:100,(bin-1)*binlength:bin*binlength),'all');
end
% plot(psth)
bar(psth_TC)
% xlim([300 500])
% xlabel('Time (sec)')
ylabel('# of spikes')

ax6 = subplot(6,5,22); title('PSTH of TC'); hold on
binlength = 100;
for bin=2:500
    psth_TC(bin-1) = sum(spikes_TC(1:100,(bin-1)*binlength:bin*binlength),'all');
end
% plot(psth)
bar(psth_TC)
% xlim([300 500])
% xlabel('Time (sec)')
ylabel('# of spikes')

ax7 = subplot(6,5,27); title('PSTH of Rt'); hold on
for bin=2:500
    psth_TR(bin-1) = sum(spikes_TR(:,(bin-1)*binlength:bin*binlength),'all');
end
% plot(psth_TR)
bar(psth_TR)
% xlim([300 500])
xlabel('Time (sec)')
ylabel('# of spikes')
% 
% % linkaxes([ax1,ax2,ax3,ax4],'x')
% fh.WindowState = 'maximized';
% 
% 
% 
% 

%% Action potentials column 3
ax8 = subplot(6,5,3); title('S example APs'); hold on
plot(timeParams.tVec,vTC5(1,:),'k'); hold on
plot(timeParams.tVec,vTC5(20,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(60,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(80,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
zoom xon

ax9 = subplot(6,5,8); title('M example APs'); hold on
plot(timeParams.tVec,vTC5(1,:),'k'); hold on
plot(timeParams.tVec,vTC5(20,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(60,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(80,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
zoom xon

ax10 = subplot(6,5,13); title('D example APs'); hold on
plot(timeParams.tVec,vTC5(1,:),'k'); hold on
plot(timeParams.tVec,vTC5(20,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(60,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(80,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
zoom xon

ax11 = subplot(6,5,18); title('IC example APs'); hold on
plot(timeParams.tVec,vTC5(1,:),'k'); hold on
plot(timeParams.tVec,vTC5(20,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(60,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(80,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
zoom xon

ax12 = subplot(6,5,23); title('TC example APs'); hold on
plot(timeParams.tVec,vTC5(1,:),'k'); hold on
plot(timeParams.tVec,vTC5(20,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(60,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTC5(80,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
zoom xon

ax13 = subplot(6,5,28); title('Rt example APs'); hold on
plot(timeParams.tVec,vTR6(1,:),'k'); hold on
plot(timeParams.tVec,vTR6(10,:)+100*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTR6(20,:)+200*ones(1,numel(timeParams.tVec)),'k'); hold on
plot(timeParams.tVec,vTR6(40,:)+300*ones(1,numel(timeParams.tVec)),'k'); hold on
zoom xon
xlabel('Time (ms)')

%% -----------------------------LFP time sries
% I post synaptic layer 1
ax14 = subplot(6,5,4); 
title('S LFP: \SigmaEPSPs'); hold on
plot(timeParams.tVec,I1)
zoom xon
ylabel('EPSC (\muA)')
% xlim([0 timeParams.tEnd])
% I post synaptic layer 2
ax15 = subplot(6,5,9); 
title('M LFP: \SigmaEPSPs');
hold on
plot(timeParams.tVec,I2)
zoom xon
% xlabel('Time (ms)')
ylabel('EPSC (\muA)')
% xlim([0 timeParams.tEnd])
% I post synaptic layer 3
ax16 = subplot(6,5,14); 
title('D LFP: \SigmaEPSPs');
hold on
plot(timeParams.tVec,I3)
zoom xon
% xlabel('Time (ms)')
ylabel('EPSC (\muA)')
% xlim([0 timeParams.tEnd])
% I post synaptic layer 4
ax17 = subplot(6,5,19); 
title('CI LFP: \SigmaIPSPs');
hold on
plot(timeParams.tVec,I4)
zoom xon
% xlabel('Time (ms)')
ylabel('IPSC (\muA)')
% xlim([0 timeParams.tEnd])
% I post synaptic layer 5
ax18 = subplot(6,5,24); 
title('Vim LFP: \SigmaEPSPs');
hold on
plot(timeParams.tVec,I5)
zoom xon
% xlabel('Time (ms)')
ylabel('EPSC (\muA)')
% xlim([0 timeParams.tEnd])
% I post synaptic layer 6
ax19 = subplot(6,5,29); 
title('Rt LFP: \SigmaIPSPs');
hold on
plot(timeParams.tVec,I6)
zoom xon
% xlabel('Time (ms)')
ylabel('IPSC (\muA)')
% xlim([0 timeParams.tEnd])
xlabel('Time (ms)')

%% --------------------------------power spectrum
fs = 1000/timeParams.dt; frange = 1:200;

[p14,f1] = pwelch(I14,fs,fs/2,frange,fs);
[p24,f2] = pwelch(I24,fs,fs/2,frange,fs);
[p34,f3] = pwelch(I34,fs,fs/2,frange,fs);
[p4,f4] = pwelch(I4,fs,fs/2,frange,fs);
[p56,f5] = pwelch(I56,fs,fs/2,frange,fs);
[p6,f6] = pwelch(I6,fs,fs/2,frange,fs);

ax20 = subplot(6,5,5); 
title('S LFP PSD');
hold on
plot(f1,p14)
ylim([0 0.002])
ylabel('PSD')

ax21 = subplot(6,5,10); 
title('M LFP PSD');
hold on
plot(f2,p24)
ylim([0 0.002])
ylabel('PSD')

ax22 = subplot(6,5,15); 
title('D LFP PSD');
hold on
plot(f3,p34)
ylim([0 0.002])
ylabel('PSD')

ax23 = subplot(6,5,20); 
title('CI LFP PSD');
hold on
plot(f4,p4)
ylim([0 0.002])
ylabel('PSD')

ax24 = subplot(6,5,25); 
title('Vim LFP PSD');
hold on
plot(f5,p56)
ylim([0 0.002])
ylabel('PSD')

ax25 = subplot(6,5,30); 
title('Rt LFP PSD');
hold on
plot(f6,p6)
ylabel('PSD')
xlabel('Frequency (Hz)')
ylim([0 0.002])

%% % function fig_AllSpikes=Fig_spikes_all_netwrok(n_s,n_m,n_d,n_IN,n_rel,n_ret,n_t,spike_times,dt,nSim)
% 
% fig_AllSpikes=figure; set(gcf,'Visible','off');
% title('Thalamo-Cortical netwrok'); hold on
% 
% % obj1=patch([0 dt*nSim dt*nSim 0],[0 0 n_ret n_ret],'b'); hold on
% % obj2=patch([0 dt*nSim dt*nSim 0],[n_ret n_ret n_ret+n_rel n_ret+n_rel],'r'); hold on
% % obj3=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel n_ret+n_rel n_ret+n_rel+n_IN n_ret+n_rel+n_IN],'y'); hold on
% % obj4=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel+n_IN n_ret+n_rel+n_IN n_ret+n_rel+n_IN+n_d n_ret+n_rel+n_IN+n_d],'g'); hold on
% % obj5=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel+n_IN+n_d n_ret+n_rel+n_IN+n_d n_ret+n_rel+n_IN+n_d+n_m n_ret+n_rel+n_IN+n_d+n_m],'m'); hold on
% % obj6=patch([0 dt*nSim dt*nSim 0],[n_ret+n_rel+n_IN+n_d+n_m n_ret+n_rel+n_IN+n_d+n_m n_ret+n_rel+n_IN+n_d+n_m+n_s n_ret+n_rel+n_IN+n_d+n_m+n_s],'c'); hold on
% % transp=.03;
% % alpha(obj1,transp); alpha(obj2,transp); alpha(obj3,transp); alpha(obj4,transp); alpha(obj5,transp); alpha(obj6,transp);
% 
% for ik=1:n_t
%     for ij=1:length(spike_times{ik,:})
%         line([spike_times{ik}(ij) spike_times{ik}(ij)],[ik-1 ik],'LineWidth',1,'Color','k')
%     end
%     % axis([0 nSim 0 n_t]); hold on;
% end
% 
% % xlim([0 dt*nSim]); zoom xon
% ylim([0 n_t])
% 
% % line([0 nSim],[n_ret n_ret],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% % line([0 nSim],[n_ret+n_rel n_ret+n_rel],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% % line([0 nSim],[n_ret+n_rel+n_IN n_ret+n_rel+n_IN],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% % line([0 nSim],[n_ret+n_rel+n_IN+n_s n_ret+n_rel+n_IN+n_s],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% % line([0 nSim],[n_ret+n_rel+n_IN+n_s+n_m n_ret+n_rel+n_IN+n_s+n_m],'LineWidth',.5,'Color','k','LineStyle','-'); hold on
% 
% 
% % dim1 = [.7 .1]; str1= 'Thalamic Reticular Nucleus (TRN)'; annotation('textbox',dim1,'String',str1,'BackgoundColor','b','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% % dim2 = [.7 .2]; str2= 'Thalamocortical Relacy Nucleus (TCR)'; annotation('textbox',dim2,'String',str2,'BackgoundColor','r','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% % dim3 = [.7 .4]; str3= 'Cortical Inhibitory Neurons (IN)'; annotation('textbox',dim3,'String',str3,'BackgoundColor','b','HorizontalAlignment','center','VerticalAlignment','middle'); drawmow
% % dim4 = [.7 .6]; str4= 'Deep Cortical Layer (D)'; annotation('textbox',dim4,'String',str4,'BackgoundColor','g','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% % dim5 = [.7 .8]; str5= 'Middle Cortical Layer (M)'; annotation('textbox',dim5,'String',str5,'BackgoundColor','y','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% % dim6 = [.7 .9]; str6= 'Superficial Cortical Layer (S)'; annotation('textbox',dim6,'String',str6,'BackgoundColor','c','HorizontalAlignment','center','VerticalAlignment','middle'); drawnow
% 
% 
% xlabel('Time (ms)');
% ylabel('Neuron #'); hold on
% set(gca,'FontSize',12,'FontWeight','bold')