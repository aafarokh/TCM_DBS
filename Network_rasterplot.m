% Network Raster Plot

spikes_S1=zeros(nN/2,50000); spikes_M=zeros(nN,50000);
spikes_S2=zeros(nN/2,50000);
spikes_D1=zeros(.7*nN,50000); spikes_IC1=zeros(nN/2,50000);
spikes_D2=zeros(.3*nN,50000); spikes_IC2=zeros(nN/2,50000);
spikes_TC=zeros(nN,50000); spikes_TR=zeros((nN/2)-(nN/10),50000);

fh = figure;

title('Network Raster Plot'); hold on

nRT = (nN/2)-(nN/10);
for j=1:nRT
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vTR6(j,i)>neurons.vp(1,6)
            spikes_TR(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[j-1 j],'LineWidth',2,'Color',"#0072BD")
            %     drawnow
        end
    end
end

for j=1:nN
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vTC5(j,i)>neurons.vp(1,5)
            spikes_TC(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[nRT+j-1 j+nRT],'LineWidth',2,'Color','#4DBEEE')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vFS4(j,i)>neurons.vp(1,4)
            spikes_IC1(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[nN+nRT+j-1 j+nN+nRT],'LineWidth',2,'Color','#EDB120')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vLTS4(j,i)>neurons.vp(1,4)
            spikes_IC2(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[(3/2)*nN+nRT+j-1 j+(3/2)*nN+nRT],'LineWidth',2,'Color','#7E2F8E')
        end
    end
end

for j=1:30
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vIB3(j,i)>neurons.vp(1,3)
            spikes_D2(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[2*nN+nRT+j-1 j+2*nN+nRT],'LineWidth',2,'Color','#A2142F')
        end
    end
end

for j=1:70
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vRS3(j,i)>neurons.vp(1,3)
            spikes_D1(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[(0.3)*nN+2*nN+nRT+j-1 j+(0.3)*nN+2*nN+nRT],'LineWidth',2,'Color','#77AC30')
        end
    end
end

for j=1:nN
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vRS2(j,i)>neurons.vp(1,2)
            spikes_M(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[3*nN+nRT+j-1 j+3*nN+nRT],'LineWidth',2,'Color','#000000')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vIB1(j,i)>neurons.vp(1,1)
            spikes_S2(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[4*nN+nRT+j-1 j+4*nN+nRT],'LineWidth',2,'Color','#D95319')
        end
    end
end

for j=1:nN/2
    % for i=length(timeParams.tVec)-2*(1000/timeParams.dt):length(timeParams.tVec)
    for i=1:length(timeParams.tVec)
        if vRS1(j,i)>neurons.vp(1,1)
            spikes_S1(j,i) = 1;
            line([timeParams.tVec(i) timeParams.tVec(i)],[(9/2)*nN+nRT+j-1 j+(9/2)*nN+nRT],'LineWidth',2,'Color','#0000FF')
        end
    end
end


ylabel('Neuron #');

ylim([0 n_tot])

xlabel('Time (ms)');
ylabel('Neuron #'); hold on
set(gca,'FontSize',12,'FontWeight','bold')