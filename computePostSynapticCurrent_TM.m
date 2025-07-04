function neurons = computePostSynapticCurrent_TM(neurons,timeParams,indSpike)
neurons.r = neurons.r + timeParams.dt*(-neurons.r./neurons.tauf) + neurons.U.*(1-neurons.r).*indSpike;
neurons.x = neurons.x + timeParams.dt*((1-neurons.x)./neurons.taud) - neurons.r.*neurons.x.*indSpike;
neurons.I = neurons.I + timeParams.dt*(-neurons.I./neurons.taus) + neurons.A.*neurons.r.*neurons.x.*indSpike;
neurons.I(isnan(neurons.I))=0;