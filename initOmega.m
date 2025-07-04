function omega = initOmega(nN,nL)
omega = zeros(nN,nL,nN,nL);
omegaIntensity = 1;% for adjusting inter-connections intensity
fac_N=1.0; fac_PD=5; %for 100 neurons
ini=0.0; fin=1; int=fin-ini;
fac=fac_PD;
%COUPLING STRENGTHs within each structure (The same in Normal and PD)
r_s     = ini + int*rand(nN,nN);
r_m     = ini + int*rand(nN,nN);
r_d     = ini + int*rand(nN,nN);
r_ins   = ini + int*rand(nN,nN);
r_ret   = ini + int*rand(nN,nN);
r_rel   = ini + int*rand(nN,nN); 

%COUPLING STRENGTHs between structures (PD)
% we use index 1 to 6 for layers S,M,D,INs,Ret,Rel
%S
aee_ss  = -1e1/fac;     omega(:,1,:,1) = aee_ss*r_s;            %S to S(was -1e-2 for IEEE paper)
aee_sm  = 1e1/fac;      omega(:,2,:,1) = aee_sm*r_s;            %M to S couplings
aee_sd  = 5e2/fac;      omega(:,3,:,1) = aee_sd*r_s;            %D to S couplings
aei_sINs= -5e2/fac;     omega(:,4,:,1) = aei_sINs*r_s;          %INs to S couplings
aei_sRet= 0/fac;        omega(:,5,:,1) = aei_sRet*r_s;          %Ret. to S couplings
aee_sRel= 0/fac;        omega(:,6,:,1) = aee_sRel*r_s;          %Rel. to S couplings

%M
aee_ms=3e2/fac;         omega(:,1,:,2)= aee_ms*r_m;             %S to M couplings
aee_m=-1e1/fac;         omega(:,2,:,2)= aee_m*r_m;              %M to M (was -1e-2 for IEEE paper)
aee_md=0/fac;           omega(:,3,:,2)= aee_md*r_m;             %D to M couplings
aei_mINs=-3e2/fac;      omega(:,4,:,2)= aei_mINs*r_m;           %INs to M couplings
aei_mRet=0/fac;         omega(:,5,:,2)= aei_mRet*r_m;           %Ret. to M couplings
aee_mRel=0/fac;         omega(:,6,:,2)= aee_mRel*r_m;           %Rel. to M couplings

%D
aee_ds=3e2/fac;         omega(:,1,:,3)=aee_ds*r_d;              %S to D couplings
aee_dm=0/fac;           omega(:,2,:,3)=aee_dm*r_d;              %M to D couplings
aee_d=-1e1/fac;         omega(:,3,:,3)= aee_d*r_d;               %D to D (was -1e-2 for IEEE paper)
aei_dINs=-7.5e3/fac;    omega(:,4,:,3)=aei_dINs*r_d;          %INs to D couplings
aei_dRet=0/fac;         omega(:,5,:,3)=aei_dRet*r_d;          %Ret. to D couplings
aee_dRel=1e1/fac;       omega(:,6,:,3)=aee_dRel*r_d;          %Rel. to D couplings

%INs
aie_inss=2e2/fac;       omega(:,1,:,4)=aie_inss*r_ins;     %S to INs couplings
aie_insm=2e2/fac;       omega(:,2,:,4)=aie_insm*r_ins;     %M to INs couplings
aie_insd=2e2/fac;       omega(:,3,:,4)=aie_insd*r_ins;     %D to INs couplings
aii_INs=-5e2/fac;       omega(:,4,:,4)=aii_INs*r_ins;        %INs to INs
aii_InsRet=0/fac;       omega(:,5,:,4)=aii_InsRet*r_ins; %Ret to INs couplings
aie_InsRel=1e1/fac;     omega(:,6,:,4)=aie_InsRel*r_ins; %Rel to INs couplings

%Ret.
aie_rets=0/fac;         omega(:,1,:,5)= aie_rets*r_ret;     %S to Ret couplings
aie_retm=0/fac;         omega(:,2,:,5)= aie_retm*r_ret;     %M to Ret couplings
aie_retd=7e2/fac;       omega(:,3,:,5)= aie_retd*r_ret;     %D to Ret couplings
aii_RetIns=0/fac;       omega(:,4,:,5)= aii_RetIns*r_ret; %INs to Ret couplings
aii_ret=-5e1/fac;       omega(:,5,:,5)= aii_ret*r_ret;        %Ret to Ret
aie_RetRel=1e3/fac;     omega(:,6,:,5)= aie_RetRel*r_ret; %Rel. Ret INs couplings

%Rel.
aee_rels=0/fac;         omega(:,1,:,6)=aee_rels*r_rel;       %S to Rel couplings
aee_relm=0/fac;         omega(:,2,:,6)=aee_relm*r_rel; 	    %M to Rel couplings
aee_reld=7e2/fac;       omega(:,3,:,6)=aee_reld*r_rel;       %D to Rel couplings
aei_RelINs=0/fac;       omega(:,4,:,6)=aei_RelINs*r_rel;   %INs to Rel couplings
aei_RelRet=-5e2/fac;    omega(:,5,:,6)=aei_RelRet*r_rel;   %Ret to Rel couplings
aee_rel=0/fac;          omega(:,6,:,6)= aee_rel*r_rel;        %Rel to Rel


 %normalize omega
 %omega = rand(nN,nL,nN,nL);
 omegaAbs= abs(omega);
 sumW = sum(omegaAbs,3);
 sumW = sum(sumW ,4);
 omega = omegaIntensity*omega./sumW*500;
 omega(isnan(omega)) = 0;
 omega(isinf(omega)) = 0;


