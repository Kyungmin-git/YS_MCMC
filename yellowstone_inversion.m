%% Yellowstone Lake tremor inversion (Y1701, April 2018)

% Required input: YL_Y1701_ELZ_201804_selected_spectra.mat
% Required variables: frequency, representative_amplitude,
%                     representative_phase (phase in radians).
% Direct custom dependencies (and their own dependencies):
%   forwardmodelgaspocket.m, generatelpparameter.m,
%   generate_excitation.m, sigma_propose.m
% Toolboxes: Signal Processing; Statistics and Machine Learning.
%
% Model vector: [T, Q, R, L, -log10(kappa), -log10(D),
%                log10(phi), r_obs]. Parameter conventions are inherited
% from the supplied forward model and proposal functions.
% Default run: 24 million iterations; requires several GB of RAM.
% Output figures and workspace checkpoints are written to Current Folder.
% The sampler is a single chain; the inactive hot-chain code was removed.
% Randomness uses MATLAB's current RNG state. Set rng(seed) before running


%% Input and observation window
input_file = 'YL_Y1701_ELZ_201804_selected_spectra.mat';
assert(isfile(input_file), 'Missing input file: %s', input_file);
load(input_file, 'frequency', 'representative_amplitude', 'representative_phase');
assert(exist('frequency', 'var') == 1 && ...
exist('representative_amplitude', 'var') == 1 && ...
exist('representative_phase', 'var') == 1, ...
'Input must contain frequency, representative_amplitude, and representative_phase.');
frequency = frequency(:).';
representative_amplitude = representative_amplitude(:).';
representative_phase = representative_phase(:).';
assert(numel(representative_amplitude) == 901 && ...
numel(frequency) == 901 && numel(representative_phase) >= 750, ...
'Expected 901 amplitude/frequency bins and at least 750 phase bins.');
zz=100;
kk=9;
segment_duration = 30;  % seconds
global LPstartt LPendt tau
LPstartt=0;
LPendt=30;
tau=LPendt-LPstartt;
%% Observation smoothing and uncertainties
data=representative_amplitude; % Using representative spectrum
N=length(data);
filtered_data=zeros(1,N);
window_width=30;
for i = 1:N
    start_idx = i; % Start of the window
end_idx = min(i + window_width - 1, N); % Ensure we don't exceed the vector length
filtered_data(i) = mean(data(start_idx:end_idx));
end
theta=representative_phase; % representative phase
noise_ratio=0.33;
Uncertainty_spect=abs(noise_ratio*filtered_data);
Uncertainty_spect_phase=abs(noise_ratio*data);
Covariance5=diag(Uncertainty_spect(16:900).^2);
inv_variance_row=1./Uncertainty_spect(16:900).^2;
icov5=diag(inv_variance_row);
sigma_row_angle=angle(Uncertainty_spect_phase*1i+data(1:901));
assert(all(isfinite(Uncertainty_spect(16:900))) && ...
all(Uncertainty_spect(16:900) > 0) && ...
all(isfinite(sigma_row_angle(16:750))) && ...
all(sigma_row_angle(16:750) > 0), ...
'Fitted bins must have finite, positive amplitude and phase uncertainties.');
Inv_phase=1./sigma_row_angle(16:750).^2;
ipcov5=diag(Inv_phase);
%% Initial physical model
T0=273.15+188;
Q0=4+2*rand;
R0=70+20*rand;
L0=12+4*rand;
exp_kappa0=10+4*rand;kappa0=10^(-exp_kappa0);
exp_D0=0.3+0.8*rand; D0=10^(-exp_D0);
exp_phi_n0=0.7*rand; phi_n0=10^(exp_phi_n0);
r_obs0=15;
YE=10^9;
nu=0.45;
sigma_k_log=-1.6;
current_sigma_k=sigma_k_log;
figure(7)
plot(frequency,filtered_data)
xlabel('Frequency (Hz)')
ylabel('Displacement (m*s)')
set(gca,'FontSize',18)
xlim([0 30])
%% Initial stochastic gas excitation
NNinit=1000;
alpha1=1;
qn_n=gamrnd(alpha1,1,1,NNinit);
qn_n=qn_n/sum(qn_n);
qn=Q0*tau*qn_n;
t0=rand(1,NNinit)*tau;
m0=[T0,Q0,R0,L0,exp_kappa0,exp_D0,exp_phi_n0,r_obs0];
%% Sampler settings and storage
n = length(m0);
niter=24000000;
mout = zeros(niter, n);
global stepsize stepsize2
stepsize=0.003*[1.5*10,0.4*25,1.5*10,1*20,0.5*4,0.3*3.5,1*3,0.5*50];
stepsize2=0.003*1;
current_sigma_k=sigma_k_log;
sigma_store=zeros(1,niter);
sigma_store(1)=current_sigma_k;
mout(1, :) = m0;
current_N=NNinit;
Q=m0(2);
currentqn=qn;
qn_initial=qn;
[Poscinit,cuuspect,cuuspect_raw,cusynangle,~,A_resin,A_excin,A_pathin,~,G_test]=forwardmodelgaspocket(m0,nu,YE,t0,qn,current_sigma_k);
data_1=abs(cuuspect_raw);
data_2=abs(cuuspect);
N=length(data_1);
filtered_cuuspect=zeros(1,N);
window_width=30;
for i = 1:N
    start_idx = i; % Start of the window
end_idx = min(i + window_width - 1, N); % Ensure we don't exceed the vector length
filtered_cuuspect(i) = mean(data_1(start_idx:end_idx));
end
cuuspect=cuuspect_raw; % Initialize the stored spectrum consistently.
current=m0;
NNmax= 10000;
lMAP = -Inf;
llMAP = -inf;
mMAP = current;
nacc = 0;
llcandidate_st=zeros(1,niter/40000);
lpcandidate_st=zeros(1,niter/40000);
llcurrent_st=zeros(1,niter/40000);
llcurrent_p_st=zeros(1,niter/40000);
lpcurrent_st=zeros(1,niter/40000);
acceptance=zeros(1,niter);
gastype_st=zeros(1,niter/2);
current_N_st=zeros(1,niter);
t0_out=zeros(NNmax,niter/40000);
qn_out=zeros(NNmax,niter/40000);
LPspect_st=zeros(900,niter/40000);
A_exc_current=A_excin;
cusynangle2=cusynangle(:).';
rk=filtered_cuuspect(16:900)-filtered_data(16:900);
diff=abs(theta(16:750)-cusynangle2(16:750));
rkp=min(diff,2*pi-diff);
llcurrent=-(1/2)*rk*icov5*rk';
llcurrent_p=-(1/2)*rkp*ipcov5*rkp';
lpcurrent=0;
rk_init=rk;
current_N_st(1)=NNinit;
Pmax=zeros(1,niter);
Pmax(1)=max(Poscinit)-min(Poscinit);
sampling_frequency=200;
order=4;low_frequency=1;high_frequency=35;
[z,p]=butter(order/2,[low_frequency high_frequency]/(sampling_frequency/2));
Poscinit2=filtfilt(z,p,Poscinit);
Pmaxf=zeros(1,niter);
Pmaxf(1)=max(Poscinit2)-min(Poscinit2);
mu_g=1.5e-5;                                                                  %gas viscosity
M=0.018;
Rg=8.3145;
lpcurrent=local_log_prior(m0);
currentt0=t0;
currentqn_n=qn_n;
llMAP=llcurrent+llcurrent_p;
mlMAP=current;
mlMAPt0=currentt0;
mlMAPqn=qn;
sigma_k_log_MAP=current_sigma_k;
mMAPP=Pmax(1);
%% Metropolis-Hastings / reversible-jump sampling
warning('off', 'MATLAB:colon:nonIntegerIndex'); %suppress warning message
for k = 2:niter
    if mod(k, 10000) == 0
        fprintf('Iteration %d / %d; acceptance %.3f\n', k, niter, nacc / (k-1));
    end
    global r rb
    gastype = NaN; % Only defined for excitation proposals.
    r=rem(k,5);
    ra=k/2;
    rb=rem(ra,9);
    if r==1 || r==2
        candidate = generatelpparameter(current);
        candidatet0=currentt0;
        candidateqn_n=currentqn_n;
        logJacobian=0;
        logproposal=0;
        candidate_sigma_k=current_sigma_k;
    elseif r==3 || r==4
        [candidatet0,candidateqn_n,gastype,logJacobian,logproposal]=generate_excitation(currentt0,currentqn_n);
        candidate=current;
        candidate_sigma_k=current_sigma_k;
    else
        [candidate_sigma_k]=sigma_propose(current_sigma_k);
        candidatet0=currentt0;
        candidateqn_n=currentqn_n;
        candidate=current;
        logJacobian=0;
        logproposal=0;
    end
    if rem(k,2)==1
        gastype_st(1,(k-1)/2)=gastype;
    else
    end
    m=candidate;
    t0=candidatet0;
    candidate_N=length(t0);
    lpcandidate = local_log_prior(m);
    D2=10^(-m(6));
    L2=m(4);
    Pex=1.5e6;
    candidateqn=candidateqn_n*candidate(2)*tau;
    [Posc,~,modelspect,casynangle2,syntheticf,A_res,A_exc,A_path,~]=forwardmodelgaspocket(m,nu,YE,t0,candidateqn,candidate_sigma_k);
    casynangle=casynangle2(:).';
    order=4;low_frequency=1;high_frequency=35;
    [z,p]=butter(order/2,[low_frequency high_frequency]/(sampling_frequency/2));
    Posc2=filtfilt(z,p,Posc);
    data=abs(modelspect);
    N=length(data);
    filtered_modelspect=zeros(1,N);
    window_width=30;
    for i = 1:N
        start_idx = i; % Start of the window
    end_idx = min(i + window_width - 1, N); % Ensure we don't exceed the vector length
    filtered_modelspect(i) = mean(data(start_idx:end_idx));
end
rk=(filtered_modelspect(16:900)-filtered_data(16:900));
diff=abs(theta(16:750)-casynangle(16:750));
if size(casynangle,1)>2
    break
end
rkp=min(diff,2*pi-diff);
llcandidate=(-1/2)*rk*icov5*rk';
llcandidate_p=(-1/2)*rkp*ipcov5*rkp';
Temp=1;
logalpha = lpcandidate+ llcandidate*(1/Temp) +llcandidate_p*(1/Temp) - lpcurrent - llcurrent*(1/Temp) -llcurrent_p*(1/Temp)+logproposal+logJacobian;
if (logalpha > 0)
    logalpha = 0;
end
logt = log(rand());
if (logt < logalpha)
    current = candidate;
    currentt0 =[];
    currentt0 = candidatet0;
    A_exc_current=A_exc;
    currentqn_n = candidateqn_n;
    current_sigma_k = candidate_sigma_k;
    current_N = candidate_N;
    cuuspect = modelspect;
    lpcurrent = lpcandidate;
    llcurrent = llcandidate;
    llcurrent_p = llcandidate_p;
    acceptance(k) =1;
    nacc = nacc + 1;
    Pmax(k)=max(Posc)-min(Posc);
    Pmaxf(k)=max(Posc2)-min(Posc2);
else
    acceptance(k) =0;
    lpcurrent = lpcurrent;
    llcurrent=  llcurrent;
    Pmax(k)=Pmax(k-1);
    Pmaxf(k)=Pmaxf(k-1);
end
mout(k,:) = current;
sigma_store(k)=current_sigma_k;
accrate = nacc / (k-1);
current_N_st(k)=current_N;
if mod(k,40000)==0
    t0_out(1:current_N,k/40000)=currentt0;
    currentqn=currentqn_n*current(2)*tau;
    qn_out(1:current_N,k/40000)=currentqn;
    llcandidate_st(1,k/40000)=llcandidate;
    lpcandidate_st(1,k/40000)=lpcandidate;
    llcurrent_st(1,k/40000)=llcurrent;
    lpcurrent_st(1,k/40000)=lpcurrent;
    llcurrent_p_st(1,k/40000)=llcurrent_p;
    LPspect_st(1:900,k/40000)=cuuspect(1:900);
    % cuuspect holds the current raw model spectrum after accepted steps.
    LPspectraw_st(1:900,k/40000)=cuuspect(1:900);
else
end
if (llcandidate+llcandidate_p > llMAP)
    llMAP = llcandidate+llcandidate_p;
    llMAP_m = llcandidate;
    mlMAP = candidate;
    mlMAPt0 = candidatet0;
    candidateqn=candidateqn_n*candidate(2)*tau;
    mlMAPqn = candidateqn;
    cuuspect_raw_MAP = modelspect;
    sigma_k_log_MAP=candidate_sigma_k;
    mMAPP=max(Posc)-min(Posc); % Pressure of the best-likelihood candidate.
end
if mod(k,400000)==0 || k==5000
    k
    figure(zz+17)
    plot(llcurrent_st)
    h=figure(zz+17);
    savefig(h,sprintf('llcurrent_%d.fig',kk),'compact')
    saveas(h,sprintf('llcurrent_%d.jpg',kk))
    close(h)
    h40=figure(zz+14);
    subplot(3,3,1)
    h11=histogram(mout(1001:k,1)-273.15,'Normalization','Probability','BinWidth',0.5);
    xlabel('\circC')
    xlim([185 190])
    title('Temperature')
    set(gca,'YTickLabel',[]);
    h11.FaceColor=[0 0.4470 0.7410];
    h11.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,2)
    h12=histogram(mout(1001:k,2),'Normalization','Probability','BinWidth',0.5);
    xlim([0 10])
    set(gca,'YTickLabel',[]);
    xlabel('Kg/s')
    title('Gas flow rate')
    h12.FaceColor=[0 0.4470 0.7410];
    h12.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,3)
    h13=histogram(mout(1001:k,3),'Normalization','Probability','BinWidth',0.1);
    xlabel('m')
    xlim([35 120])
    title('Conduit Radius')
    set(gca,'YTickLabel',[]);
    h13.FaceColor=[0 0.4470 0.7410];
    h13.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,4)
    h14=histogram(mout(1001:k,4),'Normalization','Probability','BinWidth',1);
    xlim([0 20])
    xlabel('m')
    title('Cap Thickness')
    set(gca,'YTickLabel',[]);
    h14.FaceColor=[0 0.4470 0.7410];
    h14.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,5)
    [~,edges] = histcounts(log10(10.^(-mout(1001:k,6))),50);
    h15=histogram(10.^(-mout(1001:k,6)),10.^edges,'Normalization','Probability');
    xlabel('m')
    xlim([10^(-2) 5])
    title('Gas Pocket Thick')
    set(gca,'YTickLabel',[]);
    set(gca,'xscale','log')
    xticks([10^(-2) 10^(-1) 10^(0)])
    h15.FaceColor=[0 0.4470 0.7410];
    h15.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,6)
    [~,edges] = histcounts(log10(10.^(-mout(1001:k,5))),50);
    h16=histogram(10.^(-mout(1001:k,5)),10.^edges,'Normalization','Probability');
    xlabel('m^2')
    xlim([10^(-17) 5*10^(-7)])
    xticks([10^(-15) 10^(-11) 10^(-7)])
    xlabel('m^{2}')
    title('Permeability')
    set(gca,'xscale','log')
    set(gca,'YTickLabel',[]);
    h16.FaceColor=[0 0.4470 0.7410];
    h16.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,7)
    [~,edges] = histcounts(log10(10.^(mout(1001:k,7))),50);
    h15=histogram(10.^(mout(1001:k,7)),10.^edges,'Normalization','Probability');
    xlabel('%')
    title('Porosity')
    xlim([0 30])
    set(gca,'YTickLabel',[]);
    h15.FaceColor=[0 0.4470 0.7410];
    h15.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,8)
    h17=histogram(mout(1001:k,8),'Normalization','Probability','BinWidth',1);
    xlabel('m')
    title('Distance')
    xlim([0 50])
    set(gca,'YTickLabel',[]);
    h17.FaceColor=[0 0.4470 0.7410];
    h17.EdgeColor=[0 0.4470 0.7410];
    subplot(3,3,9)
    Pmax3 = Pmax(1:k);
    Pmax3(Pmax3 <= 0) = [];   % Remove negative values
    Pmax3(Pmax3 > 10^7) = []; % Remove values greater than 10^7
    log_edges = logspace(log10(min(Pmax3)), log10(max(Pmax3)), 20); % 20 bins
    h20=histogram(Pmax3, log_edges, 'Normalization', 'Probability');
    xlabel('P (Pa)')
    title('Max Overpressure')
    yticklabels([]);
    set(gca, 'XScale', 'log'); % Ensure x-axis is log scale
    xticks([10^(1) 10^(3) 10^(5)])
    xlim([10^(0) 10^(5.5)])
    h20.FaceColor=[0.8500 0.3250 0.0980];
    h20.EdgeColor=[0.8500 0.3250 0.0980];
    h40=figure(zz+14);
    savefig(h40,sprintf('MCMC_source_%d.fig',kk),'compact')
    saveas(h40,sprintf('MCMC_source_%d.jpg',kk))
    [PoscMAP,~,spectMAPf,angleMAPf,syntheticMAPf,~,~,~,~]=forwardmodelgaspocket(mlMAP,nu,YE,mlMAPt0,mlMAPqn,sigma_k_log_MAP);
    syntheticMAPf=syntheticMAPf-mean(syntheticMAPf);
    order=4;sampling_frequency=200;low_frequency=1;high_frequency=35;
    [z,p]=butter(order/2,[low_frequency high_frequency]/(sampling_frequency/2));
    syntheticMAPf=filtfilt(z,p,syntheticMAPf);
    filtered_MAP=zeros(1,N);
    window_width=30;
    data=spectMAPf;
    for i = 1:N
        start_idx = i; % Start of the window
    end_idx = min(i + window_width - 1, N); % Ensure we don't exceed the vector length
    filtered_MAP(i) = mean(data(start_idx:end_idx));
end
figure(zz+18)
plot(frequency(1:900),filtered_data(1:900),'r',frequency(1:900), filtered_MAP(1:900),'b')
xlabel('Frequency (Hz)','FontSize',18)
ylabel('u_z (m*s)','FontSize',18)
title('Ground displacement')
xlim([0 30])
ylim([0 6*10^(-8)])
set(gca,'FontSize',18,'fontweight','bold')
h18=figure(zz+18);
savefig(h18,sprintf('fit_synthetic_%d.fig',kk),'compact')
saveas(h18,sprintf('fit_synthetic_%d.jpg',kk))
time1=(0:numel(syntheticMAPf)-1)/80; % Forward-model waveform sampling rate.
figure(zz+25)
plot(time1,syntheticMAPf,'r')
xlabel('time (s)','FontSize',18)
ylabel('Amplitude (m)','FontSize',18)
title('Ground displacement')
xlim([0 30])
set(gca,'FontSize',18,'fontweight','bold')
h22=figure(zz+25);
savefig(h22,sprintf('fit_timedomain_%d.fig',kk),'compact')
saveas(h22,sprintf('fit_timedomain_%d.jpg',kk))
figure(zz+41)
histogram(sigma_store(1001:k),'Normalization','Probability')
xlim([-2 max(sigma_store(1001:k))])
title('sigma_histogram')
set(gca,'FontSize',18,'fontweight','bold')
h41=figure(zz+41);
savefig(h41,sprintf('Sigma_store_%d.fig',kk),'compact')
saveas(h41,sprintf('Sigma_store_%d.jpg',kk))
figure(zz+34)
plot(frequency(1:900),filtered_data(1:900),'b',frequency(1:900),filtered_MAP(1:900),'r')
xlabel('Frequency (Hz)','FontSize',18)
ylabel('u_z (m*s)','FontSize',18)
title('Ground displacement')
xlim([0 30])
ylim([0 4*10^(-8)])
set(gca,'FontSize',18,'fontweight','bold')
h34=figure(zz+34);
savefig(h34,sprintf('fit_synthetic2_%d.fig',kk),'compact')
saveas(h34,sprintf('fit_synthetic2_%d.jpg',kk))
figure(zz + 23)
hold on
f=frequency;
scatter(f(31:300), angleMAPf(31:300), 60, 'r', 'filled')
scatter(f(31:300), theta(31:300), 60, 'b', 'filled')
for i = 31:300
    diff = abs(angleMAPf(i) - theta(i));
    if diff > pi
        if angleMAPf(i) > theta(i)
            plot([f(i), f(i)], [angleMAPf(i), pi], 'k-', 'LineWidth', 1.5)
            plot([f(i), f(i)], [-pi, theta(i)], 'k-', 'LineWidth', 1.5)
        else
            plot([f(i), f(i)], [theta(i), pi], 'k-', 'LineWidth', 1.5)
            plot([f(i), f(i)], [-pi, angleMAPf(i)], 'k-', 'LineWidth', 1.5)
        end
    else
        plot([f(i), f(i)], [angleMAPf(i), theta(i)], 'k-', 'LineWidth', 1.5)
    end
end
yticks([-pi 0 pi])
yticklabels({'-\pi', '0', '\pi'})
xlabel('Frequency (Hz)', 'FontSize', 18)
ylabel('Phase (radian)', 'FontSize', 18)
title('Phase Fit', 'FontSize', 18)
ylim([-pi +pi])
set(gca, 'FontSize', 18, 'FontWeight', 'bold')
h23 = figure(zz + 23);
savefig(h23, sprintf('fit_phase_%d.fig', kk), 'compact')
saveas(h23, sprintf('fit_phase_%d.jpg', kk))
hold off
figure(zz+20)
plot((0:numel(PoscMAP)-1)/80,PoscMAP)
xlabel('time (s)','FontSize',18)
ylabel('Pressure (Pa)','FontSize',18)
title('Pressure oscillation')
set(gca,'FontSize',18,'fontweight','bold')
h20=figure(zz+20);
savefig(h20,sprintf('Pressureosc_%d.fig',kk),'compact')
saveas(h20,sprintf('Pressureosc_%d.jpg',kk))
h21=figure(zz+21);
scatter(mlMAPt0,mlMAPqn)
xlabel("time(s)")
ylabel("Q (kg)")
set(gca,'FontSize',18,'fontweight','bold')
savefig(h21,sprintf('MCMC_t0qn_best_%d.fig',kk),'compact')
saveas(h21,sprintf('MCMC_t0qn_best_%d.jpg',kk))
close(h21)
acceptance1=acceptance(1:5:k);
acceptance2=acceptance(2:5:k);
acceptance3=acceptance(3:5:k);
acceptance4=acceptance(4:5:k);
acceptance5=acceptance(5:5:k);
accept1in=sum(acceptance1)/length(acceptance1);
accept2in=sum(acceptance2)/length(acceptance2);
accept3in=sum(acceptance3)/length(acceptance3);
accept4in=sum(acceptance4)/length(acceptance4);
accept5in=sum(acceptance5)/length(acceptance5);
Pmax3=Pmax(2:k);
Pmax3(Pmax3<0)=[];
Pmax3(Pmax3>10^7)=[];
h73=figure(zz+73);
histogram(Pmax3,'Normalization','Probability')
xlabel('P (Pa)')
title('overpressure')
title('Max Overpressure')
yticklabels([])
set(gca,'FontSize',18)
savefig(h73,sprintf('MCMC_overp2_%d.fig',kk),'compact')
saveas(h73,sprintf('MCMC_overp2_%d.jpg',kk))
close(h73)
h74=figure(zz+74);
histogram(Pmax3,'Normalization','Probability')
xlabel('P (Pa)')
title('overpressure')
xlim([0 10^6])
title('Max Overpressure')
yticklabels([])
set(gca,'FontSize',18)
savefig(h74,sprintf('MCMC_overp3_%d.fig',kk),'compact')
saveas(h74,sprintf('MCMC_overp3_%d.jpg',kk))
close(h74)
qn_ave=m(2)*tau/length(t0);
qns=sort(qn,'descend');
currentqn=currentqn_n*current(2)*tau;
h25=figure(zz+25);
sz=80;
scatter(currentt0,currentqn,sz,'k')
xlabel("time")
ylabel("Q (kg)")
ax=gca;
set(gca,'FontSize',20);
savefig(h25,sprintf('MCMCt0qn_%d.fig',kk),'compact')
saveas(h25,sprintf('MCMCt0qn_%d.jpg',kk))
h26=figure(zz+26);
hold on
plot(current_N_st(1:k),'color',[0 0.4470 0.7410],'LineWidth',1.2)
xlabel("iteration")
ylabel("NN")
set(gca,"Fontsize",18);
hold off
box on
saveas(h26,sprintf('gas_excitation1_%d.jpg',kk))
close(h26)
h37=figure(zz+37);
plot(llcurrent_st)
xlabel("")
ylabel("")
set(gca,'FontSize',18,'fontweight','bold')
saveas(h37,sprintf('llcurrent_%d.jpg',kk))
close(h37)
h38=figure(zz+38);
plot(llcurrent_p_st)
xlabel("")
ylabel("")
set(gca,'FontSize',18,'fontweight','bold')
saveas(h38,sprintf('llcurrent_p_%d.jpg',kk))
close(h38)
else
end
if mod(k,200000)==0 || k==10000
    save(sprintf('syntheticresultadd_%d.mat',kk), '-v7.3')
else
end
end

%% Local prior (shared by initialization and proposals)
function lp = local_log_prior(m)
if (m(1)>=458.15) && (m(1)<=463.15) && ...
    (m(2)>=0.1) && (m(2)<=100) && ...
    (m(3)>=35) && (m(3)<=120) && ...
    (m(4)>=10) && (m(4)<=20) && ...
    (m(5)>=(log10(2)+6)) && (m(5)<=17) && ...
    (m(6)>=-log10(5)) && (m(6)<=2) && ...
    (m(7)>=-1) && (m(7)<=(1+log10(3))) && ...
    (m(8)>=0) && (m(8)<=50) && (m(8)<=m(3))
    sigma7_upper = ((1+log10(3)) - (1+log10(2))) / 2.5;
    lp = ...
    -(m(4)-15)^2/(2*2^2) ...                         % m4 Gaussian
    -min(m(7),0)^2/(2*0.25^2) ...                    % m7 lower taper
    -max(m(7)-(1+log10(2)),0)^2/(2*sigma7_upper^2) ... % m7 upper taper
    -min(m(6)+log10(3),0)^2/(2*0.05^2) ...           % m6 lower taper
    -min(m(3)-50,0)^2/(2*5^2) ...                    % m3 lower taper
    -max(m(3)-100,0)^2/(2*5^2);                      % m3 upper taper
else
    lp = -Inf;
end
end
