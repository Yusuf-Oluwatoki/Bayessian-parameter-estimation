close all
clear all 
clc
%% Adding the toolbox paths.
addpath('./wafo_2017')

addpath('./mcmcstat-master')

%% Setting parameters for the spectrum realisation. 
N = 2048; Nspec = 257;
% Number of points for a realisation and for the spectrum function.
p = [3,7]; 
% Significant wave height and peak period of the spectrum.

% rng('default')
Spec = pmspec([],p); Spec.h = Inf;
 % allowing the toolbox to define the range for angular frequency and
 % gettting a spectrum realisation.
 
%figure
%plotspec(Spec)

% plot of the realised spectrum.

% R1 = spec2cov(Spec,1);%spectrum to covariance from toolbox.
%opt=simoptset(opt,'dt',0.25,'du',0.25);
opt = simoptset('Nt',N,'Nu',20, 'lalpha', 0,'ffttype','ffttime');
% Spectrum realisation and setting options for wave realisation.

% Get a realisation from the toolbox.
[w,x] = spec2ldat(Spec,opt);
type='time';
[L,Lsmooth]=ldat2lwav(w,x,type,[],0.25);

% figure
% plot(Lsmooth.t,Lsmooth.Z);      %axis([0 50 -6 6])
% ylabel("Surface elevation", 'interpreter', 'latex')
% xlabel("Time", 'interpreter', 'latex')
% axis tight
% grid on
% ax = gca; % current axes
% ax.FontSize = 10;

realisation = Lsmooth.Z';
realisation = realisation - mean(realisation);

times = Lsmooth.t;
% realisation time steps.

R1 = spec2cov(Spec,1);
%spectrum to covariance from toolbox.
n = numel(times); %number of time steps
differ = times(2)-times(1);
% time point difference for interpolation.

k = zeros(n,1);
% Loop over the time points of the realisation and interpolate the
% covariance values from the vectors (time points R1.t and covaraince values R1.R) given by the toolbox.
     for i=1:n
         q = interp1(R1.t,R1.R,(i-1)*differ,'linear','extrap');
         k(i) = q;
     end
     % Construct the covariance matrix s for the wave process
     % and add small diagonal matrix to it (i.e add white noise process).
     % Otherwise the matrix would be singular.
     ccov=toeplitz(k,k');
     s = eye(n)*0.3;
     ccov = ccov+s; ccov = 0.5*(ccov' +ccov); 
     iCov = (ccov)\eye(n);
     
     Rp = cholcov(ccov);
     
     % Realisations from covariance matrix
     X=zeros(1,n);
     m=10;
for i=1:m
    X(i,:)=Rp*randn(n,1);
end
%%
figure
plot(X(1,:)','k');
ylabel("Sea elevation [m]", 'interpreter', 'latex')
xlabel("Time [s]", 'interpreter', 'latex')
title("Wave realisation", 'interpreter', 'latex')
axis tight
x0=10; y0=10; width=550; height=350;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 22;

% Plot of first realisation out of m realisations.


%% Maximum likelihood estimation of parameters.

omega=linspace(0,45,Nspec);
% Range for spectrum's angular speed.
for i=1:m
max_lik_este(i,:) = fminsearch(@(p) -loglike(p,times,X(i,:)'),[p(1),p(2)])
end
%%
% Estimating the parameters of the first m wave realisations.
logposterior=@(p) loglike(p,times,X(10,:)')+log((p(1)-0.03)^2/(2*0.05^2)) +log((p(2)-0.03)^2/(2*0.05^2));

maximum_li=fminsearch(@(p) -loglike(p,times,X(10,:)'),[p(1),p(2)]')
maximum_ap=fminsearch(@(x) -logposterior(x),[p(1),p(2)]')
%maximum_app=fminsearch(@(p) -npost(X,p),[1,1]')

%%

figure
plotspec(Spec)
hold on
for i=1:m
plotspec(pmspec([],max_lik_este(i,:)),'r');
hold on
end

grid on;
ax = gca; % current axes
ax.FontSize = 10;
legend('peak period','Spectrum', 'Estimate-peak period','Estimate', 'interpreter', 'latex')
title('Spectrum as a fucntion of angular speed $\omega$', 'interpreter', 'latex')
%%
% a=arrayfun(@(w) S(w,p),omega);
% 
% figure
% plot(omega,a)
% hold on
% plot([omega(find(a==max(a))) omega(find(a==max(a)))],[ 0 max(a)], '--b')
% hold on
% for i=1:m
% aa=arrayfun(@(w) S(w,max_lik_este(i,:)),omega);
% plot(omega,aa)
% hold on
% plot([omega(find(aa==max(aa))) omega(find(aa==max(aa)))], [0 max(aa)], '--r')
% hold on
% end
% xlabel('Frequency')
% ylabel('Spectrum')
% grid on;
% ax = gca; % current axes
% ax.FontSize = 10;
% legend('Spectrum','peak period','Estimate', 'Estimate-peak period')
% title('Spectrum as a fucntion of angular speed \omega')

% plot of estimated parameters spectrum and true spectrum and their
% corresponding peak period.




% ACF plot of the true spectrum.

% Both Maximum likelihood estimates are most times far from the true parameters.




%% Monte Carlo Markov chain simulation
% ssfun= @(theta, data) -2*loglike(data,theta);
%realisation=mean(X,1)'; %taking mean realisation.
realisation=X(10,:)';
%max_lik_este1=mean(max_lik_este,1); % taking mean estimate.
max_lik_este1=max_lik_este(10,:); % taking mean estimate.

%%
Sppe= pmspec([],max_lik_este1);
R1 = spec2cov(Sppe,1);
[w,x,y]=autocorr(R1.R,'NumLags', 199,'NumStd',2);

% figure
% plotspec(Spec)
% hold on
% plotspec(Sppe,'r')
% grid on;
% legend('peak period','Spectrum', 'Estimate-peak period','Estimate', 'interpreter', 'latex')
% ax = gca; % current axes
% ax.FontSize = 10;
%%
mom = spec2mom(Spec);
levels=[0 1 2]*sqrt(mom(1)); % wav2slope requires absolute levels
Slope0 = wav2slope(L,levels)
tp=dat2tp(realisation);
%%
figure
lc=dat2lc(realisation);
lcplot(lc)
grid on;
xlabel("Level [m]", 'interpreter', 'latex')
ylabel("Crossing intensity", 'interpreter', 'latex')
legend('level-crossing','theoretical Guassian', 'interpreter','latex')
x0=10; y0=10; width=750; height=450;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 12;
%%
figure
plotspec(Spec,'k')
hold on
plotspec(Sppe,'r')
legend('peak period','Spectrum', 'Estimate-peak period','Estimate', 'interpreter', 'latex')
%title('Spectrum as a fucntion of angular speed $\omega$', 'interpreter', 'latex')
x0=10; y0=10; width=550; height=350;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 20;
%%
figure
qqplot(realisation)
xlabel('Standard normal quantiles' , 'interpreter', 'latex')
ylabel('Quantiles of samples' , 'interpreter', 'latex')
title("QQ-plot of sample data", 'interpreter', 'latex')
x0=10; y0=10; width=550; height=350;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 22;

%%




figure
subplot(2,1,1)
plotspec(Sppe)
grid on 
ax = gca; % current axes
ax.FontSize = 10;
% plot of estimated parameters spectrum and true spectrum.

subplot(2,1,2)
plot(x,w)
hold on
yline(y(1))
hold on 
yline(y(2))
grid on
xlabel('Lag', 'interpreter', 'latex')
ylabel('Sample Autocorrelation', 'interpreter', 'latex')
ax = gca; % current axes
ax.FontSize = 10;
%%
ssfun= @(theta, data) -2*loglike(theta,times,realisation); %defining ssfun


% Using toolbox


params= { 
    {'b_0', max_lik_este1(1),0,Inf}
    {'b_1', max_lik_este1(2),0,Inf}
    };

model.ssfun= ssfun; data=realisation; model.N=1;
RR{1}= eye(2)*1.0; RR{2}= eye(2)*0.001;
RR{3}= eye(2)*0.0001; RR{4}= eye(2)*0.0001;
options.RDR= RR; options.ntry=4;
options.nsimu=15000; options.method='dram';
[result,Chain,s2chain] = mcmcrun(model,data,params,options);
%%

chain=Chain;

% %Plot of the chains
figure
subplot(2,1,1)
plot(chain(:,1), 'Color',[0.65 0.65 0.65])
hold on
yline(mean(chain(:,1)),'k','LineWidth',3);
ylabel('$H_{\mathrm{s}}$' , 'interpreter', 'latex')
axis tight
ax = gca; % current axes
ax.FontSize = 16;

subplot(2,1,2)
plot(chain(:,2), 'Color',[0.65 0.65 0.65])
hold on
yline(mean(chain(:,2)),'k','LineWidth',3);
axis tight;
xlabel('Length of chain' , 'interpreter', 'latex')
ylabel('$T_{\mathrm{p}}$' , 'interpreter', 'latex')
x0=10; y0=10; width=750; height=400;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 16;

%% Chain statistics
chain=Chain(5000:end,:);
%chainstats(chain,result)

% [f,x1] = ksdensity(chain(:,1)); 
% figure
% plot(x1,f);

figure
subplot(2,1,1)
[w,x,y]=autocorr(chain(:,1),'NumLags', 399,'NumStd',2);
plot(x,w,'k')
hold on
yline(y(1),'r')
hold on 
yline(y(2),'r')
grid on
ylabel('ACF $H_{\mathrm{s}}$' , 'interpreter', 'latex')
ax = gca; % current axes
ax.FontSize = 16;
subplot(2,1,2)
[w,x,y]=autocorr(chain(:,2),'NumLags', 399,'NumStd',2);
plot(x,w,'k')
hold on
yline(y(1),'r')
hold on 
yline(y(2),'r')
grid on
xlabel('Lag' , 'interpreter', 'latex')
ylabel('ACF $T_{\mathrm{p}}$' , 'interpreter', 'latex')
x0=10; y0=10; width=750; height=400;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 16;

figure
histogram(chain(:,1),'FaceColor',[0 0 0])
xlabel('$H_{\mathrm{s}}$', 'interpreter', 'latex')
hold on
xline(p(1),'r', 'Linewidth', 3)
axis tight;
ax = gca; % current axes
ax.FontSize = 12;
x0=10;
y0=10;
width=350;
height=350;
set(gca,'ytick',[]);
set(gcf,'position',[x0,y0,width,height])

figure
histogram(chain(:,2), 'FaceColor',[0 0 0])
xlabel('$T_{\mathrm{p}}$' , 'interpreter', 'latex')
hold on
xline(p(2),'r', 'Linewidth', 3)
axis tight;
view([90 -90])
ax = gca; % current axes
ax.FontSize = 12;
x0=10;
y0=10;
width=350;
height=350;
set(gca,'ytick',[]);
set(gcf,'position',[x0,y0,width,height])


%%
figure
plot(chain(:,1), chain(:,2),'.', 'color', [0.65 0.65 0.65])
hold on
plot(p(1),p(2), 'or','LineWidth',2, 'MarkerSize', 7)
hold on
plot(mean(chain(:,1)),mean(chain(:,2)),  'ok','LineWidth',2, 'MarkerSize', 7)
xlabel('$H_{\mathrm{s}}$' , 'interpreter', 'latex')
ylabel('$T_{\mathrm{p}}$' , 'interpreter', 'latex')
legend('Chain pairs','True parameter','Estimated parameter' , 'interpreter', 'latex')
axis tight;
x0=10; y0=10; width=350; height=350;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 12;

%%

% 
% figure
% subplot(2,1,1)
% qqplot(chain(:,1))
% subplot(2,1,2)
% qqplot(chain(:,2))




% [f,x1] = ksdensity(chain(:,1)); 
% figure
% plot(x1,f);
%%

figure
plotspec(Spec,'k')
hold on
for i=1:5000
pp=pmspec([],chain(i,:));
plot(pp.w,pp.S,'color','[0.65 0.65 0.65]')
hold on
end
plotspec(Spec,'k')
legend('peak period','Spectrum', 'Estimate-peak period','Estimate', 'interpreter', 'latex')
%title('Spectrum as a fucntion of angular speed $\omega$', 'interpreter', 'latex')
axis tight
x0=10; y0=10; width=550; height=350;
set(gcf,'position',[x0,y0,width,height])
ax = gca; ax.FontSize = 22;
%%
rmpath(genpath('./mcmcstat-master'))
rmpath(genpath('mcmcstat-master'))
%%
pd=makedist('rayleigh', 'b', 1.68037)
pd2=makedist('Weibull', 'a', 14.2, 'b', 4.2)
pd3=makedist('GeneralizedExtremeValue', 'k', -0.271918, 'sigma', 2.3948, 'mu',-0.849724)


% figure
% hist(realisation)
figure
subplot(2,2,1)
h=qqplot(realisation,pd)
set(h(1),'markeredgecolor',[0 0 0]);
xlabel('Quantiles of rayleigh distribution' , 'interpreter', 'latex')
ylabel('Quantiles of input samples' , 'interpreter', 'latex')
title("")
grid on;
axis square
ax = gca; % current axes
ax.FontSize = 8;
subplot(2,2,2)
h=qqplot(realisation,pd2)
set(h(1),'markeredgecolor',[0 0 0]);
xlabel('Quantiles of weibull distribution' , 'interpreter', 'latex')
ylabel('Quantiles of input samples' , 'interpreter', 'latex')
title("")
grid on;
axis square
ax = gca; % current axes
ax.FontSize = 8;
subplot(2,2,3)
h=qqplot(realisation, pd3)
set(h(1),'markeredgecolor',[0 0 0]);
xlabel('Quantiles of GEV distribution' , 'interpreter', 'latex')
ylabel('Quantiles of input samples' , 'interpreter', 'latex')
title("")
grid on;
axis square
ax = gca; % current axes
ax.FontSize = 8;
subplot(2,2,4)
h=qqplot(realisation)
set(h(1),'markeredgecolor',[0 0 0]);
xlabel('Standard normal quantiles' , 'interpreter', 'latex')
ylabel('Quantiles of input samples' , 'interpreter', 'latex')
title("")
grid on;
axis square
x0=10;
y0=10;
width=650;
height=450;
set(gcf,'position',[x0,y0,width,height])
ax = gca; % current axes
ax.FontSize = 8;
%%
% Chain statistics
% global R;
% global IR;
% R(1,:,:)= eye(2)*0.1;
% IR(1,:,:)= inv(squeeze(R(1,:,:)));
% 
% for i=2:3
%     R(i,:,:)= R(i-1,:,:)*sqrt(0.5);
%     IR(i,:,:)= squeeze(IR(i-1,:,:))*sqrt(1/0.5);
% end

N=2500;
t=cputime;
dramc=chain;
dramacf= autocorr(dramc(:,1),'Numlags',20);
EssDram= N/(1+2*sum(dramacf(2:end)))
OESDram= EssDram/t



%% Functions

% Likelihood function.
function ll = loglike(p,times,realisation)
    if (any(p<0))
        ll = -Inf;
        return
    end

    Spec = pmspec([],[p(1) p(2)]); Spec.h = Inf;
    R1 = spec2cov(Spec,1);
         n = numel(times);
    differ = times(2)-times(1);
    
    k = zeros(n,1);
    % Loop over the time points of the realisation and interpotalte the
    % covariance values from the vectors (time points R1.t and covaraince values R1.R) given by the toolbox.
     for i=1:n
         q = interp1(R1.t,R1.R,(i-1)*differ,'linear','extrap');
         k(i) = q;
     end
     % Construct the covariance matrix s for the wave process
     % and add small diagonal matrix to it (i.e add white noise process).
     % Otherwise the matrix would be singular.
     ccov=toeplitz(k,k');
     s = eye(n)*0.3;
     ccov = ccov+s; ccov = 0.5*(ccov' +ccov); 
     iCov = (ccov)\eye(n);
     
     R = cholcov(ccov);

     % Log-determinant of the full covariance.
     logdete = 2*sum(log(diag(R)));
     
    
 % Overall log-likelihood for the measurement.
 ll = -0.5*realisation'*iCov*realisation - 0.5*(n*log((2*pi))+logdete);

end

% spectrum function

function spec = S(w,p)
    if (w<=eps)
        spec= 0;
    else
        spec = (5*(p(1))^2./(p(2).*(w./p(2)).^(5))).*exp((-5/4).*(w./p(2)).^(-4));
    end
end

