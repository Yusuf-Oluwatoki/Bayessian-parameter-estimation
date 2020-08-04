close all
clear all 
clc
%% Adding the toolbox paths.

addpath(genpath('../wafo'))
addpath(genpath('./mcmcstat-master'))

%% Setting parameters for the spectrum realisation. 

N = 5500; Nspec = 256;
% Number of points for a realisation and for the spectrum function.

p = [10,12]; 
% Significant wave height and peak period of the spectrum.
%omega=linspace(0,3.2,Nspec);
% Range for spectrum's angular speed.
% rng('default')
Spec = pmspec([],p); Spec.h = 8;
opt = simoptset('Nt',N,'Nu', 20, 'lalpha', 0.8,'ffttype','ffttime');
% Spectrum realisation and setting options for wave realisation.

R1 = spec2cov(Spec,1); %spectrum to covariance from toolbox.

% Get a realisation from the toolbox.
[realisation2,realisedX] = spec2ldat(Spec,opt);

% realisation time steps.
times = realisation2.t;
 
figure
plot(realisation2.Z(1,:))
axis tight
grid on;
% title('Wave realisations from covariance matrix')
xlabel('$t$', 'interpreter','latex')
ylabel('$W(0,t)$', 'interpreter','latex')
ax = gca; % current axes
ax.FontSize = 10;


% plot of first realisation out of a Nu realisations.
%% Level crossing
figure
lc=dat2lc(realisation2.Z(1,:));
lcplot(lc)
grid on;
ax = gca; % current axes
ax.FontSize = 10;
legend('level-crossing','theoretical Guassian')
%%
% [L2,L02] = ldat2lwav(realisation2,realisedX,'time');
% mom = spec2mom(Spec);
% levels=2*sqrt(mom(1));
% Slope2 = wav2slope(L02,levels);
% 
% figure
% plotedf(Slope2.up{1}); hold on
% plotedf(-Slope2.down{1},'r-.');

% %% Error bars of ten parameters
% N = 248; Nspec = 256;
% 
% H=[6, 8, 7, 9, 10, 5]; T=[7, 8, 12, 9, 10, 13];
%  for i=1:length(H)
%      p=[H(i),T(i)];
%      Spec=pmspec([],p); Spec.h = 8;
%      opt = simoptset('Nt',N,'Nu', 20, 'lalpha', 0.8,'ffttype','ffttime');
%      [realisation2,realisedX] = spec2ldat(Spec,opt,'iseed', rng('default'));
%         times = realisation2.t;
% 
% for j=1:5
% max_lik_este(j,:) = fminsearch(@(p) -loglike(p,times,realisation2.Z(j,:)'),p);
% end
% max_lik_este3(i,1)=sqrt(mean((H(i)-max_lik_este(:,1)').^2))
% max_lik_este3(i,2)=sqrt(mean((T(i)-max_lik_este(:,2)').^2))
% 
%  end
%  err = max_lik_este3(:,1)';
%  err2 = max_lik_este3(:,2)';
%  %%
% figure
% subplot(2,1,1)
% errorbar(1:6,H,err)
% axis([0 7 2 12.5])
% ylabel('\theta_1')
% grid on;
% ax = gca; % current axes
% ax.FontSize = 10;
% subplot(2,1,2)
% errorbar(1:6,T,err2)
% axis([0 7 5 16])
% xlabel('Number of parameters')
% ylabel('\theta_2')
% grid on;
% ax = gca; % current axes
% ax.FontSize = 10;
%% Maximum likelihood estimation of parameters.


for i=1:1
max_lik_este2(i,:) = fminsearch(@(p) -loglike(p,times,realisation2.Z(i,:)'),[10,12])
end
% Estimating the parameters of the first ten wave realisations.


figure
plot(Spec.S)
hold on
for i=1:10
 Sppe= pmspec([],max_lik_este2(i,:));
plot(Sppe.S)
hold on
end
xlabel('Frequency')
ylabel('Spectrum')
grid on;
ax = gca; % current axes
ax.FontSize = 10;
title('Spectrum as a fucntion of angular speed \omega')
% plot of estimated parameters spectrum and true spectrum.


figure
% plot(R1.t,R1.R)
title('Covariance function as a function of lag \tau')
[w,x,y]=autocorr(R1.R,'NumLags', 199,'NumStd',2);
plot(x,w)
hold on
yline(y(1))
hold on 
yline(y(2))
grid on
xlabel('Lag')
ylabel('Sample Autocorrelation')
ax = gca; % current axes
ax.FontSize = 10;

% ACF plot of the true spectrum.

% First parameter is mostly diverging from its true value. Second parameter of the maximum likelihood 
% estimate is most times greater than the correct one.

%% Monte Carlo Markov chain simulation
% ssfun= @(theta, data) -2*loglike(data,theta);
realisation=mean(realisation2.Z,1)'; %taking mean realisation.
max_lik_este=mean(max_lik_este2,1); % taking mean estimate.
ssfun= @(theta, data) -2*loglike(theta,times,realisation); %defining ssfun


% Using toolbox


params= { 
    {'b_0', max_lik_este(1),0,Inf}
    {'b_1', max_lik_este(2),0,Inf}
    };

model.ssfun= ssfun; data=realisation; model.N=1;
RR{1}= eye(2)*0.0001; RR{2}= eye(2)*0.001;
RR{3}= eye(2)*0.01; RR{4}= eye(2)*0.1;
% RR{1}= eye(2)*1.0; RR{2}= eye(2)*0.001;
% RR{3}= eye(2)*0.0001; RR{4}= eye(2)*0.0001;
options.RDR= RR; options.ntry=4;
options.nsimu=10000; options.method='dram';
[result,chain,s2chain] = mcmcrun(model,data,params,options);


% %Plot of the chains
figure
subplot(2,1,1)
plot(chain(:,1))
ylabel('\theta_1')
grid on;
ax = gca; % current axes
ax.FontSize = 10;
subplot(2,1,2)
plot(chain(:,2))
axis tight;
grid on;
%title('(b)')
xlabel('Length of chain')
ylabel('\theta_2')
ax = gca; % current axes
ax.FontSize = 10;

%%

figure
subplot(2,1,1)
[w,x,y]=autocorr(chain(:,1),'NumLags', 399,'NumStd',2);
plot(x,w)
hold on
yline(y(1))
hold on 
yline(y(2))
grid on
ylabel('Sample Autocorrelation')
ax = gca; % current axes
ax.FontSize = 10;
subplot(2,1,2)
[w,x,y]=autocorr(chain(:,2),'NumLags', 399,'NumStd',2);
plot(x,w)
hold on
yline(y(1))
hold on 
yline(y(2))
grid on
xlabel('Lag')
ylabel('Sample Autocorrelation')
ax = gca; % current axes
ax.FontSize = 10;
%% Chain statistics

chainstats(chain,result)

% [f,x1] = ksdensity(chain(:,1)); 
% figure
% plot(x1,f);

figure
subplot(2,2,1)
histogram(chain(:,1))
xlabel('$H_{\mathrm{s}}$', 'interpreter', 'latex')
grid on;
axis square
ax = gca; % current axes
ax.FontSize = 10;
subplot(2,2,4)
histogram(chain(:,2))
xlabel('$T_{\mathrm{p}}$' , 'interpreter', 'latex')
grid on;
view([90 -90])
axis square
ax = gca; % current axes
ax.FontSize = 10;
subplot(2,2,3)
mcmcplot(chain, 1:2, [],'pairs')
xlabel('$H_{\mathrm{s}}$' , 'interpreter', 'latex')
ylabel('$T_{\mathrm{p}}$' , 'interpreter', 'latex')
grid on;
axis square
ax = gca; % current axes
ax.FontSize = 10;


figure
subplot(2,1,1)
qqplot(chain(:,1))
subplot(2,1,2)
qqplot(chain(:,2))
%%
rmpath(genpath('./mcmcstat-master'))
rmpath(genpath('mcmcstat-master'))
pd=makedist('rayleigh', 'b', 1.68037)
pd2=makedist('Weibull', 'a', 14.2, 'b', 4.2)
pd3=makedist('GeneralizedExtremeValue', 'k', -0.271918, 'sigma', 2.3948, 'mu',-0.849724)


% figure
% hist(realisation)
figure
subplot(2,2,1)
qqplot(realisation,pd)
subplot(2,2,2)
qqplot(realisation,pd2)
subplot(2,2,3)
qqplot(realisation, pd3)
subplot(2,2,4)
qqplot(realisation)

%%
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

