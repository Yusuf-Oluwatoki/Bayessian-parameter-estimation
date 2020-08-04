close all
clear all
clc
addpath(genpath('C:\Users\user\Dropbox\Thesis\wafo_2017\wafo'))
addpath(genpath('C:\Users\user\Dropbox\mcmcstat-master'))

%% Pierson Moskowitz spectrum

S = pmspec([],[6.5 10]); S.h = 8;
figure
plotspec(S,'k')
title('Spectral density', 'interpreter','latex')
grid on
x0=10;
y0=10;
width=600;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 12;
%% Realisation from the spectrum using the tool box
figure

rng default
opt = simoptset('Nt',2256,'dt',0.15,'lalpha', 0.8, 'ffttype', 'ffttime');
[w,x] = spec2ldat(S,opt) ; % Keep [w,x]
horizon=x.u+x.Z(:,128);
verti=w.Z(:,128);
plot(horizon,verti)
axis tight
% axis([0 500 -10 10]) % Keep the figure
xlabel('Horizontal component')
ylabel('Vertical component')
title('simulated space waves for P–M orbital spectrum at h=8, \alpha = 0.8')
hold on

%%
figure
verti(:,1)=1:2048;
verti(:,2)=w.Z(:,128);
NIT = 3, paramt = [0 10 51];
dtyex = spec2tpdf(S,[],'Tt',paramt,0,NIT);
[T, index] = ldat2lwav(verti);
histgrm(T,25,1,1), hold on
pdfplot(dtyex)
axis([0 10 0 0.35]), hold off
%%
ww=linspace(0,50,256);

[H,Tp]=pierson(w,ww);



%Using toolbox

params= { 
    {'b_0', H,0,Inf}
    {'b_1', Tp,0,Inf}
    };

model.modelfun= modelfun; data=ww; model.N=1;
RR{1}= eye(2)*1.0; RR{2}= eye(2)*0.001;
RR{3}= eye(2)*0.0001; RR{4}= eye(2)*0.0001;
options.RDR= RR; options.ntry=4;
options.nsimu=100; options.method='dram';
[result,chain,s2chain] = mcmcrun(model,data,params,options);



figure
subplot(2,1,1)
plot(chain(:,1))
subplot(2,1,2)
plot(chain(:,2))

figure
for i=1:7
rng default
S = pmspec(1.5,[chain(i,1) chain(i,2)]); S.h = 8;
opt = simoptset('Nt',256,'dt',0.125,'Nu',256*8,'du',0.25,'lalpha', 0.8);
[w,x] = spec2ldat(S,opt) ; % Keep [w,x]
horizon(:,i)=x.u+x.Z(:,128);
verti(:,i)=w.Z(:,128);
plot(S.S)
hold on
title('P-M spectral density for \omega =1.5 and different parameters')
end


%%
figure
for i=1:7

plot(horizon, verti)
axis([0 500 -10 10]) % Keep the figure
xlabel('Horizontal component')
ylabel('Vertical component')
title('simulated space waves for P–M orbital spectrum at h=8, \alpha = 0.8')
hold on
end
%%
global R;
global IR;
R(1,:,:)= eye(2)*0.1;
IR(1,:,:)= inv(squeeze(R(1,:,:)));

for i=2:3
    R(i,:,:)= R(i-1,:,:)*sqrt(0.5);
    IR(i,:,:)= squeeze(IR(i-1,:,:))*sqrt(1/0.5);
end

N=2500;
t=cputime;
dramc=chain;
dramacf= autocorr(dramc(:,1),'Numlags',20);
EssDram= N/(1+2*sum(dramacf(2:end)))
OESDram= EssDram/t


    
function [H, Tp]= pierson(w)
H= 4*sqrt(w.mom(1));
Tp=2*pi*w.mom(2);
modelfun=@(data,theta)(5*(H)^2/(Tp.*(data./Tp).^(5))).*exp((-5/4).*(data./Tp).^(-4));
end