clear all
clc
addpath('C:\Users\user\Dropbox\Thesis\wafo_2017\wafo')
%%
% S = pmspec(1.5,[6.5 10]); 
% figure
% plotspec(S)

S = pmspec(2.5,[6.5 10]); S.h = 8;
figure
plotspec(S,'k')
%title('P-M spectral density','interpreter','latex')
%grid on
x0=10;
y0=10;
width=400;
height=200;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 8;

%%
rng default
opt = simoptset('Nt',256,'dt',0.125,'Nu',256*8,'du',0.25,'lalpha', 0);
[w,x] = spec2ldat(S,opt) ; % Keep [w,x]
figure
subplot(2,1,1)
plot(x.u+x.Z(:,128),w.Z(:,128),'k','Linewidth',2);
axis([0 500 -10 10]) % Keep the figure
% xlabel('Horizontal component')
ylabel('Vertical [m] ','interpreter','latex')
% title('simulated space waves for P–M orbital spectrum at h=8, \alpha = 0')
%grid on;
%title('(a)')
ax = gca; % current axes
ax.FontSize = 8;

opt = simoptset('Nt',256,'dt',0.125,'Nu',256*8,'du',0.25, 'lalpha', 0.8);
[w,x] = spec2ldat(S,opt) ; % Keep [w,x]
subplot(2,1,2)
plot(x.u+x.Z(:,128),w.Z(:,128),'k','Linewidth',2);
axis([0 500 -10 10]) % Keep the figure
xlabel('Horizontal component [m]','interpreter','latex')
ylabel('Vertical [m]','interpreter','latex')
% axis tight;
%grid on;
%title('(b)')
%grid on
x0=10;
y0=10;
width=500;
height=150;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 8;
%
S = pmspec(1.5,[6.5 10]); S.h = inf;
opt = simoptset('Nt',256,'dt',0.125,'Nu',256*8,'du',0.25, 'lalpha', 0);
[w,x] = spec2ldat(S,opt) ; % Keep [w,x]
figure
subplot(2,1,1)
plot(x.u+x.Z(:,128),w.Z(:,128),'k','Linewidth',2);
axis([0 500 -10 10]) % Keep the figure
ylabel('Vertical [m]','interpreter','latex')
% title('simulated space waves for P–M orbital spectrum at h=8, \alpha = 0')
%grid on;
%title('(a)')
ax = gca; % current axes
ax.FontSize = 8;

opt = simoptset('Nt',256,'dt',0.125,'Nu',256*8,'du',0.25, 'lalpha', 0.8);
[w,x] = spec2ldat(S,opt) ; % Keep [w,x]

subplot(2,1,2)
plot(x.u+x.Z(:,1),w.Z(:,1),'k','Linewidth',2);
axis([0 500 -10 10]) % Keep the figure
xlabel('Horizontal component [m]','interpreter','latex')
ylabel('Vertical [m]','interpreter','latex')
%axis tight;
%grid on;
%title('(b)')
%grid on
x0=10;
y0=10;
width=500;
height=150;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 8;

%%
% Lmax = 500;
% SS = dat2spec(x.Z,Lmax);
% plotspec(SS); axis([0 5 0 0.7])
ww(:,1)=1:2048;
ww(:,2)=w.Z(:,1);
Lmax = 200;
R1 = spec2cov(S,1);
Rest = dat2cov(ww,Lmax);
figure
covplot(R1,Lmax,[],'-'), hold off

%%
% % XX=x.u+x.Z(:,128); YY=w.Z(:,128);
% 
% 
% 
% %%
%%
% clear all;
% close all;
% %%ccov=cov(w.Z);
% n=500;
% X=zeros(1,n);
% X(1)= randn()*10;
% lambda=0.067; sigma2=0.015;
% 
% for j=2:n
%     X(j)=X(j-1)*exp(-lambda)+ sqrt(sigma2/ (2*lambda)*(1-exp(-2*lambda)))*randn();
% end
% Lmax = 500;
% xx(:,1)=1:500;
% xx(:,2)=X;
% me=mean(xx(:,2));
% sd= std(xx(:,2));
% xx(:,2)=xx(:,2)-me;
% lc = dat2lc(xx);
% 
% figure
% plotflag = 2;
% lcplot(lc,plotflag,0,sd)
% 
% SS = dat2spec(xx,Lmax);
% figure
% plotspec(SS); axis([0 5 0 0.7])
% 
% %%
% Lmax1 = 200;
% Lmax2 = 50;
% SS1 = dat2spec(xx,Lmax1);
% SS2 = dat2spec(xx,Lmax2);
% figure
% plotspec(SS1,[],'-.'), hold on
% plotspec(SS2), hold off
% 
% 
% 
% % Lmax = 80;
% R1 = spec2cov(SS1,1);
% Rest = dat2cov(xx,Lmax1);
% figure
% covplot(R1,Lmax1,[],'.'), hold on
% covplot(Rest), hold off