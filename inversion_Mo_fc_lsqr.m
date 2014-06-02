function [Mw,fc,std_Mw_fc,mod,kvg1]=inversion_Mo_fc_lsqr(freq,log10_acc_obs,sig_log10_acc_obs,Mw_start);
% Define parameters and starting point for inversion
rad=0.55;
vs=3500;
rho=2800;
%  Mw_start=3.0;
fc_start=10^(2-0.5*Mw_start);

%  A=load('data_corrected.in');
%  %  freq=A(:,1);
%  %  log10_acc_obs=log10(A(:,2));
%  %  sig_log10_acc_obs=A(:,3);
%  freq=A(find(A(:,4)==1),1);
%  log10_acc_obs=log10(A(find(A(:,4)==1),2));
%  sig_log10_acc_obs=A(find(A(:,4)==1),3);

% Here x is frequency
leasqrfunc = @(x, p) p(1)+log10(((2*pi*x).^2)./(1+(x./10^(p(2))).^2));
leasqrdfdp = @(x, f, p, dp, func) [1.0*ones(size(x)), (2*x.^2)./(10^(2*p(2))+x.^2)];
% Set params(1)=log10(M0)+log10(cste)
p(1)=1.5*Mw_start+9.1+log10(2*rad/(4*pi*rho*vs^3));
% Set params(2)=log10(fc)
p(2)=log10(fc_start);
p=p';

%  log10_acc_obs=leasqrfunc(freq,p);
[freq,i]=sort(freq);
log10_acc_obs=log10_acc_obs(i);

figure(1)
semilogx(freq,log10_acc_obs,'.')

F = leasqrfunc;
dFdp = leasqrdfdp; % exact derivative
dp = [0.0001; 0.0001];
pin = [0.0; 0.0]; 
stol=0.001; niter=100;

%  figure(2)
%  global verbose;
%  verbose = 1;

[f1, p1, kvg1, iter1, corp1, covp1, covr1, stdresid1, Z1, r21] = ...
    leasqr (freq, log10_acc_obs, pin, F, stol, niter, sig_log10_acc_obs, dp, dFdp);

Mw=(p1(1)-log10(2*rad/(4*pi*rho*vs^3))-9.1)/1.5;
fc=10^(p1(2));
std_Mw_fc=sqrt(diag(covp1));
std_Mw_fc(1)=std_Mw_fc(1)/1.5;
std_Mw_fc(2)=std_Mw_fc(2)*fc;

mod=10.^leasqrfunc(freq,p1);
figure(1)
hold on
semilogx(freq,leasqrfunc(freq,p),'k')
semilogx(freq,leasqrfunc(freq,p1),'r')
legend('data','Starting model','Inverted model')
hold off

%  fid_out1=fopen('Mw_fc.out','w');
%  fprintf(fid_out1,'%s %3.1f\n','Mw=',Mw);
%  fprintf(fid_out1,'%s %6.2f\n','fc=',fc);
%  fclose(fid_out1);
%  fid_out2=fopen('model_source_acc.out','w');
%  fprintf(fid_out2,'%f %e \n',[freq 10.^leasqrfunc(freq,p1)]');
%  fclose(fid_out2);
%  %  fid_out3=fopen('rmse.out','w');
%  %  fprintf(fid_out3,'%e \n',rmse);
%  %  fclose(fid_out3);
