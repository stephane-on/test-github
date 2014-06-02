%  %function model_coda(amp_env_coda,v,b,gamma,time_from_Origin,dist)
%  
%  time_from_Origin=[0:0.01:0.01*58442];
%  dist=324.09;
addpath(genpath('~/octave'),genpath('~/prog/octave'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    TEST
ifreq=5;
time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
dist=S.DIST;
max_cod_env_smoothed=max(coda_env_smoothed(:,ifreq));

time_for_model=[S.DELTA:S.DELTA:max(time)-tpeak_env_coda(ifreq)];
v=[3.45 96 35];
b=[-0.0075 0.10 500];
gamma=[0.1 -10 6];
v_r=v(1)-v(2)/(v(3)+dist);
b_r=b(1)-b(2)/(b(3)+dist);
gamma_r=gamma(1)-gamma(2)/(gamma(3)+dist);
model=exp(b_r*time_for_model).*time_for_model.^-gamma_r;

plot(time.-S.O,coda_env_smoothed(:,ifreq))
hold on
gamma_r=0;
for i=1:10
   gamma_r=i*0.02;
   model=exp(b_r*time_for_model).*time_for_model.^-gamma_r;
   plot(time_for_model+tpeak_env_coda(ifreq)-S.O,model-max(model)+max_cod_env_smoothed,'r')
end
%  plot(time_from_Origin,model.-max(model)+max_cod_env_smoothed,'r')
gamma_r=gamma(1)-gamma(2)/(gamma(3)+dist);
b_r=-0.1;
model=exp(b_r*time_for_model).*time_for_model.^-gamma_r;
   plot(time_for_model+tpeak_env_coda(ifreq)-S.O,model-max(model)+max_cod_env_smoothed,'g')
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the whole coda
plot(time,coda_env_smoothed(:,ifreq))
hold on
% plot the part usable for computing coda parameters
plot(time(tpeak_env_coda(ifreq)/S.DELTA:tmax(ifreq)/S.DELTA),coda_env_smoothed(tpeak_env_coda(ifreq)/S.DELTA:tmax(ifreq)/S.DELTA,ifreq),'r')
hold off
coda_for_fit=coda_env_smoothed(tpeak_env_coda(ifreq)/S.DELTA:tmax(ifreq)/S.DELTA,ifreq);
time_for_fit=[0.01:S.DELTA:(length(coda_for_fit)+1)*S.DELTA]';
%  coda_normalised=coda_for_fit./max(coda_for_fit);
%  gamma_r=0;
for i=1:10
   gamma_r=i*0.1
   [coef_fit, struct_fit] = polyfit(time_for_fit,log10(10.^(coda_for_fit).*time_for_fit.^gamma_r),1);
   struct_fit.normr
   plot(time_for_fit,log10(10.^(coda_for_fit).*time_for_fit.^gamma_r))
   hold on
   plot(time_for_fit,struct_fit.yf,'r')
end
hold off

