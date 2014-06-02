function [tpeak_env_coda,tmax,coda_env,coda_env_smoothed,noise_level,norder]=calc_coda_envelope(S,S2,f1,f2)
% tpeak_env_coda: time of the peak in the smoothed coda envelope
% tmax: is the time with S/N > 3 for t > tpeak_env_coda
%       if no S/N > 3 exist, tmax=S.DELTA
addpath(genpath('~/octave'),genpath('~/prog/octave'));

time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
% Demean
data=detrend(S.DATA1,0);
data2=detrend(S2.DATA1,0);
% Detrend
data=detrend(data,1);
data2=detrend(data2,1);
% taper
data=data.*tukeywin(length(data),0.05);
data2=data2.*tukeywin(length(data2),0.05);

% Reference data from origin time S.O, correct for S.B if necessary
% Test if picked arrival times exist (S.A and S.T0) otherwise use theoretical times S.T1 and S.T2 for P and S
NO=round((S.O-S.B)/S.DELTA);
if ~isnan(S.A), NP=round((S.A-S.B)/S.DELTA);, elseif ~isnan(S.T1), NP=round((S.T1-S.B)/S.DELTA);, else, disp('P arrival not defined. Exit.'), return, end

Fnyquist=1/(2*S.DELTA);

coda_env_smoothed=[];
coda_env=[];
tmax=[];
tpeak_env_coda=[];
noise_level=[];
norder=[];
for i=1:length(f1)
   w1=f1(i)/Fnyquist;
   w2=f2(i)/Fnyquist;
   % Filter instabilities depending on order and frequencies (for lower frequencies, instabilities
   % appear for low order (2-4) while for higher frequencies for order (4-8)
   % Test using the gain output from butter, apparently when it goes below 1e-10 instabilities
   % appear.
   % Goal order 4, but in case gain is <1e-10, decrease order.
   norder(i)=4;
   [z,p,k] = butter(norder(i),[w1 w2]);
   while k < 1e-9
      norder(i)=norder(i)-1;
      if norder(i) < 2, disp('problem with the filter order, exit'), return, end
      [z,p,k] = butter(norder(i),[w1 w2]);
   end
   [B,A]=zp2tf(z,p,k);
   
   X=hilbert(filter(B,A,data));
   X2=hilbert(filter(B,A,data2));
   X=log10(abs(X));
   X2=log10(abs(X2));
   X_average=(X+X2)./2;
   
   noise_level(i)=mean(X_average(NO:NP),1);
   data_smoothed=smooth(X_average,100);
   Npeak_env_coda=find(data_smoothed-max(data_smoothed) >= 0);
   tpeak_env_coda(i)=Npeak_env_coda*S.DELTA;
   % Tests on the level of the signal-to-noise ratio
   % If toto equals Npeak_env_coda, it means that the peak of the coda envelope
   % is in the level of noise, no good data
   % If toto equals 0, the whole length is good
   toto=0;
   for j=Npeak_env_coda:length(X_average)
       if data_smoothed(j) < noise_level(i)+log10(3.0), toto=j;, break, end
   end
   if toto == Npeak_env_coda, tmax(i)=S.DELTA;, elseif toto == 0, tmax(i)=length(X_average)*S.DELTA, else, tmax(i)=j*S.DELTA;, end
   coda_env=[coda_env X_average];
   coda_env_smoothed=[coda_env_smoothed data_smoothed];
   
end

endfunction