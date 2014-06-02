
clear all
close all
addpath(genpath('~/octave'),genpath('~/prog/octave'));

path_name='/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3/2012/2012_12_19_04_54_37_700/';
% List all the files for that directory
tmp=dir([path_name "*HHE*.sac"]);
% Search the maximum sampling rate in order to define the maximum common frequency available
max_sampling_rate=0;
for i=1:size(tmp,1)
   S=readsac([path_name tmp(i).name]);
   if S.DELTA > max_sampling_rate, max_sampling_rate=S.DELTA;, end
end
f_max_common=1/(2*max_sampling_rate);
% Define the frequency bands for the coda analysis
% f is the central frequency. f1=f/2 and f2=f*2
% start with f=0.1 Hz, stop when f2 gets larger than f_max_common/1.25
f=0.1/2;
i=0;
while f < f_max_common/(2.0*1.25)
   i=i+1;
   f=f*2;
   f1(i)=f/sqrt(2.0);
   f2(i)=f*sqrt(2.0);
%     fprintf(1,'%s %f %s %f %s %f\n','fmin=',f1,' - f=',f,' - fmax=',f2)
end
% Compute the coda enveloppes for each station
for i=1:size(tmp,1)
   clear f_nameH1 f_nameH2 S S2
   f_nameH1=tmp(i).name;
   f_nameH2=strrep(f_nameH1,'HHE','HHN');
   S=readsac([path_name f_nameH1]);
   S2=readsac([path_name f_nameH2]);
   clear tpeak_env_coda tmax coda_env coda_env_smoothed noise_level norder
   [tpeak_env_coda,tmax,coda_env,coda_env_smoothed,noise_level,norder]=calc_coda_envelope(S,S2,f1,f2);
   plot_coda_envelopes(S,S2,f1,f2,tpeak_env_coda,tmax,coda_env,coda_env_smoothed,noise_level,norder)
   system(['cp filtered_time_series.jpg ' f_nameH1 '_filtered_ts.jpg'])
   system(['cp coda_envelopes.jpg ' f_nameH1 '_coda_env.jpg'])
end
