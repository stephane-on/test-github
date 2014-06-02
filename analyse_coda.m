
clear all
close all
addpath(genpath('~/octave'),genpath('~/prog/octave'));
disp('Dummy message')

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
=======
clear all
close all
addpath(genpath('~/octave'),genpath('~/prog/octave'));
disp('Needs some octave-forge extra packages');
disp('Modify again the file in branch')
disp('c est vraiment l enfer git')

% Loop to cycle through 2012 folder
dirlist=glob("/home/michel/Earthquakes_sac_SC3/2012/2012_*/");
for j=1:length(dirlist)
   path_name=dirlist{j,1};
   % path_name='/home/michel/Spectral_analysis/Prog_octave/2012_03_31_16_09_51_100/';
   % List all the files for that directory
   tmp=glob([path_name "*HHE*.sac"]);
   % Search the maximum sampling rate in order to define the maximum common frequency available
   max_sampling_rate=0;
   for i=1:length(tmp)
      S=readsac(tmp{i,1});
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
   end
end
% Creating Velocity for frequencies file
pname=('/home/michel/Earthquakes_sac_SC3/2012/');
fid_out=fopen([pname 'Vel_f.txt'],'w');
fprintf(fid_out,'%s\n','file_name                         Distance     Vel_f1      Vel_f2      Vel_f3      Vel_f4      Vel_f5      Vel_f6      Vel_f7      Vel_f8      Vel_f9');
fid_outd=fopen([pname 'Dist_f.txt'],'w');
fprintf(fid_outd,'%s\n','Distance');
fid_outv=fopen([pname 'Vel_f6.txt'],'w');
fprintf(fid_outv,'%s\n','Vel_f for frequency 6');

% Loop to write the data
for j=1:length(dirlist)
   path_name=dirlist{j,1};
   tmp=glob([path_name "*HHE*.sac"]);
   % Compute the coda enveloppes for each station
   for i=1:length(tmp)
      clear f_nameH1 f_nameH2 S S2 
      f_nameH1=tmp{i,1};
      f_nameH2=strrep(f_nameH1,'HHE','HHN');
      S=readsac(f_nameH1);
      S2=readsac(f_nameH2);
      clear tpeak_env_coda tmax coda_env coda_env_smoothed noise_level norder
      [tpeak_env_coda,tmax,coda_env,coda_env_smoothed,noise_level,norder]=calc_coda_envelope(S,S2,f1,f2);
      plot_coda_envelopes(S,S2,f1,f2,tpeak_env_coda,tmax,coda_env,coda_env_smoothed,noise_level,norder)
      system(['cp filtered_time_series.jpg ' f_nameH1 '_filtered_ts.jpg']);
      system(['cp coda_envelopes.jpg ' f_nameH1 '_coda_env.jpg']);
   % calculating Vel_f
      clear Vel_f G
      Vel_f=(S.DIST./tpeak_env_coda);
      G=find((tmax-(S.DELTA*ones(length(f1),1))')==0);
      Vel_f(G)=0*ones(length(G),1);
   % printing out
      formato=[];
      for g=1:length(f1)
         formato=[formato, ' %10.6f '];
      end
      fprintf(fid_out,['%s %10.5f ',formato,'\n'],substr(tmp{i,1}, 39, 32),S.DIST,Vel_f);
      fprintf(fid_outd,'%f/n ',S.DIST);
      fprintf(fid_outv,'%f/n ',Vel_f(6));
   end
end
fclose(fid_out);
fclose(fid_outd);
fclose(fid_outv);
