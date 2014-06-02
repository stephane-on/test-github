function manual_t_picks
clear
format long
close all
addpath(genpath('~/octave'));
% Read filenames from list_sac_files.txt
% Define time windows for noise, P-, and S-waves automatically and allows repicking if windows are not fine.
% Results are written in data_dir/spectra/list_picks.txt. Backup file is created if list_picks.txt already exists.

[fname,pname]=uigetfile('list_*.txt','Select a list of filenames for one event','/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3');
fid_duration=fopen(strcat(pname,'distance_duration.txt'),'w');
if exist(strcat(pname,'spectra'),'dir') == 0
    mkdir(pname,'spectra')
end

fich_t_picks=[pname 'spectra/list_picks.txt'];
if exist(fich_t_picks,'file') == 2
    backup_file_pick(fich_t_picks);
    fid_picks_t=fopen(fich_t_picks,'r');
    info_peaks_t=textscan(fid_picks_t,'%s %s %s %f %f %f %f %f %f');
    fclose(fid_picks_t);
    figure(1)
    clf
    plot_all_time_series(pname,info_peaks_t)    
end
fich_t_picks_new=[pname 'spectra/list_picks.txt_new'];
fid_picks_t_new=fopen(fich_t_picks_new,'w');

% Data are assumed to be in m/s
conv=1
answer=inputdlg('Conversion factor from data unit to m/s (if data are in m/s, conv=1, if data are in nm/s, conv=1e-9).','Enter unit converion factor',1,{'1'});
conv=str2double(answer{1});

list_sacfiles=read_file_list_string([pname fname]);
for i=1:size(list_sacfiles,1)
    sacfile=deblank(list_sacfiles(i,:));
    fprintf(1,'%s\n',[pname sacfile]);
    %% Read SAC data
    S=readsac([pname sacfile]);
    %% Test if components is recognised from S.KCMPNM or from filename
    if exist('comp') == 1
        clear comp
    end
    comp=get_comp_name(deblank(S.KCMPNM),sacfile);
    if comp == -99
       disp('component not recognized'),disp([pname sacfile])
       % If component not recognized, jump to next file
       continue
    end
    % Test is necessary data are not empty: S exists means reading sac is ok, and header variables
    % S.DIST and S.DELTA are defined
    if (~isempty(S) && isnan(S.A) ~=1 && isnan(S.T0) ~= 1)
        clear time data_tmp data
        time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
        %% demean, detrend
        data_tmp=detrend(S.DATA1*conv,'constant');
        data=detrend(data_tmp,'linear');
        clear data_tmp
        data=data.*tukeywin(length(data),0.05);
        clear Fnyquist w1 w2 A B data_filter1
        Fnyquist=1/(2*S.DELTA);
        w1=1.0/Fnyquist;
        % WARNING USING w2 larger i.e. 0.999999999999 induce some strange behaviour of the signal
        % the amplitude seems to exponentially increase for large times
        w2=0.9999;
        [B,A] = butter(4,[w1 w2]);
        data_filter1=filter(B,A,data);
        %% Define windows to compute signal and noise spectra
        % Use low f for filter =1 Hz, with 0.1, the Husid plot is not good
        
        if exist(info_peaks_t) == 0
           for i=1:length(info_peaks_t)
               if strcmp(sacfile,char(info_peaks_t{1,3}(i)))
                   tnbeg=info_peaks_t{1,4}(i);
                   tnend=info_peaks_t{1,5}(i);
                   tpbeg=info_peaks_t{1,6}(i);
                   tpend=info_peaks_t{1,7}(i);
                   tsbeg=info_peaks_t{1,8}(i);
                   tsend=info_peaks_t{1,9}(i);
                   break
               end
           end
        else
            [tnbeg,tnend,tpbeg,tpend,tsbeg,tsend,durationS2]=get_windows(time,data_filter1,S.B,S.A,S.T0,S.DELTA,S.DIST);
            fprintf(fid_duration,'%f %f %s %s\n',S.DIST,durationS2,S.KSTNM,sacfile);
        end
        
        figure(2)
        clf
        plot(time,filter(B,A,S.DATA1),'r')
        hold on
        plot([time(ceil(tpbeg/S.DELTA)) time(ceil(tpbeg/S.DELTA))],[max(filter(B,A,S.DATA1)) min(filter(B,A,S.DATA1))],'k')
        plot([time(ceil(tpend/S.DELTA)) time(ceil(tpend/S.DELTA))],[min(filter(B,A,S.DATA1)) max(filter(B,A,S.DATA1))],'k')
        plot([time(ceil(tsbeg/S.DELTA)) time(ceil(tsbeg/S.DELTA))],[max(filter(B,A,S.DATA1)) min(filter(B,A,S.DATA1))],'b')
        plot([time(ceil(tsend/S.DELTA)) time(ceil(tsend/S.DELTA))],[min(filter(B,A,S.DATA1)) max(filter(B,A,S.DATA1))],'b')
        plot([time(ceil(tnbeg/S.DELTA)) time(ceil(tnbeg/S.DELTA))],[max(filter(B,A,S.DATA1)) min(filter(B,A,S.DATA1))],'g')
        plot([time(floor(tnend/S.DELTA)) time(floor(tnend/S.DELTA))],[min(filter(B,A,S.DATA1)) max(filter(B,A,S.DATA1))],'g')
        
        pause(1)
        ButtonName = questdlg('Do you want to redo the picking?','Yes','No');
        if strcmp(ButtonName,'Yes') == 1
            [x,y,button]=ginput(6);
            fprintf(1,'%s %s %s %f %f %f %f %f %f\n',S.KSTNM,comp,sacfile,x(1),x(2),x(3),x(4),x(5),x(6));
            fprintf(fid_picks_t_new,'%s %s %s %f %f %f %f %f %f\n',S.KSTNM,comp,sacfile,x(1),x(2),x(3),x(4),x(5),x(6));
        else
            fprintf(1,'%s %s %s %f %f %f %f %f %f\n',S.KSTNM,comp,sacfile,tnbeg,tnend,tpbeg,tpend,tsbeg,tsend);
            fprintf(fid_picks_t_new,'%s %s %s %f %f %f %f %f %f\n',S.KSTNM,comp,sacfile,tnbeg,tnend,tpbeg,tpend,tsbeg,tsend);
        end
    end
end
fclose(fid_picks_t_new);
system(['mv ' fich_t_picks_new ' ' fich_t_picks]);
endfunction

function [tnbeg,tnend,tpbeg,tpend,tsbeg,tsend,durationS2]=get_windows(time,data_filter1,S_B,S_A,S_T0,S_DELTA,S_DIST)
        tsbeg=(S_T0-S_B)-5;
        durationS2=compute_5_95_nrj(time,data_filter1,(S_T0-S_B),S_DELTA);
        % Model de duree fait rapidement a la main
        % voir distance_duration.txt
        if S_DIST <= 50
            durmax=20;
            durmin=5;
        elseif S_DIST <= 500
            durmax=150;
        else
            durmax=150+0.3*(S_DIST-500);
        end
        if durationS2 <= durmax
            tsend=(S_T0-S_B)+durationS2;
        else
            tsend=(S_T0-S_B)+durmax;
        end
        if (S_A-S_B)-25 > 0
            tnbeg=(S_A-S_B)-25;
            tnend=(S_A-S_B)-5;
        else
            if (S_A-S_B)-5 > 10
                tnbeg=S_DELTA; % if 0.0 is used, the function ceil(tsbeg/S_DELTA) will return 0 which is not acceptable as matrix indice
                tnend=(S_A-S_B)-5;
            else
                if (S_DELTA*length(data_filter1)-10) > tsend
                   tnbeg=S_DELTA*length(data_filter1)-10;
                else
                   tnbeg=tsend+S_DELTA
                end
                tnend=S_DELTA*length(data_filter1);
            end
        end
        tpbeg=-99.0;
        if (S_A-S_B) > 5
            tpbeg=(S_A-S_B)-5;
        elseif (S_A-S_B) > 2.5
            tpbeg=(S_A-S_B)-2.5;
        else
            tpbeg=(S_A-S_B);
        end
        if tpbeg == -99.0
            fprintf(1,'%s\n','the beginning of the P window is not defined, paused')
            pause
        end
        tpend=tsbeg-5;
        if tpend < tpbeg
            fprintf(1,'%s\n','the end time of P window is smaller than begin time, paused')
            pause
            tpend=tsbeg-S_DELTA;
        end
        if (tsbeg < 0 || tsend < 0 || tnbeg < 0 || tnend < 0 || tsbeg > ...
                S_DELTA*length(data_filter1) || tsend > S_DELTA*length(data_filter1) ...
                || tnbeg > S_DELTA*length(data_filter1) || tnend > S_DELTA*length(data_filter1))
            fprintf(1,'%s\n','the S wave window does not fit record time')
        end
endfunction