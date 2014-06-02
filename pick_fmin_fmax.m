function pick_fmin_fmax
format long
clear
%  close all
addpath(genpath('~/octave'));

pname=uigetdir('/home/stephane/DATA/ON/Earthquakes/','Select an event directory (which includes spectra/ sub-directory):');
pname=strrep(strcat(pname,'/'),'//','/');

if exist(strcat(pname,'spectra'),'dir') == 0
    disp('Selected directory does not include spectra/ sub-directory. Stopping.')
    return
end

%% Pick fmin in log-log, and fe and fmax in lin-log space
%% Load the picks or open a new file
fich_picks_f=[pname 'spectra/list_picks_f.txt'];
if exist(fich_picks_f,'file') == 2
    backup_file_pick(fich_picks_f);
    fid_picks_f=fopen(fich_picks_f,'r');
    info_peaks_f=textscan(fid_picks_f,'%s %f %f %f %s');
    fclose(fid_picks_f);
else
    info_peaks_f=[];
end
fich_picks_f_new=[pname 'spectra/list_picks_f.tmp'];
if exist(fich_picks_f_new,'file') == 2
    system(['rm -f ' fich_picks_f_new]);
end

fid_picks_t=fopen([pname 'spectra/list_picks.txt'],'r');
info_peaks_t=textscan(fid_picks_t,'%s %s %s %f %f %f %f %f %f');
fclose(fid_picks_t);

figure(1)
clf
plot_all_time_series(pname,info_peaks_t)

figure(2)
clf
figure(3)
clf

for i=1:length(info_peaks_t{1,1})
    %% Read the time series
    sacfile=deblank(char(info_peaks_t{1,3}(i)))
    % Need the distance, read sac file and sampling rate
    S=readsac([deblank(pname) sacfile]);
    Fnyquist=1/(2*S.DELTA);
    %% Read the spectra
    %-------------------------------------------------------
    % multi-taper spectra
    C=load([pname 'spectra/' sacfile '.spc_n']);
    D=load([pname 'spectra/' sacfile '.spc_s']);
    index=find(D(:,1)>0);
    freq2=D(index,1);
    acc2=D(index,2).*(2*pi*freq2);
    noise2=spline(C(:,1),C(:,2),freq2).*(2*pi*freq2);
    %-------------------------------------------------------
    A=load([pname 'spectra/' sacfile '.spc']);
    freq=A(:,1);
    % Convert velocity to acceleration
    acc=A(:,2).*(2*pi*freq);
    noise=A(:,3).*(2*pi*freq);
    % S/N >= 5
    index=find(smooth(acc,20) >= 5*smooth(noise,20));

    %% plot raw signals
    % log-log plot
    figure(2)
    loglog(freq,noise,'Color',[211/255 211/255 211/255])
    hold on
    loglog(freq,acc,'Color',[0.68 0.85 0.90])
    loglog([Fnyquist Fnyquist],[max(acc) min(noise)],'b')
    loglog([Fnyquist/2 Fnyquist/2],[max(acc) min(noise)],'b')
    loglog(freq2,noise2,'k')
    loglog(freq2,acc2,'b')
    %% Plot smoothed signals
    %loglog(freq,smooth(noise,20),'Color',[190/255 190/255 190/255])
    xlim([min(freq) max(freq)])
    title([sacfile ' - ' num2str(S.DIST) ' km'],'Interpreter', 'none')
    %loglog(freq,smooth(acc,20),'b')
    acc_smooth=smooth(acc,20);
    %loglog(freq(index),acc_smooth(index),'k.')

    % lin-log plot
    figure(3)
    semilogy(freq,noise,'Color',[211/255 211/255 211/255])
    hold on
    semilogy(freq,acc,'Color',[0.68 0.85 0.90])
    semilogy([Fnyquist Fnyquist],[max(acc) min(noise)],'b')
    semilogy([Fnyquist/2 Fnyquist/2],[max(acc) min(noise)],'b')
    semilogy(freq2,noise2,'k')
    semilogy(freq2,acc2,'b')
    %% Plot smoothed signals
    %semilogy(freq,smooth(noise,20),'Color',[190/255 190/255 190/255])
    xlim([min(freq) max(freq)])
    %semilogy(freq,smooth(acc,20),'b')
    %semilogy(freq(index),acc_smooth(index),'k.')
    
    %% Load the picks and test if current file already has been picked
    [fmin,fmax,fe]=set_f_picks(info_peaks_f,sacfile,index,freq);
    
    figure(2)
    loglog([fmin fmin],[max(acc) min(noise)],'c')
    loglog([fe fe],[max(acc) min(noise)],'r')
    if fmax > 0
        loglog([fmax fmax],[max(acc) min(noise)],'g')
    end
    figure(3)
    semilogy([fmin fmin],[max(acc) min(noise)],'c')
    semilogy([fe fe],[max(acc) min(noise)],'r')
    semilogy([fmax fmax],[max(acc) min(noise)],'g')
    if fmax < 0 || fmax < fmin
        xlim([fmax  max(freq)])
    end
    xlim([0 50])
    pause(1)
    
    clear ButtonName
    ButtonName = questdlg('Redo the picking?','Yes','No');
    if strcmp(ButtonName,'Yes') == 1
        clear x y button
        figure(2)
        [x,y,button]=ginput(1);
        fmin=x(1);
        clear x y button
        figure(3)
        [x,y,button]=ginput(2);
        fe=x(1);
        fmax=x(2);
    end
    
    fid_picks_f_new=fopen(fich_picks_f_new,'a');
    fprintf(fid_picks_f_new,'%s %f %f %f %s\n',[sacfile '.spc'],fmin,fmax,fe,'No');
    fclose(fid_picks_f_new);
    
    figure(2)
    hold off    
    figure(3)
    hold off
    
end
system(['mv ' fich_picks_f_new ' ' fich_picks_f]);
endfunction

function [fmin,fmax,fe]=set_f_picks(info_peaks_f,sacfile,index,freq)
    clear fmin fmax fe
    if ~isempty(info_peaks_f)
%          sacfile
%          info_peaks_f{1,1}
        for k=1:length(info_peaks_f{1,1})
            if strcmp(strcat(sacfile,'.spc'),deblank(char(info_peaks_f{1,1}(k)))) == 1
                strcat(sacfile,'.spc')
                fmin=info_peaks_f{1,2}(k);
                fmax=info_peaks_f{1,3}(k);
                fe=info_peaks_f{1,4}(k);
                break
            end
        end
        if (isempty(fmin) || isempty(fmax) || isempty(fe))
            disp([sacfile ' not found in list_picks_f.txt'])
            break
        end
    else
        disp('case no picks f')
        if isempty(index)
            % Case no data with S/N >= 5, set dummy picks
            % Need fmax < fmin in order to exclude the data in the next steps (see analyse_spectra)
            fmax=min(freq)
            fmin=fmax*1.1;
            fe=fmax*1.2;
        else
            % New data set fmin, fe, and fmax
            % Use min frequency with S/N >= 5
            fmin=min(freq(index));
            % Use max frequency with S/N >= 5
            fmax=max(freq(index));
            % Use fmin+20%
            fe=fmin*1.20;
        end
    end
endfunction
