function compute_kappa
format long
clear
addpath(genpath('~/octave'));

pname=uigetdir('/home/stephane/DATA/ON/Earthquakes/','Select an event directory (which includes spectra/ sub-directory):');
pname=strrep(strcat(pname,'/'),'//','/');

if exist(strcat(pname,'spectra'),'dir') == 0
    disp('Selected directory does not include spectra/ sub-directory. Stopping.')
    return
end

%------------------- Plot time series
fid_picks_t=fopen([pname 'spectra/list_picks.txt'],'r');
info_peaks_t=textscan(fid_picks_t,'%s %s %s %f %f %f %f %f %f');
fclose(fid_picks_t);
figure(1)
clf
plot_all_time_series(pname,info_peaks_t)
%-------------------

fid_picks_f=fopen([pname 'spectra/list_picks_f.txt'],'r');
info_peaks_f=textscan(fid_picks_f,'%s %f %f %f %s');
fclose(fid_picks_f);

fich_kappa=[pname 'spectra/kappa.out'];
if exist(fich_kappa,'file') == 2
    backup_file_pick(fich_kappa);
end
fid_kappa=fopen(fich_kappa,'w');

for i=1:length(info_peaks_f{1,1})
    spcfile=deblank(char(info_peaks_f{1,1}(i)));
    sacfile=strrep(spcfile,'.spc','');
    fmin=info_peaks_f{1,2}(i);
    fmax=info_peaks_f{1,3}(i);
    fe=info_peaks_f{1,4}(i);
    % Need the distance, read sac file and sampling rate
    S=readsac([deblank(pname) sacfile]);
    %% Read the spectra - Convert velocity to acceleration
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
    acc=A(:,2).*(2*pi*freq);
    noise=A(:,3).*(2*pi*freq);

    figure(2)
    clf
    %% plot raw signals
    % log-log plot
    subplot(1,2,1)
    loglog(freq,noise,'Color',[211/255 211/255 211/255])
    hold on
    loglog(freq,acc,'Color',[0.68 0.85 0.90])
%      loglog(freq2,noise2,'k')
    loglog(freq2,acc2,'b')
    xlim([min(freq) max(freq)])
    title([sacfile ' - ' num2str(S.DIST) ' km'],'Interpreter', 'none')
    xlabel('frequency (Hz)')
    ylabel('FAS')
    acc_smooth=smooth(acc,20);
    
    % lin-log plot
    subplot(1,2,2)
    semilogy(freq,noise,'Color',[211/255 211/255 211/255])
    hold on
    semilogy(freq,acc,'Color',[0.68 0.85 0.90])
    semilogy(freq2,noise2,'k')
    semilogy(freq2,acc2,'b')
    xlim([min(freq) max(freq)])
    xlabel('frequency (Hz)')
    ylabel('FAS')

    % Data between fe and fmax
    index_kappa=find(freq >= fe & freq < fmax);
    if isempty(index_kappa) == 0
        subplot(1,2,1)
        loglog(freq(index_kappa),acc(index_kappa),'r')
        subplot(1,2,2)
        semilogy(freq(index_kappa),acc(index_kappa),'r')
    end
    pause(1)

    clear ButtonName
    if isempty(index_kappa) == 0
%          figure(2)
%          subplot(1,2,2)
%          hold on
%          semilogy(freq(index_kappa),acc_smooth(index_kappa),'g')
        
        ButtonName = questdlg('Use this data to compute kappa?','Yes','No');
        if strcmp(ButtonName,'Yes') == 1
            % Kappa computed from smoothed data
            [kappa_ls, int_ls, kappa_rob, int_rob] = kappa1file([freq(index_kappa) acc_smooth(index_kappa)]);
            % Kappa computed from raw data
%              [kappa_ls, int_ls, kappa_rob, int_rob] = kappa1file([freq(index_kappa) acc(index_kappa)])
            fprintf(fid_kappa,'%f %f %f %f %f %s\n',S.DIST,kappa_ls,int_ls,kappa_rob,int_rob,[sacfile '.spc']);
%              fprintf(1,'%f %f %f %f %f %s\n',S.DIST,kappa_ls,int_ls,kappa_rob,int_rob,[sacfile '.spc']);
%              pause(1)
            
            subplot(1,2,1)
            loglog([min(freq(index_kappa)) max(freq(index_kappa))],[exp(int_rob-kappa_rob*pi*min(freq(index_kappa))) exp(int_rob-kappa_rob*pi*max(freq(index_kappa)))],'r')
            subplot(1,2,2)
            semilogy([min(freq(index_kappa)) max(freq(index_kappa))],[exp(int_rob-kappa_rob*pi*min(freq(index_kappa))) exp(int_rob-kappa_rob*pi*max(freq(index_kappa)))],'r')

            figure(3)
            plot(S.DIST,kappa_ls,'r+')
            plot(S.DIST,kappa_rob,'bo')
            xlim([0 1500])
            ylim([-0.05 0.3])
            hold on
        end
        figure(2)
        subplot(1,2,1)
        hold off
        subplot(1,2,2)
        hold off
    else
        % If no data between fe and fmax, no kappa computed
        ButtonName='No';
    end
end
fclose(fid_kappa);
