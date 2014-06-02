format long
clear all
close all

fid_duration=fopen('distance_duration.txt','a');

[fname,pname]=uigetfile('list_*.txt','Select a list of filenames for one event','/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3');

answer=inputdlg('Conversion factor between data unit and m/s','Data Unit',1,{'1e-9'});
conv=str2double(answer);

mkdir(pname,'spectra')
fich_picks=[pname '/spectra/list_picks.txt'];
if exist(fich_picks,'file') == 2
    fid_picks=fopen(fich_picks,'r');
    write_picks='n';
else
    fid_picks=fopen(fich_picks,'w');
    write_picks='y';
end

list_sacfiles=read_list([pname fname]);

figure(1)
clf
iHHZcomp=0;
for i=1:size(list_sacfiles,1)
    fprintf(1,'%s\n',[pname list_sacfiles(i,:)]);
    %% Read SAC data
    S=readsac([pname list_sacfiles(i,:)]);
    %% Bidouille pour les noms de fichier spectre, dans le cas des donnees de Marcelo, les noms de fichiers
    %% et les noms de composantes dans les headers sac ne sont pas homogenes.
    %% Hypothese fichiers .2.sac.vel et .5.sac.vel sont N, et fichiers .3.sac.vel et .6.sac.vel sont E
    %% De plus si comp=radial je modeifie pour nord et transversal pour est
    if strcmp(S.KCMPNM(1:2),'HH') ~= 1
        if strcmp(S.KCMPNM,'vertical') == 1
            comp='z';
        elseif strcmp(S.KCMPNM,'north') == 1
            comp='n';
        elseif strcmp(S.KCMPNM,'east') == 1
            comp='e';
        else
            % If component name not recognised, look at filename
            pos=strfind(list_sacfiles(i,:),'.');
            test=list_sacfiles(i,pos(length(pos)-1)-1:pos(length(pos)-1)-1);
            if (strcmp(test,'1') == 1 || strcmp(test,'4') == 1 || strcmp(test,'z') == 1 || strcmp(test,'Z') == 1)
                comp='z';
            elseif (strcmp(test,'2') == 1 || strcmp(test,'5') == 1 || strcmp(test,'n') == 1 || strcmp(test,'N') == 1 || strcmp(test,'r') == 1)
                comp='n';
            elseif (strcmp(test,'3') == 1 || strcmp(test,'6') == 1 || strcmp(test,'e') == 1 || strcmp(test,'E') == 1 || strcmp(test,'t') == 1)
                comp='e';
            else
                uiwait(msgbox('Component name in file name not recognized','Warning','warn','modal'))
            end
        end
    else
        if strcmp(S.KCMPNM,'HHZ') == 1
            comp='z';
        elseif strcmp(S.KCMPNM,'HHN') == 1
            comp='n';
        elseif strcmp(S.KCMPNM,'HHE') == 1
            comp='e';
        end
    end
    if (~isempty(S) && isnan(S.A) ~=1 && isnan(S.T0) ~= 1)
        %fprintf(1,'%s\t%f\t%f\t%f\t%f\t%f\n',list_sacfiles(i,:),S.DELTA,S.B,S.T0,S.E,S.DELTA*length(S.DATA1));
        time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
        %% demean, detrend
        data_tmp=detrend(S.DATA1,'constant');
        data=detrend(data_tmp);
        clear data_tmp
        Fnyquist=1/(2*S.DELTA);
        w1=1.0/Fnyquist;
        %w2=min(0.99999999999,10.0/Fnyquist);
        w2=0.99999999999;
        [B,A] = butter(4,[w1 w2]);
        data_filter=filter(B,A,data);
        %% Define windows to compute signal and noise spectra
        % Use low f for filter =1 Hz, with 0.1, the Husid plot is not good
        tsbeg=(S.T0-S.B);
        %         tsend=(S.T0-S.B)+20;
        durationS2=compute_5_95_nrj(time,data_filter,(S.T0-S.B),S.DELTA);
        fprintf(fid_duration,'%f %f\n',S.DIST,durationS2);
        %% Change the low f of filter for FFT
        w1=0.1/Fnyquist;
        clear B A data_filter
        [B,A] = butter(4,[w1 w2]);
        data_filter=filter(B,A,data);
        % Model de duree fait rapidement a la main
        % voir distance_duration.txt
        if S.DIST <= 50
            durmax=20;
        elseif S.DIST <= 500
            durmax=150;
        else
            durmax=150+0.3*(S.DIST-500);
        end
        if durationS2 <= durmax
            tsend=(S.T0-S.B)+durationS2;
        else
            tsend=(S.T0-S.B)+durmax;
        end
        if (S.A-S.B)-25 > 0
            tnbeg=(S.A-S.B)-25;
            tnend=(S.A-S.B)-5;
        else
            if (S.A-S.B)-5 > 10
                tnbeg=S.DELTA; % if 0.0 is used, the function ceil(tsbeg/S.DELTA) will return 0 which is not acceptable as matrix indice
                tnend=(S.A-S.B)-5;
            else
                tnbeg=S.DELTA*length(data)-10;
                tnend=S.DELTA*length(data);
            end
        end
        tpbeg=-99.0;
        if (S.A-S.B) > 0
            tpbeg=(S.A-S.B);
        end
        if tpbeg == -99.0
            fprintf(1,'%s\n','the beginning of the P window is not defined, paused')
            pause
        end
        tpend=tsbeg-5;
        if tpend < tpbeg
            fprintf(1,'%s\n','the end time of P window is smaller than begin time, paused')
            pause
            tpend=tsbeg-S.DELTA;
        end
        if (tsbeg < 0 || tsend < 0 || tnbeg < 0 || tnend < 0 || tsbeg > ...
                S.DELTA*length(data) || tsend > S.DELTA*length(data) ...
                || tnbeg > S.DELTA*length(data) || tnend > S.DELTA*length(data))
            fprintf(1,'%s\n','the S wave window does not fit record time')
        end
        if strcmp(write_picks,'y') == 1
            fprintf(fid_picks,'%s %s %s %f %f %f %f %f %f\n',S.KSTNM,comp,list_sacfiles(i,:),tnbeg,tnend,tpbeg,tpend,tsbeg,tsend);
        end
        %% Plot time series for vertical components of all records and filtered time series between 1 and 10 Hz

        figure(1)
        if (strcmp(comp,'z') == 1)
            iHHZcomp=iHHZcomp+1;
            %subplot(ceil(size(list_sacfiles,1)/9),3,iHHZcomp)
            figure(iHHZcomp)
            plot(time,data,'b')
            hold on
            plot(time,data_filter,'r')
            plot([time(ceil(tpbeg/S.DELTA)) time(ceil(tpbeg/S.DELTA))],[max(data) min(data)],'k')
            plot([time(ceil(tpend/S.DELTA)) time(ceil(tpend/S.DELTA))],[min(data) max(data)],'k')
            plot([time(ceil(tsbeg/S.DELTA)) time(ceil(tsbeg/S.DELTA))],[max(data) min(data)],'k:')
            plot([time(ceil(tsend/S.DELTA)) time(ceil(tsend/S.DELTA))],[min(data) max(data)],'k:')
            plot([time(ceil(tnbeg/S.DELTA)) time(ceil(tnbeg/S.DELTA))],[max(data) min(data)],'k--')
            plot([time(ceil(tnend/S.DELTA)) time(ceil(tnend/S.DELTA))],[min(data) max(data)],'k--')
            title([S.KSTNM ' - ' S.KCMPNM ' - ' num2str(S.DIST) ' km'],'Interpreter', 'none')
            xlabel('time (s)')
            ylabel('velocity (m/s)')
            hold off
        end
        %% Test calcul mr
        %% Compute mr
%         figure(3)
%         clf
        % !!! Warning conversion nm/s -> micron/s
        % Two ways to compute the mR (Assumpcao 1983), using displacement
        % or using velocity (Marcelo RSB meeting march 2013)
        % amplitude displacement
%         plot(time,cumtrapz(time,data_filter,'b')
%         pause
%         tmp=ginput(2);
%         per=abs(tmp(3)-tmp(1))*S.DELTA;
%         amp=abs(tmp(2));
%         mr(i)=log10(1e-3*amp/per)+2.3*log10(S.DIST)-1.48
        % Amplitude velocity
%         plot(time,data_filter,'b')
%         pause
%         tmp=ginput(2);
%         amp=abs(tmp(2)-tmp(4))/2;
%         mr(i)=log10(1e-3*amp)+2.3*log10(S.DIST)-2.28;
        %% Compute Fourier spectra
        % Fourier transform of velocity
        NdatNoise=length(data(ceil(tnbeg/S.DELTA):ceil(tnend/S.DELTA)));
        NdatS=length(data(ceil(tsbeg/S.DELTA):ceil(tsend/S.DELTA)));
        if NdatS >= NdatNoise
            NFFT=2^nextpow2(NdatS);
        else
            NFFT=2^nextpow2(NdatNoise);
        end
        [f2,spec2,df,norm]=fourier_transform(data_filter(ceil(tsbeg/S.DELTA):ceil(tsend/S.DELTA)),S.DELTA,NFFT,1,length(data(ceil(tsbeg/S.DELTA):ceil(tsend/S.DELTA))));
        % Take out f=0 Hz
        f=f2(2:length(f2));
        vel_spec_S=abs(spec2(2:NFFT/2+1));
        [fn2,specn2,df,norm]=fourier_transform(data_filter(ceil(tnbeg/S.DELTA):ceil(tnend/S.DELTA)),S.DELTA,NFFT,1,length(data(ceil(tnbeg/S.DELTA):ceil(tnend/S.DELTA))));
        vel_spec_N=abs(specn2(2:NFFT/2+1));
        if strcmp(S.IDEP,'IACC') == 1
            fprintf(1,'%s\n','Data are not velocity data');
            vel_spec_S=abs(spec2(2:NFFT/2+1))./(2*pi*f)';
            vel_spec_N=abs(specn2(2:NFFT/2+1))./(2*pi*f)';
        end
        % Data are in nm/s
        
        vel_spec_S=conv*vel_spec_S;
        vel_spec_N=conv*vel_spec_N;
        fid_spectra=fopen([pname '/spectra/' deblank(list_sacfiles(i,:)) '.spc'],'w');
        fprintf(fid_spectra,'%e %e %e\n',[f' vel_spec_S vel_spec_N]');
        fclose(fid_spectra);

    end
end
fclose(fid_picks);
saveas(gcf,[pname '/spectra/plot_time_series'],'fig');
saveas(gcf,[pname '/spectra/plot_time_series'],'jpeg');
fclose(fid_duration);
figure(1)
hold off
