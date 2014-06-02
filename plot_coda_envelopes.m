function plot_coda_envelopes(S,S2,f1,f2,tpeak_env_coda,tmax,coda_env,coda_env_smoothed,noise_level,norder)

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

figure(1)
subplot(length(f1)+1,2,1)
plot(time,data,"-;comp1;")
subplot(length(f1)+1,2,2)
plot(time,data2,"-;comp2;")

Fnyquist=1/(2*S.DELTA);
time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';

for i=1:length(f1)
   w1=f1(i)/Fnyquist;
   w2=f2(i)/Fnyquist;
   [B,A] = butter(norder(i),[w1 w2]);
   
   figure(1)
   legend=[num2str(f1(i)) "-" num2str(f2(i)) " Hz"];
   subplot(length(f1)+1,2,1+2*i)
   color_tmp=[rand rand rand];
   plot(time,filter(B,A,data),'Color',color_tmp,["-;" legend ";"])
   subplot(length(f1)+1,2,2+2*i)
   plot(time,filter(B,A,data2),'Color',color_tmp)

   if noise_level(i) < 0, ymin=1.1*noise_level(i);, else, ymin=0.9*noise_level(i);, end
   if max(coda_env(:,i)) < 0, ymax=0.9*max(coda_env(:,i));, else, ymax=1.1*max(coda_env(:,i));, end
   figure(2)
   subplot(length(f1),1,i)
   plot(time,coda_env(:,i),'Color',[0.86 0.86 0.86],"-;signal;",time,coda_env_smoothed(:,i),'color',color_tmp,"-;smoothed data;",[1.11*S.O 0.99*0.95*max(time)],[noise_level(i) noise_level(i)],'color','k','LineWidth',1,"-;average pre-P noise;",[tmax(i) tmax(i)],[ymin ymax],"-;signal>3*noise;",[tpeak_env_coda(i) tpeak_env_coda(i)],[ymin ymax],"-;peak coda env.;")
   ylim([ymin ymax])

end

figure(1)
%  saveas(gcf,['./filtered_time_series.jpg'],'jpeg');
print -djpeg filtered_time_series.jpg

figure(2)
%  saveas(gcf,['./coda_envelopes.jpg'],'jpeg');
print -djpeg coda_envelopes.jpg

endfunction