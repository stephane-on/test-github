function [f,spec,df,norm]=fourier_transform(amplitude,dt,NFFT,imin,imax)

B=fft(amplitude(imin:imax),NFFT);
df=1/(dt*NFFT);
f=([1:NFFT/2+1]-1)*df;
%norm=dt*sqrt(size(amplitude,2));
%norm=dt*size(amplitude,2);
norm=dt;
spec=B.*norm;
