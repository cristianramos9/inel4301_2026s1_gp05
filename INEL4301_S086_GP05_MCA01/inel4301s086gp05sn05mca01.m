%inel4301s000gp00sn00mca01
%Prof. Domingo Antonio Rodriguez
%Name of Student, sn01, DTMF_WORD, DITS, %XX of Contribution
%Name of Student, sn02, DTMF_WORD, DITS, %XX of Contribution
%Name of Student, sn03, DTMF_WORD, DITS, %XX of Contribution
%Name of Student, sn04, DTMF_WORD, DITS, %XX of Contribution
%Cristian Ramos Ramos, sn05, dtmfau, DITS: 1 and 5, %20 of Contribution
%**************************Prof. D. Rodriguez***********************
%*********DSB-SC COMMUNICATIONS SYSTEMS AND GAUSSIAN NOISE**********
%*********ANALOG AND DIGITAL COMMUNICATIONS SIMULATIONS*************
%*******************************************************************
clear all                            %Clear working space
close all
%*******************************************************************
disp('******************************')
disp('Begin Simulation for Underwater Communications.')
disp('******************************')
%*******************************************************************
%********************Wanted Signal Input/Output*********************
%*******************************************************************
[sig,Fs]=audioread('star_story.wav');  %Get overall sampling frequency
save dsig.txt sig -ascii               %Save signal as dsig.txt file
load -ascii dsig.txt;                %Load ascii file dsig.txt
wsiz=length(dsig);                   %Get the length of the "dsig" signal
%*******************************************************************
%**************************Parameter Setting************************
%*******************************************************************
Ts=1/Fs;                          %Sampling time or temporal resolution
Nq=wsiz;                          %Number of samples of wanted signal
Vq=Nq*Ts;                         %Observation window or signal duration
fincq=1/Vq;                       %Frequency resolution of wanted signal
fkq=-(Fs/2):fincq:+(Fs/2)-fincq;  %Frequency axis of wanted signal
tq=0:Ts:Vq-Ts;                    %Time axis of wanted signal
%
fV=4000;                          %Cut-off frequency of wanted signal
Wn=(2*fV/Fs);                     %Normalized cut-off frequency
Lp=(1/567);                       %Proportionality factor
Nh=2*floor(Lp*Nq)+1;              %Low-pass filter length
hL=fir1(Nh,Wn);                   %Impulse response of low-pass filter
dsigL=conv(dsig,hL);              %Filtered wanted signal up to fV=4000 Hz
%
f1=4400;                             %Input interference frequency
f2=4900;                             %Input interference frequency
f3=5400;                             %Input interference frequency
g1=cos(2*pi*f1*tq);                  %Unwanted input signal 
g2=cos(2*pi*f2*tq);                  %Unwanted input signal
g3=cos(2*pi*f3*tq);                  %Unwanted input signal
g=(1/80)*g1+(1/160)*g2+(1/240)*g3;   %Sum of unwanted signals
gmax=max(g);                         %Maximum value of sum
g=transpose((1/gmax)*g);             %Normalized interference signal
%********************************************************************
%*********************SIGNALS AND SYSTEMS MODELING*******************
fc=8640;                        %Carrier frequency
c=transpose((cos(2*pi*fc*tq))); %Carrier signal
wsig=dsigL(1:Nq,1);             %Wanted input signal
xm=wsig+g;                      %Sum of wanted and unwanted signals
xc=xm.*c;                       %Modulated wanted signal plus interference
%*********************************************************************
%*************************BAND-PASS FILTER****************************
Nb=1001;                           %Order of band-pass filter
Wb=[(2*4640)/Fs  (2*12640)/Fs];    %Vector of normilized cut-off freqs.
hB=fir1(Nb-1,Wb,'bandpass');       %Impulse response of band-pass filter
%*********************************************************************
%*************************CHANNEL INPUT SIGNAL************************
yci=conv(xc,hB);                    %Removal of unwanted interferences 
Nyci=length(yci);                   %Order of channel input signal
tyci=0:Ts:(Nyci-1)*Ts;              %Time axis for yci(t)
fryci=1/(Nyci*Ts);                  %Frequency resolution for Yci(f)
fayci=-(Fs/2):fryci:+(Fs/2)-fryci;  %Frequency axis of wanted signal
%*********************************************************************
%*************************CHANNEL NOISE*******************************
nsig=randn(size(yci));          %Generation of AWGN signal: noise signal
nnsig=(1/max(abs(nsig)))*nsig;  %Normalized noise signal
nATT=sqrt(1/24000);                %Noise attenuation about -44 dB SNR
nsig=nATT*nnsig;                %Attenuated noise present in the channel
%*********************************************************************
%*************************CHANNEL OUTPUT******************************
yco=yci+nsig;           %Channel output signal with Gaussian noise
tqn=0:Ts:(length(yco)-1)*Ts;        %New time axis for output signal
cn=transpose((cos(2*pi*fc*tqn)));   %New carrier signal
yd=yco.*(2*cn);                     %Demodulation stage
%*********************************************************************
%*************************LOW-PASS FILTER*****************************
Nlow=663;              %Designed low-pass filter order
Wlow=(2*fV)/Fs;         %Normalized cut-off frequency
hlow=fir1(Nlow-1,Wlow); %FIR low-pass impulse response function
taxis_hlow=0:Ts:(length(hlow)-1)*Ts;
fres_hlow=1/((length(hlow))*Ts);
faxis_hlow=(-Fs/2):fres_hlow:(+Fs/2)-fres_hlow;
%*********************************************************************
%*************************RECEIVER OUTPUT*****************************
xr=conv(yd,hlow);       %Low pass filtering at the receiver
taxis_xr=0:Ts:(length(xr)-1)*Ts;
fres_xr=1/((length(xr))*Ts);
faxis_xr=(-Fs/2):fres_xr:(+Fs/2)-fres_xr;
%********************************SPECTRA*****************************
%********************************************************************
fwsig=fft(wsig);                %Spectrum of wanted signal (Fund. Reg.)
sfwsig=fftshift(fwsig);         %Spectrum of wanted signal (Princ. Reg.)
asfwsig=abs(sfwsig);            %Magnitude of spectrum of wanted signal
%
fxc=fft(xc);                  %Spectrum of modd. signal (Fund. Reg.)
sfxc=fftshift(fxc);           %Spectrum of modd. signal (Princ. Reg.)
asfxc=abs(sfxc);              %Magnitude of spectrum of modd. signal
%
fyci=fft(yci);                  %Spectrum of input signal (Fund. Reg.)
sfyci=fftshift(fyci);           %Spectrum of input signal (Princ. Reg.)
asfyci=abs(sfyci);              %Magnitude of spectrum of input signal
%
fyco=fft(yco);                  %Spectrum of output signal (Fund. Reg.)
sfyco=fftshift(fyco);           %Spectrum of output signal (Princ. Reg.)
asfyco=abs(sfyco);              %Magnitude of spectrum of output signal
%
fyd=fft(yd);                    %Spectrum of output signal (Fund. Reg.)
sfyd=fftshift(fyd);             %Spectrum of output signal (Princ. Reg.)
asfyd=abs(sfyd);                %Magnitude of spectrum of output signal
%
fnsig=fft(nsig);                %Spectrum of noise signal (Fund. Reg.)
sfnsig=fftshift(fnsig);         %Spectrum of noise signal (Princ. Reg.)
asfnsig=abs(sfnsig);            %Magnitude of spectrum of noise signal
%
fhlow=fft(hlow);                %Spectrum of IRF (Fund. Region)
sfhlow=fftshift(fhlow);         %Spectrum of IRF (Princ. Region)
asfhlow=abs(sfhlow);            %Magnitude of spectrum IRF
%
fxr=fft(xr);                    %Spectrum of xr(t) (Fund. Region)
sfxr=fftshift(fxr);             %Spectrum of xr(t) (Princ. Region)
asfxr=abs(sfxr);                %Magnitude of spectrum IRF
%******************************PLOTS*********************************
plot(tq,wsig);
axis([0 (length(wsig)*Ts) min(wsig)  max(wsig)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Wanted signal s(t) in the Time Domain')
grid
%
figure
plot(tq,g);
axis([0 (length(g)*Ts) min(g)  max(g)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Interference signal g(t) in the Time Domain')
grid
%
figure
plot(tq,xm);
axis([0 (length(xm)*Ts) min(xm)  max(xm)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Modulating signal xm(t) in the Time Domain')
grid
%
figure
plot(tq,xc)
axis([0 (length(xc)*Ts) min(xc)  max(xc)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Modulated Signal xc(t)')
grid
%
figure
plot(tyci,nnsig)
axis([0 (length(nnsig)*Ts) min(nnsig)  max(nnsig)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Normalized Noise Signal nn(t)')
grid
%
figure
plot(tyci,yci)
axis([0 (length(yci)*Ts) min(yci) max(yci)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Input Signal yci(t)')
grid
%
figure
plot(tyci,yco)
axis([0 (length(yco)*Ts) min(yco)  max(yco)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Channel Output Signal yco(t)')
grid
%
figure
plot(fkq,asfwsig)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum S(f) of Wanted Signal s(t)')
grid
%
figure
plot(fkq,asfxc)
xlabel('Freq in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum Xc(f)of Modulated Signal xc(t)')
grid
%
figure
plot(fayci,asfyci)
xlabel('Freq in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum Yci(f) of Channel Input Signal yci(t)')
grid
%
figure
plot(fayci,asfyco)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum Yco(f) of Channel Output Signal yco(t)')
grid
%
figure
plot(fayci,asfyd)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum Yd(f) of Demodulator Output yd(t)')
grid
%
figure
plot(taxis_hlow,hlow);
axis([0 (length(hlow)*Ts) min(hlow)  max(hlow)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Impulse Response Function of Receiver''s Low-Pass Filter')
grid
%
figure
plot(faxis_hlow,asfhlow)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Impulse Response Function hlow(t)')
grid
%
figure
plot(taxis_xr,xr);
axis([0 (length(xr)*Ts) min(xr)  max(xr)]);
xlabel('Time in Seconds')
ylabel('Amplitude')
title('Receiver''s Low-Pass Filter Output Signal xr(t)')
grid
%
figure
plot(faxis_xr,asfxr)
xlabel('Frequency in Hertzs')
ylabel('Magnitude')
title('Magnitude of Spectrum of Receiver''s Output xr(t)')
grid
%
%*******************************OUTPUTS******************************
sound(wsig,Fs)  %Sound of Wanted Signal
disp('Sound of Wanted Signal')
pause(20)
sound(g,Fs)     %Sound of Interference Signal
disp('Sound of Interference Signal')
pause(20)
sound(xm,Fs)    %Sound of Modulating Signal
disp('Sound of Modulating Signal')
pause(20)
sound(xc,Fs)    %Sound of Modulated Signal
disp('Sound of Modulated Signal')
pause(20)
sound(yci,Fs)   %Sound of Channel's Input Signal
disp('Sound of Channel''s Input Signal')
pause(20)
sound(yco,Fs)   %Sound of Channel's Output Signal
disp('Sound of Channel''s Output Signal')
pause(20)
sound(xr,Fs)   %Sound of Channel's Output Signal
disp('Sound of Receiver''s Output Signal')
%********************************************************************
disp('******************************')
disp('Simulation ended succesfully.')
disp('******************************')
