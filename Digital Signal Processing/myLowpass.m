%%%%%%%%%%%%%%%%%%%%%%%%% P R O J E C T in Matlab
% Stathopoulos Manos
% Axelos Christos

% Ylopoihthikan ola ta meri tis ergasia ektos apo to upoerwtima sto B meros
% tis 1is askisis gia N = 512, afou to sxima den itan to epithimito

%%%%%%%%%%%%%%%%%%% M E R O S  A ' SXEDIASMOS FILTRWN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A1. FIR Filter
% - Efarmozoume ti methologia gia to FIR filtro, to opoio prokuptei 5ou
% vathmou
% - Apta dedomena vriskoume tin suxnotita apokopis (ws+wp)/2
% - Prokuptei vathmos filtrou N = 40 kai einai megaluteros aptou IIR
% afto simainei oti to filtro exei megaluteri kathisteri sto na katevei apto
% ena sto miden
clear; clear all; close all;
fsample = 16000;
As = 40;Ap = 1;%db

fs = 0.15; fp = 0.05;
ws=0.3*pi; wp=0.1*pi;

% find min filter length M from hamming's equation:
N = ceil( 4/(fs-fp) );
%discrete time fir-hamming-lowpass coefficients
filtCoe = fir1(N, (wp+ws)/2, 'low', hamming(N+1));

%% A2. IIR Filter
% - Kataskevi IIR filtrou me vasi tou digrammikou metasximatismou, gnwstis 
% methodologias
% - Prokuptei filtro vathmou N = 5 poly mikrotero apto FIR
% s=jW...s=e^jw
WP = 2*fsample*tan(wp/2);%convert to analog
WS = 2*fsample*tan(ws/2);
[N,Wn] = buttord(WP, WS, Ap, As, 's');
[z, p] = butter(N, Wn, 'low', 's');
[num, den] = bilinear(z, p, fsample);


%% A3.Impulse response
hold on; figure (1);  title('Impulse Response'); grid('ON');
xlabel('Samples'); ylabel('Amplitude'); 

[firY, firX] = impz(filtCoe);
[iirY, iirX] = impz(num, den);
stem(firX, firY, 'g', 'filled');
hold on;
stem(iirX, iirY, 'b', 'filled');
legend('FIR', 'IIR');
hold off;


%% A4.Step response
figure (2); hold on; title('Step Response'); grid('ON');
xlabel('Samples'); ylabel('Amplitude');

[firY, firX] = stepz(filtCoe);
[iirY, iirX] = stepz(num, den);
stem(firX, firY, 'g', 'filled');
hold on;
stem(iirX, iirY, 'b', 'filled');
legend('FIR', 'IIR');
hold off;


%% A5. Frequency response H(e^jw)
% - Sxima ilopoieitai se DB
% - An koitousame tin metavoli tou platous apo 1 se miden paratiroume oti
% to IIR filtro kanei tin metavasi aftin se stenotero euros suxnotitwn,
% pragma pou einai epithimito. To FIR/hamming omws exei megaluteri
% statherothta opws dixnei to sxima twn DB, afou stin xeiroteri 
% periptwsi ftanei ta -70dB
% - Oi suxnotites tou sximatos einai ston aksona tou x kai einai
% kanonikopoihmenes kata Nyquist, diladi kata w/pi

%frequency domain FIR-IIR
[hFir,wFir]=freqz(filtCoe,1,256);
[hIir,wIir]=freqz(num,den,256);

%amplitute-gain in dB
hFirDB = 20*log10(abs(hFir));
hIirDB = 20*log10(abs(hIir));

%plot
figure(3); hold on; title('FIR-IIR Frequency Response');
grid('ON'); ylabel('Gain');xlabel('Normalised frequency(f)');
plot(wFir/(pi), hFirDB, 'g');%normalized frequency(Nyquist)
hold on;
plot(wIir/(pi), hIirDB,'b');%Normalized frequency(Nyquist)
legend('FIR', 'IIR');
hold off;
%% A6. Group delay
% -d/dw{ -tan^-1(  Im(H(e^jw))/Re(H(e^jw)) }
% paratiroume oti gia fsample <1641.25 oso i normalised freq teinei sto 1
% to group delay teinei sto + apeiro enw apo tin oriaki syxnotita fsample 
% fsample=1641.5 allazei to prosimo tou group delay kai teinei sto - apeiro


[gdFir, wFir] = grpdelay(filtCoe);%default 8192 samples
[gdIir, wIir] = grpdelay(num, den);%default 8192 samples
%plot
figure(4); hold on; title('Group delay');
grid('ON');xlabel('Normalised frequency(f)');
plot(wFir/pi,gdFir, 'g');
hold on;
plot(wIir/pi, gdIir,'b');
legend('FIR', 'IIR');
hold off;

%% A7. Zeros/Poles
% - Afto afora mono to IIR, afou to FIR den exei polous kai ara
% einai efstathi

[b, a] = zplane(z, p);
hold on; title('IIR Zeros/Poles');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   M E R O S   B ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%** Ypologismos DFT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Edw upologizoume tous DFT mesw tis fft gia diaforetika N

%% B1.1 Input Signal 
%Configure signal x[n] = A1 * cos(w1*n) + A2*cos(w2*n)
A1 = 1;
A2 = 0.50;

l1 = length('stathopoulos');
l2 = length('axelos');

w1 = pi*mod( (10/7.5)*( max(l1,l2)/(l1+l2) ), 1 );
w2 = mod( w1 + pi/4, pi);
Dw = abs(w1-w2);

%% B1.2 Windowing and FFT

Nbig = 2^14;
Nbins = 16000;

%L=16
rectw = rectwin(16);
hamw = hamming(16);


%% B1.3 DFT with windows

figure(5); hold on; title('FFT of signal x[n]');
grid('ON'); xlabel('frequency(bins)'); ylabel('Magnitude(db)');

for n = 1:16  
  x16(n) = (A1*cos(w1*n) + A2*cos(w2*n) ) *rectw(n);  
end
plot( 20*log10(abs(fftshift(fft(x16, Nbins)))), 'b' ); hold on;

%L=64
rectw = rectwin(64);
for n = 1:64
  x64(n) = (A1*cos(w1*n) + A2*cos(w2*n) ) *rectw(n); 
end
plot( 20*log10( abs(fftshift( fft(x64, Nbins)) )), 'g'); hold on;
% 
% %L=512
% rectw = rectwin(512);
% for n = 1:512  
%   x512(n) = (A1*cos(w1*n) + A2*cos(w2*n) ) *rectw(n);  
% end
% plot( 20*log10( fftshift(abs(fft(x512, Nbins)) )) , 'r');
%exw problima edw

L=2^14;
rectw = rectwin(2^14);
for n = 1:2^14
  xbig(n) = (A1*cos(w1*n) + A2*cos(w2*n) ) *rectw(n);  
end
plot( 20*log10(abs(fftshift( fft(xbig, Nbins)))), 'c'); 

for n = 1:16  
  x16(n) = (A1*cos(w1*n) + A2*cos(w2*n) ) *hamw(n);  
end
plot( 20*log10(abs(fftshift( fft(x16, Nbins)))), 'y' );


legend('rect,L=16', 'rect,L=64', 'rect,L=16384', 'ham,L=16');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% MEROS B2  Voice recording, spectrogram %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fs: sampling frequency
% noverlap: number of samples each segment overlaps
% F: frequency interpolation
% den exei sumperilifthei o thorivos sto fasmatograma
%%%%% for fsample = 22KHz
% wind1 = hamming(2200); % points per segment for 22KHZ
% wind2 = hamming(220);
% noverlap = 110;

%%%%% for fsample = 16KHz %%%%%%%%%%
freqs = 0:10:5000;
fsample = 16000;
wind1 = hamming(1600); % points per segment for 16KHZ
wind2 = hamming(160); % points per segment, oso afksanetai exoume kaluteri poiotita
noverlap = 80; % overlap with 80 points per segment

%% B2.1 Record 20secs 
% -Se periptwsi pou den exoume ixografisei kati
% r = audiorecorder(fsample, 16, 1);
% recordblocking(r, 20);     % speak into microphone...
% p = play(r);   % listen to complete recording
% mySpeech = getaudiodata(r, 'double'); 
% audiowrite('voice_no_filter16KHZ.wav', mySpeech, fsample); 

% -Se periptwsi pou exoume idi ixografisei kati
[mySpeech, fsample] = audioread('voice_no_filter16KHZ.wav');
%% B2.2 Sound spectrograms, without lowpass
% - Xrisimopoiw tin sinartisi imagsc wste na metaferw to sxima se morfi 
% pdf. Alliws mporousa apla na xrisimopoihsw tin spectrogram xwris outputs
% - Gia L = 1600, paratiroume oti to dinatotero sima vrisketai sto diastima
% suxnotitwn [70, 180]. 
% - Sima iparxei kai se psiloteres sixnotites(ana oktava diplasiazetai i
% suxnotita tou simatos) alla me mikrotero platos 
% - An meiwthei to megethos tou parathirou se 160, den mporoume na analusoume
% tis suxnotites tou deigmatos to idio efkola
% - Stin prwti periptosi pou exoume megalutero mikos parathirou, iparxei
% kalyteri analusi(perissotera deigmata ana segment) kai etsi to ovelap
% epidra se mikrotero pososto


figure(9); hold on; title('Time domain'); plot(mySpeech); 
xlabel('n samples'); hold off;

% window1
figure(10); hold on; title('SEGMENT SIZE = 1600 using hamming window');     
[S, F, T, P] = spectrogram(mySpeech, wind1, noverlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time'); ylabel('Frequency(Hz)'); title('Without filter and window L = 1600');
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

% window2
figure (11); hold on; title('SEGMENT SIZE = 160 using hamming window');
[S, F, T, P] = spectrogram(mySpeech, wind2, noverlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time'); ylabel('Frequency(Hz)'); title('Without filter and window L = 160');
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

%% B2.3 Sound spectrograms with lowpass(IIR)
% - Edw ilopoieitai to xamiloperato filtro
% - Syxnotita deigmatolipsias xrisimopoiw ta 16KHz
%  
% - Meta tin efarmogi toy katwperatou filtrou paratiroume oti ektos
% tou oti kovontai oi psiles sixnotites mikrainei kai to platos tou DFT
% epomenws meiwnetai i entasi. Afto ofeiletai ston metasximatismo Fourier
% afou an efarmosoume ena xamiloperato filtro Nvathmou, oso megalutero
% einai to N, toso mikrotero platos pairnoume
% - Dokimasame proeraitika na enisxusw to sima pollaplasiazontas to sto pedio
% tou xronou me to 2 alla fainetai pws to teliko sima den einai toso katharo

%speechFilt = filter(filtCoe, 1, mySpeech); dialeksa telika to IIR
speechFilt = filter(num, den, mySpeech);

audiowrite('voice_with_filter16KHZ.wav', speechFilt, fsample);
audiowrite('voice_amplified16KHZ.wav', 2 .*speechFilt, fsample);

% window1
figure (12); hold on; title('10ms Hamming with 5ms shift');
[S, F, T, P] = spectrogram(speechFilt, wind1, noverlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time'); ylabel('Frequency(Hz)'); title('With Lowpass and window with L = 1600');
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

% window2
figure (13); hold on; title('10ms Hamming with 5ms shift'); colorbar();
[S, F, T, P] = spectrogram(speechFilt, wind2, noverlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time'); ylabel('Frequency(Hz)'); title('With Lowpass and window L = 160');
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

