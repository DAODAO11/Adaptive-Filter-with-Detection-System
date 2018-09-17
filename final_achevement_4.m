clear all;
[x,Fs]=audioread('E:\recording\2018.03.02\2018.03.02.2nd.1.wav');
x1=x(:,1);   

L = length(x1);       % Length of signal
T = 1/Fs;             % Sampling period
t = (0:L-1)*T;        % Time vector
sound(x1,Fs);
figure(1);
stem(t,x1);
title('time plot')
xlabel('t (seconds)')
ylabel('amplitude')

Y = fft(x1);
P2 = abs(Y);
P1 = P2(1:(L/2));
P1(1:end) = 2*P1(1:end);
f = Fs*(0:(L/2-1))/L;
figure(2);
stem(f,P1);
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%short-time Fourier transform
s = spectrogram(x1);
figure(3)
spectrogram(x1,Fs,'yaxis')
% spectrogram(x1,'yaxis')
title('spectrogram')

nsc = floor(L/4.5);
nov = floor(nsc/2);
S=nextpow2(nsc);
nff = max(256,2^nextpow2(nsc));
D = spectrogram(x1,hamming(nsc),nov,nff);
maxerr = max(abs(abs(D(:))-abs(s(:))));

ns = 8;
ov = 0.5;
lsc = floor(L/(ns-(ns-1)*ov));
B = spectrogram(x1,lsc,floor(ov*lsc),nff);
C = max(abs(abs(B(:))-abs(s(:))));

% Hamming window 1
figure(4)
spectrogram(x1,hamming(512),500,512,Fs,'yaxis');
% set(gca,'xdir','reverse')
title('short time fourier transform')

% Hamming window 2
figure(5)
spectrogram(x1,1024,1000,1024,Fs,'yaxis')
detection = spectrogram(x1,1024,1000,1024,Fs);
title('short time fourier transform')

% Logarithmic
figure(6)
spectrogram(x1,1024,1000,[],Fs,'yaxis')
title('short time fourier transform with Logarithmic frequency scale')
ax = gca;
ax.YScale = 'log';

% Blackman window
figure(7)
spectrogram(x1,blackman(128),60,128,Fs,'yaxis')
xlabel('time')
ylabel('frequency (Hz)')
title('Blackman window')
ax = gca;
ax.YDir = 'reverse';

detection1=abs(detection);
[m,n]=size(detection1);
q=8000*2*m/Fs;
q=floor(q);
p=15000*2*m/Fs;
p=floor(p);
su1= rand(n,1);
for i=1:1:n
for j=q:1:p
    if j==q
     su1(i,1)=detection1(j,i);
    else
     su1(i,1)= su1(i)+detection1(j,i); 
    end
end
end



wlen = 1024;
% hop = wlen/4;
hop=wlen/4;
nfft = wlen;
[stft, fvector, t_stft] =stft(x1, wlen, hop, nfft, Fs);
stfts=size(stft);
stfts2=stfts(1,2);
stfts1=stfts(1,1);

u=(L-1000)/24;
u=floor(u);


det1 = real(stft);
det2 = imag(stft);
% det=det';
% signal = x1;
% signal=signal';
% stft2=(L-1000)/24;
% stft2=floor(stft2);
d=floor(5000*2/Fs*stfts1);
signal=x1;


w=1;
i=1;
while(w<n)
 
    
% for i=1:n
while(i<n)    
  Su= su1(i,1);
  if Su>20
      K=i;
      for j=i:1:n
      Su= su1(j,1);
      if Su<20
      if su1(j+10,1)<20
      J=j;
           
% K1=K*L/u;
% K1=round(K1);
% K1=K1-0.05*Fs;
% J1=J*L/u;
% J1=round(J1);
% J1=J1+0.05*Fs;
K1=floor(K*stfts2/u);
K1=K1-floor(0.01*stfts2*Fs/L);
J1=floor(J*stfts2/u);
J1=J1+floor(0.08*stfts2*Fs/L);


% y=fir1(48,0.22,'low');
% F=conv(signal(K1:J1),y);
% signal(K1:J1)=F(1:(J1-K1+1));

for b=1:1:d
      
% signal = con2seq(signal);    % Convert concurrent vectors to sequential vectors
% Xi = signal(K1-50:K1-1);             % set delay for input signal
% X = signal(K1:J1);              % X is the input signal
% timex = t(K1:J1);
% T = signal(K1:J1);              % T is the target,which should hav the same size as input 

% G=det1(b,:);
A=det1(b,:);
% A=A';
A=con2seq(A);  
Xi=A(2*K1-J1-11:2*K1-2-J1);
X = A(2*K1-1-J1:K1-1); 
% timex = t(K1:J1);
T = A(K1:J1);


net = linearlayer(1:10,0.001);
% view(net)
[net1,Y1] = adapt(net,X,T,Xi);   % y is the output of the system
 
% net = linearlayer(1:50,0.2);
% view(net)
% [net,Y1] = adapt(net,Y,Y,Xi);   % y is the output of the system

% signal = seq2con(signal);
A = cell2mat(A);
Y2 = cell2mat(Y1);
A(K1:J1)=Y2;
det1(b,:)=A;

% g=t(1:J1-K1+1);
% h=G(K1:J1);
% figure(8);
% % ,g,Y2,'r'
% plot(g,h,'g',g,Y2,'r')
% xlabel('Time');
% ylabel('magnitude');
% legend('Target','Output');
% title('Output and Target Signals');



A1=det2(b,:);
A1=con2seq(A1);  
Xi1=A1(2*K1-J1-11:2*K1-2-J1);
X1 = A1(2*K1-1-J1:K1-1); 
T1 = A1(K1:J1);

net = linearlayer(1:10,0.001);
[net2,Y11] = adapt(net,X1,T1,Xi1);   % y is the output of the system

A1 = cell2mat(A1);
Y21 = cell2mat(Y11);
A1(K1:J1)=Y21;
det2(b,:)=A1;

end

for b=d+1:1:stfts1
      
% signal = con2seq(signal);    % Convert concurrent vectors to sequential vectors
% Xi = signal(K1-50:K1-1);             % set delay for input signal
% X = signal(K1:J1);              % X is the input signal
% timex = t(K1:J1);
% T = signal(K1:J1);              % T is the target,which should hav the same size as input 

% G=det1(b,:);
A=det1(b,:);
% A=A';
A=con2seq(A);  
Xi=A(2*K1-J1-11:2*K1-2-J1);
X = A(2*K1-1-J1:K1-1); 
% timex = t(K1:J1);
T = A(K1:J1);


net = linearlayer(1:10,0.00001);
% view(net)
[net3,Y1] = adapt(net,X,T,Xi);   % y is the output of the system
 
% net = linearlayer(1:50,0.2);
% view(net)
% [net,Y1] = adapt(net,Y,Y,Xi);   % y is the output of the system

% signal = seq2con(signal);
A = cell2mat(A);
Y2 = cell2mat(Y1);
A(K1:J1)=Y2;
det1(b,:)=A;

% g=t(1:J1-K1+1);
% h=G(K1:J1);
% figure(8);
% % ,g,Y2,'r'
% plot(g,h,'g',g,Y2,'r')
% xlabel('Time');
% ylabel('magnitude');
% legend('Target','Output');
% title('Output and Target Signals');

A1=det2(b,:);
A1=con2seq(A1);  
Xi1=A1(2*K1-J1-11:2*K1-2-J1);
X1 = A1(2*K1-1-J1:K1-1); 
T1 = A1(K1:J1);

net = linearlayer(1:10,0.00001);
[net4,Y11] = adapt(net,X1,T1,Xi1);   % y is the output of the system

A1 = cell2mat(A1);
Y21 = cell2mat(Y11);
A1(K1:J1)=Y21;
det2(b,:)=A1;

end

det3=det1+det2*1i;

det11=det3(:,K1:J1);
% [x2,t2] = istft(det11,1024,24,Fs);
[x_istft, t_istft] = istft(det11, wlen, hop, nfft, Fs);
x21=x_istft(1,385:length(x_istft)-384);
x21=real(x21);
x21=2*x21;
K2=K1*L/stfts2;
K2=round(K2);
J2=J1*L/stfts2;
J2=floor(J2);
signal(K2:J2-256)=x21(1:J2-K2-255);

% K1=K*L/u;
% K1=round(K1);
% K1=K1-0.05*Fs;
% J1=J*L/u;
% J1=round(J1);

w=J+1;
i=J;      
      break
      else
      J=0;
      end
      else 
      J=0 ;   
      end
      end
       break;
  else
      K=0;
  end
i=i+1;
end


w=w+1;
end

% signal=signal';
% x1=x1'; 
% signal=abs(signal);

figure(9)
plot(t,x1,'g',t,signal,'r')
xlabel('Time');
ylabel('magnitude');
legend('Target','Output');
title('Output and Target Signals');


figure(10)
E = x1-signal;
plot(t,E,'r')
hold off
xlabel('Time');
ylabel('magnitude');
title('Error Signal');

% % Y1 = seq2con(Y1);
% % signal = seq2con(signal);
% % signal_new = cell2mat(signal);
% signal_new =signal;
% Y_new = cell2mat(Y1);
% signal_new(K1:J1)=Y_new;
% signal_new=signal_new';
% sound(signal_new,Fs);

figure(11);
plot(t,signal)
xlabel('Time');
ylabel('magnitude');
title('Actual Output After filtering');

% percentage=x1(242551:277830);
% for i=1:1:35280
% percentage(i)=E(i)/x1(242550+i);
% end
% figure(12);
% stem(timex,percentage)

sound(signal,Fs);
signal=signal';
figure(12)
spectrogram(signal,1024,1000,1024,Fs,'yaxis')
detection2 = spectrogram(signal,1024,1000,1024,Fs);