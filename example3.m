% This is an example of how to use the Multi-Pitch Estimation Toolbox for
% pitch estimation on a real-life signal, in this case a speech signal.

% Load file
addpath pitch/
filename=['fa.wav'];
directory='./';
[input,fs,nbits]=wavread(strcat(directory,filename));
input=input(:,1);

% Set parameters
F=2^14;
f_est=[];
t_est=[];
N=2*round(0.03*fs/2);
H=N/2;
pos=1:N;
f0_min=100;
f0_max=400;
w0_lim=2*pi*[f0_min f0_max]/fs;

% Process file
while pos(end)<=length(input),

    x=hilbert(input(pos));

    [w0,L]=joint_orth(x,w0_lim,F,length(x)/2);
    f0=w0/2/pi*fs
    
    f_est=[f_est f0];
    t_est=[t_est pos(1)];

    pos=pos+H;
end

% Plot spectrogram and estimates
figure(1);clf;
subplot(2,1,1);specgram(input,1024,fs,256,128);
set(gca,'FontSize',16)
xlabel('Time [ms]','Fontsize',16);
ylabel('Frequency [Hz]','Fontsize',16)
box on
subplot(2,1,2);h=plot(t_est/fs,f_est,'x');
set(gca,'FontSize',16)
xlabel('Time [ms]','Fontsize',16);
ylabel('Fundamental Frequency [Hz]','Fontsize',16)
axis([0 length(input)/fs f0_min f0_max]);