close all
load m_ih.mat
load m_a1.mat
load m_a2.mat
load f_ih.mat
load f_a2.mat
load f_a1.mat

d = daq("ni");

d.Rate = 6400; 

sub_time = 10;

t = (0:1/d.Rate:sub_time)';
t(end,:) = [];
t(end,:) = []; %clarify w jded

%ih = addinput(d,"Dev1", "ai0", "IEPE");
%a1 = addinput(d,"Dev1", "ai1", "IEPE"); %8694
%a2 = addinput(d,"Dev1","ai2","IEPE"); %8696
ih = addinput(d,"Dev1", "ai0", "Accelerometer");
a1 = addinput(d,"Dev1", "ai1", "Accelerometer"); %8694
a2 = addinput(d,"Dev1","ai2","Accelerometer"); %8696

ih.Sensitivity = 0.01;
a1.Sensitivity = 0.105;
a2.Sensitivity = 0.1017;

%tot_time = 10; %10 seconds
%sub_time = 8; % 8 seconds
%L_offset = 10;%number of samples before the peak to offset

% ih.Sensitivity = .002208;
% a1.Sensitivity = .00105;
% a2.Sensitivity = .00922;                                             
 
m_ih = [];
m_a1 = [];
m_a2 = [];
f_ih = [];
f_a1 = [];
f_a2 = [];


%m_ih(:,1) = 0:1/d.Rate:sub_time; 

close all
 
disp('press any key to continue') 
pause 

start(d,"continuous")
Dev1_1 = [];




data = read(d,seconds(10));
Dev1_1 = [Dev1_1; data];

stop(d)

%[M,I] = max(Dev1_1.Dev1_ai0);

%x1 = I - 20;
%xx1 = x1 + d.Rate * 8;

%finding zero before the peak

%x1 = ;
%xx1 = x1 + d.Rate*8;

d_ih = data.Dev1_ai0(:,:);
d_a1 = data.Dev1_ai1(:,:);
d_a2 = data.Dev1_ai2(:,:);

d_ih(end,:) = [];
d_a1(end,:) = [];
d_a2(end,:) = [];



figure(1);
tiledlayout('vertical')

nexttile
plot(t, d_ih)
title('ih')
xlabel("Time (s)")
ylabel("Amplitude")
hold on

nexttile
plot(t, d_a1)
title('a1')
xlabel("Time (s)")
ylabel("Amplitude")

nexttile
plot(t, d_a2)
title('a2')
xlabel("Time (s)")
ylabel("Amplitude")


hold on

figure(2);


%fs = d.Rate; %6400

fs = d.Rate; 
dt = 1/fs; 
N = length(t);
T = dt*N; 
df = 1/T; 
f = (-N/2:((N/2)-1))*df;

y_ih = d_ih;
yHat_ih = fft(y_ih);
yHat_ih = fftshift(yHat_ih);

y_a1 = d_a1;
yHat_a1 = fft(y_a1);
yHat_a1 = fftshift(yHat_a1);

y_a2 = d_a2;
yHat_a2 = fft(y_a2);
yHat_a2 = fftshift(yHat_a2);


tiledlayout('vertical')

nexttile
semilogy(f, abs(yHat_ih))
title('ih')
xlabel("Frequency (Hz)")
ylabel("Amplitude")
hold on

nexttile
semilogy(f, abs(yHat_a1))
title('a1')
xlabel("Frequency (Hz)")
ylabel("Amplitude")

nexttile
semilogy(f, abs(yHat_a2))
title('a2')
xlabel("Frequency (Hz)")
ylabel("Amplitude")

hold off


post =input("ano gusto mo? ('S','D'): ",'s');

if post == 'S' || post == 's'
    m_ih(:,end + 1) = d_ih;
    m_a1(:,end + 1) = d_a1;
    m_a2(:,end + 1) = d_a2;

    f_ih(:,end + 1) = yHat_ih;
    f_a1(:,end + 1) = yHat_a1;
    f_a2(:,end + 1) = yHat_a2;

    save m_ih.mat m_ih
    save m_a1.mat m_a1
    save m_a2.mat m_a2

    save f_ih.mat f_ih
    save f_a1.mat f_a1
    save f_a2.mat f_a2


elseif post == 'D' || post == 'd'
    save m_ih.mat m_ih
    save m_a1.mat m_a1
    save m_a2.mat m_a2

    save f_ih.mat f_ih
    save f_a1.mat f_a1
    save f_a2.mat f_a2

else
end
clearvars;
