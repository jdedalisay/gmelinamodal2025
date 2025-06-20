%% FRF analyzer

clear; clc;

load m_ih.mat
load m_a1.mat
load m_a2.mat

%% shifting and removing stuff
for idx = 1:24
    ref = m_ih(:,idx);
    [~, jdx] = max(ref);
    while ref(jdx) > 0.1 % threshold
        jdx = jdx - 1;
    end
    jdx = jdx - 1;
    ref = ref(jdx:end);
    m_a1d(:,idx) = [m_a1(jdx:end,idx); zeros(64000-length(ref),1)];
    m_a2d(:,idx) = [m_a2(jdx:end,idx); zeros(64000-length(ref),1)];
    m_ihd(:,idx) = [ref; zeros(64000-length(ref),1)];
end
m_ihd(30:end,:) = 0; % cleaning input signal
m_a1d(27800:end,:) = 0; % removing some weird vibration at some point
m_a2d(27800:end,:) = 0; % removing some weird vibration at some point

% getting displacement data from acceleration
fs = 6400;
for idx = 1:24
    m_v1d(:,idx) = cumtrapz(m_a1d(:,idx))/fs;
    m_x1d(:,idx) = cumtrapz(m_v1d(:,idx))/fs;
    m_v2d(:,idx) = cumtrapz(m_a2d(:,idx))/fs;
    m_x2d(:,idx) = cumtrapz(m_v2d(:,idx))/fs;
end

% filtering acceleration data
fnyq = 3200;
[b, a] = butter(6,35/fnyq,'high');
for idx = 1:size(m_a1d,2)
    atemp = [-flip(m_a1d(:,idx));(m_a1d(:,idx))];
    vtemp = [-flip(m_v1d(:,idx));(m_v1d(:,idx))];
    xtemp = [-flip(m_x1d(:,idx));(m_x1d(:,idx))];
    atemp = filtfilt(b,a,atemp);
    vtemp = filtfilt(b,a,vtemp);
    xtemp = filtfilt(b,a,xtemp);
    m_a1df(:,idx) = atemp(end-64000+1:end);
    m_v1df(:,idx) = vtemp(end-64000+1:end);
    m_x1df(:,idx) = xtemp(end-64000+1:end);
end

for idx = 1:size(m_a2d,2)
    atemp = [-flip(m_a2d(:,idx));(m_a2d(:,idx))];
    vtemp = [-flip(m_v2d(:,idx));(m_v2d(:,idx))];
    xtemp = [-flip(m_x2d(:,idx));(m_x2d(:,idx))];
    atemp = filtfilt(b,a,atemp);
    vtemp = filtfilt(b,a,vtemp);
    xtemp = filtfilt(b,a,xtemp);
    m_a2df(:,idx) = atemp(end-64000+1:end);
    m_v2df(:,idx) = vtemp(end-64000+1:end);
    m_x2df(:,idx) = xtemp(end-64000+1:end);
end

sub_time = 10; 
nSamples = size(m_ihd,1); 
fs = nSamples / sub_time; 
fnyq = fs/2;

totalTrials = size(m_ihd,2); 
nDOFs = 8; 
trials = totalTrials/nDOFs;

m_a1_rs = reshape(m_a1df,[nSamples*trials,nDOFs]);
m_v1_rs = reshape(m_v1df,[nSamples*trials,nDOFs]);
m_x1_rs = reshape(m_x1df,[nSamples*trials,nDOFs]);
m_a2_rs = reshape(m_a2df,[nSamples*trials,nDOFs]);
m_v2_rs = reshape(m_v2df,[nSamples*trials,nDOFs]);
m_x2_rs = reshape(m_x2df,[nSamples*trials,nDOFs]);
m_ih_rs = reshape(m_ihd,[nSamples*trials,nDOFs]);

%% testing matching, the excitation and response must align-ish

idx = 1:8; % test idx
clf;
nexttile
plot(m_a1_rs(:,idx));
hold on;
plot(m_ih_rs(:,idx));
nexttile
plot(m_a2_rs(:,idx));
hold on;
plot(m_ih_rs(:,idx));

%% modalfrf
clf;
[frf, freqs] = modalfrf(m_ih_rs,m_x1_rs,fs,nSamples,'Measurement','rovinginput','Sensor','dis');
[temp, ~] = modalfrf(m_ih_rs,m_x2_rs,fs,nSamples,'Measurement','rovinginput','Sensor','dis');
frf(:,2,:) = temp;
for idx = 1:8
semilogy(freqs,abs(frf(:,:,idx)));
hold on;
end
%xlim([0,100]);

%%
clf;
modalsd(frf,freqs,fs,'MaxModes',50);


%%
clf;
phf = [47.5, 63.4, 70.7, 160.3,377.1,568.8,621.5,...
    1052.8, 1114.9,1255.8,1576.8,1821.4,1981.3]; %ask saang plot nakuha
modalfit(frf,freqs,fs,50,'PhysFreq',phf,'DriveIndex',[1,4],'FitMethod','lsce','FreqRange',[0,2000]);
[~,~,ms] = modalfit(frf,freqs,fs,50,'PhysFreq',phf,'DriveIndex',[1,4],'FitMethod','lsce','FreqRange',[0,2000]);
%%
idx = 1;

plot(abs(ms(:,idx)).*sign(real(ms(:,idx))));
%%
clf;
FRF = frf;
f = freqs;
nodeids = 1:8;
fmax = 100;
FRF2 = permute(FRF,[1 3 2]);
newplot;
waterfall(f(f<fmax),nodeids,imag(FRF2(f<fmax,nodeids))');
yticks(nodeids);
xlabel('Frequency (Hz)')
ylabel('Location')
zlabel('Imaginary Part of FRF') 
ax = gca;
set(ax,'YDir','Reverse');
zlim([-1,1]*5e-5);
xlim([10,fmax]);

%% ms through imaginary
clear ms
testfr = [45.2,62.1, 156.2, 366.5, 536.8, 610.4, 1025.1, 1079.1, 1219.6, 1540];
for jdx = 1:10
    fr = testfr(jdx)
    [~,idx] = min(abs(freqs-fr));
    resp = frf(idx,1,:);
    resp = reshape(resp,1,[]);
    resp
    ms(jdx,:) = imag(resp);
    %ms
end
clf
plot(ms')
legend
ms = ms';
%% MAC
clc; 
close all;
m = zeros(8,10);
m1 = readmatrix('m1.txt');
m2 = readmatrix('m2.txt');
m3 = readmatrix('m3.txt');
m4 = readmatrix('m4.txt');
m5 = readmatrix('m5.txt');
m6 = readmatrix('m6.txt');
m7 = readmatrix('m7.txt');
m8 = readmatrix('m8.txt');
m9 = readmatrix('m9.txt');
m10 = readmatrix('m10.txt');
m(:,1) = m1(:,2);
m(:,2) = m2(:,2);
m(:,3) = m3(:,2);
m(:,4) = m4(:,2);
m(:,5) = m5(:,2);
m(:,6) = m6(:,2);
m(:,7) = m7(:,2);
m(:,8) = m8(:,2);
m(:,9) = m9(:,2);
m(:,10) = m10(:,2);
% v = real(ms(:,1:10));
v = ms;
m_mode = 1;
v_mode = 1;

nexttile;
stem(m(:,m_mode)/norm(m(:,m_mode)));
hold on;
stem(v(:,v_mode)/norm(v(:,v_mode)));
m = m/norm(m);
v = v/norm(v);
mac = ((abs(m(:,m_mode)'*v(:,v_mode)))^2)/((m(:,m_mode)'*m(:,m_mode))*(v(:,v_mode)'*v(:,v_mode)));
disp(['Ansys Mode: ', num2str(m_mode)])
disp(['Experiment Mode: ', num2str(v_mode)])
disp(['MAC: ', num2str(mac)])
legend(['Ansys Mode: ', num2str(m_mode)],['Experiment Mode: ', num2str(v_mode)]);

[num_dofs_m, num_modes_m] = size(m);
[num_dofs_v, num_modes_v] = size(v);

if num_dofs_m ~= 8
    error('The simulation data (m) must have 8 DOFs (rows).');
end

if num_dofs_v ~= 8
    error('The experimental data (v) must have 8 DOFs (rows).');
end

if num_modes_m ~= 10
    error('The simulation data (m) must have 10 modes (columns).');
end

if num_modes_v ~= 10
    error('The experimental data (v) must have 10 modes (columns).');
end

num_modes = 10; % Number of modes
num_dofs = 8;   % Number of DOFs
MAC = zeros(num_modes, num_modes);

% Calculate MAC for each mode pair
for i = 1:num_modes
    for j = 1:num_modes
        % Extract mode shapes (columns from the matrices)
        phi_exp = v(:, i); % Experimental mode shape i (column vector)
        phi_sim = m(:, j); % Simulation mode shape j (column vector)
        MAC(i, j) = abs(phi_exp' * phi_sim)^2 / ((phi_exp' * phi_exp) * (phi_sim' * phi_sim));
    end
end

disp('MAC Matrix:');
disp(MAC);

% Visualize MAC matrix with colormap
nexttile;
imagesc(MAC);
colorbar;
title('Modal Assurance Criterion (MAC)');
xlabel('Ansys Modes');
ylabel('Experiment Modes');
xticks(1:num_modes);
yticks(1:num_modes);
colormap jet; % Use the 'jet' colormap (or any other colormap you prefer)

% Add text annotations to the MAC matrix for better readability (optional)
for i = 1:num_modes
    for j = 1:num_modes
        text(j, i, sprintf('%.2f', MAC(i, j)), 'HorizontalAlignment', 'center', 'Color', 'white');
    end
end

