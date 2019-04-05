% close all;
% clear all;
% clc;

% Parameters
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
RO = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;
Cn = 0.00001;

% V = [ V1; V2; V3; V4; V5; IL; In]
Vin=1;
G=zeros(6);
C=zeros(6);

%% V1
G(1,:)=[1 0 0 0 0 0]; % V1
C(1,:)=[0 0 0 0 0 0]; % V1

%% V2
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1]; 
C(2,:)=[-C1 +C1 0 0 0 0];

%% V3
G(3,:)=[0 0 1/R3 0 0 -1]; 
C(3,:)=[0 0 Cn 0 0 0]; 

%% V4
G(4,:)=[0 0 -1*alpha/R3 1 0 0]; 
C(4,:)=[0 0 0 0 0 0]; 

%% V5
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0]; 
C(5,:)=[0 0 0 0 0 0];

%% V6
G(6,:)=[0 -1 1 0 0 0]; 
C(6,:)=[0 0 0 0 0 L1]; 


C

%% Q3
% C(2,:)=[-C1 +C1 0 0 0 0 0];
% F=[Vin; 0; 0; 0; 0; 0; 0];
% Vinit=G\F;

C(2,:)=[-C1 +C1 0 0 0 0];
timesteps=1000;
fulltime=1;
stepsize=fulltime/timesteps;
fs_Q2=1/stepsize;
n_Q2=timesteps+1;
fD_Q2=(-(n_Q2-1)/2:(n_Q2-1)/2)*(fs_Q2/n_Q2);
% initial
Data_Q2=zeros(3,3,timesteps+1);
Vold=[1; 0; 0; 0; 0; 0];

input_Q3=zeros(3,timesteps+1);
f_B=1/0.03;
for inputtype=1:3
    for ii=1:timesteps
        input_Q3(inputtype,ii)=exp(-(ii*stepsize-0.1)^2/(2*0.03^2));
    end
end

for inputtype=1:3

    if(inputtype==1)
        Cn_mod=Cn;
    elseif(inputtype==2)
        Cn_mod=Cn*100;
    else
        Cn_mod=Cn*1000;
    end
    C(3,:)=[0 0 Cn_mod 0 0 0]; 
    Vold=[0; 0; 0; 0; 0; 0];

    for ii=1:timesteps
        
        Vin=input_Q3(inputtype,ii);

        In=0.001*randn();
        F=[Vin; 0; -In; 0; 0; 0];
        Data_Q2(inputtype,1,ii+1)=ii*stepsize;
        Data_Q2(inputtype,2,ii+1)=Vin;
        
        A=C/stepsize+G;
        V=(A)\(C*Vold/stepsize+F);
        Data_Q2(inputtype,3,ii+1)=V(5);

        Vold=V;
    end
    
end
C(3,:)=[0 0 Cn 0 0 0]; 

jj=1;
Data_Q2_1=zeros(3,timesteps+1);
Data_Q2_1(:,:)=Data_Q2(jj,:,:);
figure(11)
hold on;
plot(Data_Q2_1(1,:),Data_Q2_1(2,:));
plot(Data_Q2_1(1,:),Data_Q2_1(3,:));
hold off;
legend('V_i','V_O');
title('Figure 11: Voltages for Gaussian Pulse Input with Noise');
ylabel('Voltage (V)');
xlabel('time (s)');

figure(12)
X=fft(Data_Q2_1(2,:));
Y=fft(Data_Q2_1(3,:));
hold on;
plot(fD_Q2,fftshift(abs(X)));
plot(fD_Q2,fftshift(abs(Y)));
hold off;
legend('V_i','V_O');
title('Figure 12: Frequency Domain for Gaussian Pulse Input with Noise');
ylabel('Magnitude');
xlabel('frequency (Hz)');

jj=1;
Data_Q2_1=zeros(3,timesteps+1);
Data_Q2_1(:,:)=Data_Q2(jj,:,:);
figure(13)
subplot(3,1,1)
hold on;
plot(Data_Q2_1(1,:),Data_Q2_1(2,:));
plot(Data_Q2_1(1,:),Data_Q2_1(3,:));
hold off;
legend('V_i','V_O');
title({'\fontsize{14}Figure 13: Voltages with Noise for Different Cn'; '\fontsize{10}Cn=10\mu F'});
ylabel('Voltage (V)');
xlabel('time (s)');

jj=2;
Data_Q2_2=zeros(3,timesteps+1);
Data_Q2_2(:,:)=Data_Q2(jj,:,:);
subplot(3,1,2)
hold on;
plot(Data_Q2_2(1,:),Data_Q2_2(2,:));
plot(Data_Q2_2(1,:),Data_Q2_2(3,:));
hold off;
legend('V_i','V_O');
title('Cn=1mF');
ylabel('Voltage (V)');
xlabel('time (s)');

jj=3;
Data_Q2_3=zeros(3,timesteps+1);
Data_Q2_3(:,:)=Data_Q2(jj,:,:);
subplot(3,1,3)
hold on;
plot(Data_Q2_3(1,:),Data_Q2_3(2,:));
plot(Data_Q2_3(1,:),Data_Q2_3(3,:));
hold off;
legend('V_i','V_O');
title('Cn=10mF');
ylabel('Voltage (V)');
xlabel('time (s)');


% figure(14)
% subplot(3,1,1)
% X=fft(Data_Q2_1(2,:));
% Y=fft(Data_Q2_1(3,:));
% hold on;
% plot(fftshift(abs(X)));
% plot(fftshift(abs(Y)));
% hold off;
% legend('V_i','V_O');
% title('Figure 8: Frequency Domain for Step Input');
% ylabel('Magnitude');
% xlabel('frequency (Hz)');
% 
% %figure(9)
% subplot(3,1,2)
% X=fft(Data_Q2_2(2,:));
% Y=fft(Data_Q2_2(3,:));
% hold on;
% plot(fftshift(abs(X)));
% plot(fftshift(abs(Y)));
% hold off;
% legend('V_i','V_O');
% title('Figure 9: Frequency Domain for Sinusoidal Input with f=33.3Hz');
% ylabel('Magnitude');
% xlabel('frequency (Hz)');
% 
% %figure(10)
% subplot(3,1,3)
% X=fft(Data_Q2_3(2,:));
% Y=fft(Data_Q2_3(3,:));
% hold on;
% plot(fftshift(abs(X)));
% plot(fftshift(abs(Y)));
% hold off;
% legend('V_i','V_O');
% title('Figure 10: Frequency Domain for Gaussian Pulse Input');
% ylabel('Magnitude');
% xlabel('frequency (Hz)');



%% Q3d Change Step
% C(2,:)=[-C1 +C1 0 0 0 0 0];
% F=[Vin; 0; 0; 0; 0; 0; 0];
% Vinit=G\F;

C(2,:)=[-C1 +C1 0 0 0 0];
timesteps_ALT=10000;
fulltime=1;
stepsize_ALT=fulltime/timesteps_ALT;
fs_Q3=1/stepsize_ALT;
n_Q3=timesteps_ALT+1;
fD_Q3=(-(n_Q3-1)/2:(n_Q3-1)/2)*(fs_Q3/n_Q3);
% initial
Data_Q2_ALT=zeros(3,3,timesteps_ALT+1);
Vold=[1; 0; 0; 0; 0; 0];

input_Q3_ALT=zeros(3,timesteps_ALT+1);
f_B=1/0.03;
for inputtype=1:1
    for ii=1:timesteps_ALT
        input_Q3_ALT(inputtype,ii)=exp(-(ii*stepsize_ALT-0.1)^2/(2*0.03^2));
    end
end

for inputtype=1:1
    C(3,:)=[0 0 Cn 0 0 0]; 

    Vold=[0; 0; 0; 0; 0; 0];

    for ii=1:timesteps_ALT
        
        Vin=input_Q3_ALT(inputtype,ii);

        In=0.001*randn();
        F=[Vin; 0; -In; 0; 0; 0];
        Data_Q2_ALT(inputtype,1,ii+1)=ii*stepsize_ALT;
        Data_Q2_ALT(inputtype,2,ii+1)=Vin;
        
        A=C/stepsize_ALT+G;
        V=(A)\(C*Vold/stepsize_ALT+F);
        Data_Q2_ALT(inputtype,3,ii+1)=V(5);

        Vold=V;
    end
    
end
C(3,:)=[0 0 Cn 0 0 0]; 


figure(15)
subplot(2,1,1)
jj=1;
Data_Q2_1=zeros(3,timesteps+1);
Data_Q2_1(:,:)=Data_Q2(jj,:,:);
hold on;
plot(Data_Q2_1(1,:),Data_Q2_1(2,:));
plot(Data_Q2_1(1,:),Data_Q2_1(3,:));
hold off;
legend('V_i','V_O');
title({'\fontsize{14}Figure 14: Voltages with Noise for Different Timesteps'; '\fontsize{10}1 000 timesteps'});
%title('Figure 14: Voltages for 1 000 timesteps');
ylabel('Voltage (V)');
xlabel('time (s)');

subplot(2,1,2)
jj=1;
Data_Q2_1_ALT=zeros(3,timesteps_ALT+1);
Data_Q2_1_ALT(:,:)=Data_Q2_ALT(jj,:,:);
hold on;
plot(Data_Q2_1_ALT(1,:),Data_Q2_1_ALT(2,:));
plot(Data_Q2_1_ALT(1,:),Data_Q2_1_ALT(3,:));
hold off;
legend('V_i','V_O');
title('10 000 timesteps');
ylabel('Voltage (V)');
xlabel('time (s)');

% figure(16)
% subplot(2,1,1)
% X=fft(Data_Q2_1(2,:));
% Y=fft(Data_Q2_1(3,:));
% hold on;
% plot(fD_Q2,fftshift(abs(X)));
% plot(fD_Q2,fftshift(abs(Y)));
% hold off;
% legend('V_i','V_O');
% title('Figure 8: Frequency Domain for Step Input');
% ylabel('Magnitude');
% xlabel('frequency (Hz)');
% 
% subplot(2,1,2)
% X_ALT=fft(Data_Q2_1_ALT(2,:));
% Y_ALT=fft(Data_Q2_1_ALT(3,:));
% hold on;
% % plot(fftshift(abs(X_ALT)));
% % plot(fftshift(abs(Y_ALT)));
% plot(fD_Q3,fftshift(abs(X_ALT)));
% plot(fD_Q3,fftshift(abs(Y_ALT)));
% hold off;
% legend('V_i','V_O');
% title('Figure 8: Frequency Domain for Step Input');
% ylabel('Magnitude');
% xlabel('frequency (Hz)');






