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

% V = [ V1; V2; V3; V4; V5; IL]
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
C(3,:)=[0 0 0 0 0 0]; 

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
G


%% Q1
%% Part b
numsteps=20;
Data=zeros(3,numsteps);
Data(1,:)=linspace(-10,10,numsteps);
for ii=1:numsteps
    % omega=0 for DC
    F=[Data(1,ii); 0; 0; 0; 0; 0];
    V=G\F;
    Data(2,ii)=V(5);
    Data(3,ii)=V(3);
end

figure(1)
hold on;
plot(Data(1,:),Data(2,:));
plot(Data(1,:),Data(3,:));
hold off;
legend('V_O','V_3');
title('Figure 1: Voltage as a Function of Input Voltage');
ylabel('Voltage (V)');
xlabel('V_i (V)');


%% Part c
numsteps=2000;
Data2=zeros(2,numsteps);
Data2(1,:)=linspace(0,500,numsteps);
Vin=1;
for ii=1:numsteps
    % omega=0 for DC
    omega=Data2(1,ii);
    F=[Vin; 0; 0; 0; 0; 0];
    V=(G+1j*omega*C)\F;
    Data2(2,ii)=V(5);
    %Data2(3,ii)=Data2(2,ii)/Data2(1,ii);
end

figure(2)
plot(Data2(1,:),real(Data2(2,:)));
title('Figure 2: Output Voltage as a Function of Angular Frequency');
legend('V_O')
ylabel('V_O (V)');
xlabel('\omega (radians/s)');

figure(3)
semilogx(Data2(1,:),20*log10(real(Data2(2,:))./Vin));
title('Figure 3: Gain as a Function of Angular Frequency');
legend('dB(V_O/V_i)');
ylabel('Gain (dB)')
xlabel('\omega (radians/s)');


%% Part d
Vin=1;
omega=pi;
count=10000;
Data3=zeros(1,count);
for ii=1:count
    C2=C1+randn()*0.05;
    C(2,:)=[-C2 +C2 0 0 0 0];
    F=[Vin; 0; 0; 0; 0; 0];
    V=(G+1j*omega*C)\F;
    Data3(1,ii)=V(5);
end

figure(4);
hist(real(Data3),50);
title('Figure 4: Gain for perturbations in C');
ylabel('Count')
xlabel('Gain');


%% Q2
C(2,:)=[-C1 +C1 0 0 0 0];
timesteps=1000;
fulltime=1;
stepsize=fulltime/timesteps;
fs_Q2=1/stepsize;
n_Q2=timesteps+1;
fD_Q2=(-(n_Q2-1)/2:(n_Q2-1)/2)*(fs_Q2/n_Q2);
% initial
Data_Q2=zeros(3,3,timesteps+1);
Vold=[0; 0; 0; 0; 0; 0];

input_Q2=zeros(3,timesteps+1);
f_B=1/0.03;
for inputtype=1:3
    for ii=1:timesteps
        if(inputtype==1)
            if(ii*stepsize<0.03)
                input_Q2(inputtype,ii)=0;
            else
                input_Q2(inputtype,ii)=1;
            end
        elseif(inputtype==2)
            input_Q2(inputtype,ii)=sin(2*pi*f_B*ii*stepsize);
        else
            input_Q2(inputtype,ii)=exp(-(ii*stepsize-0.1)^2/(2*0.03^2));
        end
    end
end

for inputtype=1:3
    Vold=[0; 0; 0; 0; 0; 0];
    for ii=1:timesteps
        
        Vin=input_Q2(inputtype,ii);
        F=[Vin; 0; 0; 0; 0; 0];
        Data_Q2(inputtype,1,ii+1)=ii*stepsize;
        Data_Q2(inputtype,2,ii+1)=Vin;
        
        A=C/stepsize+G;
        V=(A)\(C*Vold/stepsize+F);
        Data_Q2(inputtype,3,ii+1)=V(5);

        Vold=V;
    end
    
end

jj=1;
Data_Q2_1=zeros(3,timesteps+1);
Data_Q2_1(:,:)=Data_Q2(jj,:,:);
figure(5)
hold on;
plot(Data_Q2_1(1,:),Data_Q2_1(2,:));
plot(Data_Q2_1(1,:),Data_Q2_1(3,:));
hold off;
legend('V_i','V_O');
title('Figure 5: Voltages for Step Input');
ylabel('Voltage (V)');
xlabel('time (s)');


jj=2;
Data_Q2_2=zeros(3,timesteps+1);
Data_Q2_2(:,:)=Data_Q2(jj,:,:);
figure(6)
hold on;
plot(Data_Q2_2(1,:),Data_Q2_2(2,:));
plot(Data_Q2_2(1,:),Data_Q2_2(3,:));
hold off;
legend('V_i','V_O');
title('Figure 6: Voltages for Sinusoidal Input with f=33.3Hz');
ylabel('Voltage (V)');
xlabel('time (s)');


jj=3;
Data_Q2_3=zeros(3,timesteps+1);
Data_Q2_3(:,:)=Data_Q2(jj,:,:);
figure(7)
hold on;
plot(Data_Q2_3(1,:),Data_Q2_3(2,:));
plot(Data_Q2_3(1,:),Data_Q2_3(3,:));
hold off;
legend('V_i','V_O');
title('Figure 7: Voltages for Gaussian Pulse Input');
ylabel('Voltage (V)');
xlabel('time (s)');


figure(8)
X=fft(Data_Q2_1(2,:));
Y=fft(Data_Q2_1(3,:));
hold on;
plot(fD_Q2,fftshift(abs(X)));
plot(fD_Q2,fftshift(abs(Y)));
hold off;
legend('V_i','V_O');
title('Figure 8: Frequency Domain for Step Input');
ylabel('Magnitude');
xlabel('frequency (Hz)');

figure(9)
X=fft(Data_Q2_2(2,:));
Y=fft(Data_Q2_2(3,:));
hold on;
plot(fD_Q2,fftshift(abs(X)));
plot(fD_Q2,fftshift(abs(Y)));
hold off;
legend('V_i','V_O');
title('Figure 9: Frequency Domain for Sinusoidal Input with f=33.3Hz');
ylabel('Magnitude');
xlabel('frequency (Hz)');

figure(10)
X=fft(Data_Q2_3(2,:));
Y=fft(Data_Q2_3(3,:));
hold on;
plot(fD_Q2,fftshift(abs(X)));
plot(fD_Q2,fftshift(abs(Y)));
hold off;
legend('V_i','V_O');
title('Figure 10: Frequency Domain for Gaussian Pulse Input');
ylabel('Magnitude');
xlabel('frequency (Hz)');



