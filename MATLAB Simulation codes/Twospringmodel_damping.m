clear all
close all
clc

pars.mdroplet=1;
pars.kdroplet=1;
zeta=0.75;
pars.c=2*zeta;

pars.fo=sqrt(pars.kdroplet/pars.mdroplet)/(2*3.14);

x0stylus=-1;
x0droplet=-1;

F=[];
Alpha=[];
Lambda=[];
tbyT=[];
Tmax_vel=[];
Tmax_comp=[];
FoF=[];
ACC_drop=[];
ACC_stylus=[];
Fmin=[];
i=0.00001;
for i=0.0001:0.0025:10
    i
pars.mstylus=2000*pars.mdroplet; %%pars.mstylus=0.2*pars.mdroplet for SLF, pars.mstylus=2000*pars.mdroplet for sharpshooter; 
pars.kstylus=2000*pars.kdroplet/(i); %%pars.kstylus=10*pars.kdroplet/(i) for SLF, pars.kstylus=2000*pars.kdroplet/(i) for sharpshooter;
pars.fstylus=sqrt(pars.kstylus/(pars.mdroplet+pars.mstylus))/(2*3.14); %%denominaotr: droplet+stylus mass

fmin=min(pars.fo,pars.fstylus);
if fmin==pars.fo
Fmin=[Fmin 1];
else
Fmin=[Fmin 2];
end
Tmin=1/fmin;%Get the period
numcycles=0.5;
time=linspace(0,numcycles*Tmin,30000);

%% Simulation
BC=[x0stylus 0 x0droplet 0];                                               %input displacement or initial displacement of the web
[t,dynamics] = ode45(@supertwospring,time,BC,[],pars);
dt=t(2)-t(1);
% Kinematics from Simulation
xstylus = dynamics(:,1);%Stylus displacement
vstylus = dynamics(:,2);%Stylus velocity
xdroplet = dynamics(:,3);%Droplet displacement
vdroplet = dynamics(:,4);%Droplet velocity

x_dtos=xdroplet-xstylus;
v_dtos=vdroplet-vstylus;

% Extract details
MaxAbsTimeLoc=find(vdroplet==max(vdroplet));%Corresponds to take-off (When abs velocity is max)
tej=t(MaxAbsTimeLoc);
tejbyT=tej*pars.fstylus;
% Params at Max Abs Velocity 
alpha=(max(vdroplet)/max(vstylus))^2;

Max_vel_loc=find(vstylus==max(vstylus));
tmax_vel=t(Max_vel_loc(1));

[~,Max_comp_loc]=findpeaks(-x_dtos);
tmax_comp=t(Max_comp_loc(1))

%% Geometrical Analysis - Compression Extension cycle
% Acc_drop= gradient(vdroplet,dt);
% Acc_stylus= gradient(vstylus,dt);
% 
% Acc_drop_eje=Acc_drop(MaxAbsTimeLoc)/max(Acc_stylus);
% Acc_stylus_eje=Acc_drop(MaxAbsTimeLoc)/max(Acc_stylus);
% 
% ACC_drop=[ACC_drop_eje Acc_drop_eje];
% ACC_stylus=[ACC_stylus_eje Acc_stylus_eje];




F=[F pars.fstylus];



Alpha=[Alpha alpha];
Lambda=[Lambda sqrt(alpha)];
tbyT=[tbyT tejbyT];
Tmax_vel=[Tmax_vel tmax_vel];
Tmax_comp=[Tmax_comp tmax_comp];

fof=pars.fo/pars.fstylus;

FoF=[FoF fof];

end


%%
figure(1)
subplot(2,1,1)

plot(t,xstylus,'linewidth',3,'color',[0 0 0])
hold on
plot(t,xdroplet,'linewidth',3,'color',[0 0 1])
set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor
xlabel('Time(s)')
ylabel('${x}$(m)', 'Interpreter','latex','fontweight','bold')
title('Absolute Kinematics')

subplot(2,1,2)

plot(t,vstylus,'linewidth',3,'color',[0 0 0])
hold on
plot(t,vdroplet,'linewidth',3,'color',[0 0 1])
set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor
xlabel('Time(s)')
ylabel('$\dot{x}$(m/s)', 'Interpreter','latex','fontweight','bold')




legend

%%%%%%%%%%% Relative Motion
figure(2)
subplot(2,1,1)
plot(t,xstylus,'linewidth',3,'color',[0 0 0])
hold on
plot(t,x_dtos,'linewidth',3,'color',[1 0 0])

set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor
xlabel('Time(s)')
ylabel('${x}$(m)', 'Interpreter','latex','fontweight','bold')
title('Relative Kinematics')

subplot(2,1,2)
plot(t,vstylus,'linewidth',3,'color',[0 0 0])
hold on
plot(t,v_dtos,'linewidth',3,'color',[1 0 0])

set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor
xlabel('Time(s)')
ylabel('$\dot{x}$(m/s)', 'Interpreter','latex','fontweight','bold')



%Phase plots
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

subplot(2,1,1)
plot(xdroplet,vdroplet,'color',[0 0 1],'linewidth',3)
xlabel('${x}$(m) - drop ', 'Interpreter','latex','fontweight','bold')
ylabel('$\dot{x}$ (m/$s$) - droplet', 'Interpreter','latex','fontweight','bold')
set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor

subplot(2,1,2)
plot(xstylus,vstylus,'color',[1 0 0],'linewidth',3)
xlabel('${x}$(m) - stylus', 'Interpreter','latex','fontweight','bold')
ylabel('$\dot{x}$ (m/$s$) - sylus', 'Interpreter','latex','fontweight','bold')
set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor


figure % figure 4_ (fo/f vs. T_max,comp-T_max,vel)*F absolute value
delt_0=(Tmax_comp-Tmax_vel).*F; %delt_0: magnitude of dimensionless time difference between peak ang velocity and max droplet compression moments
plot(2*FoF,abs(delt_0),'color',[1 0 0],'linewidth',3)
hold on
ylim([0 1])
xlabel('$f_o/f$(m)', 'Interpreter','latex','fontweight','bold')
ylabel('$\dot{x}$ (m/$s$) - sylus', 'Interpreter','latex','fontweight','bold')
set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor

figure % figure 5_ fo/f vs. lambda(=velocity_max,droplet/velocity_max,stylus) 

plot(2*FoF,Lambda,'color',[0 0 0],'linewidth',3)
xlim([0, max(2*FoF)])
hold on

xlabel('$f_o/f$(m)', 'Interpreter','latex','fontweight','bold')
ylabel('$\dot{x}$ (m/$s$) - sylus', 'Interpreter','latex','fontweight','bold')
set(gca,'linewidth',1,'Fontsize',16)
grid on ;grid minor

% figure 5 데이터 저장
x = 2*FoF(:);     % x 데이터
y = Lambda(:);    % y 데이터

T = table(x, y);
writetable(T, 'figure5_data_mstylus=200mdroplet,kstylus=200kdroplet.csv');   % save as a CSV file


%%
function dynamics=supertwospring(t,bc,pars)

%Initial condition
xstyl=bc(1);
vstyl=bc(2);
xdroplet=bc(3);
vdroplet=bc(4);

%Inertia
mstylus=pars.mstylus;
mdroplet=pars.mdroplet;

% Forces %
% Spring
Fk_stylus= pars.kstylus*(xstyl);
Fk_droplet= pars.kdroplet*(xdroplet-xstyl);

% Damping 

Fd_droplet=pars.c*(vdroplet-vstyl);



%% Equation
vstylus=vstyl;
astylus=(-Fk_stylus-Fk_droplet-Fd_droplet)/pars.mstylus;
vdrop=vdroplet;
adrop=(-Fk_droplet-Fd_droplet)/pars.mdroplet;

%% Output
dynamics=[vstylus;astylus;vdrop;adrop];
end