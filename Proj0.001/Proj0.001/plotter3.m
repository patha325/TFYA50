clear;
clc;
s=importdata('L:\TFYA50\Proj0.001\Proj0.001\energytemp.txt');
% DATA LINE 1
%fs2 << steps << " " << step_out << " " <<time_step<< " " << 0<< " " << 0<< " " << 0<< " " << 0<< " " << 0<<endl;

time=1:s(1,3):s(1,2);
time2=1:s(1,3):s(1,2)-s(1,4);
%fs2 << total_energy_vector[i] << " " << pot_energy_vector[i] << " " << kin_energy_vector[i]<< " " << temperature_vector[i]<< " " <<pressure_vector[i]<< " " <<MSD_vector[i]<< " " <<Debye_temp_vector[i]<< " " <<Diff_coeff_vector[i] <<endl;

% total_energy = zeros(size(s,1));
% potential_energy = zeros(size(s,1));
% kinetic_energy = zeros(size(s,1));
% temp = zeros(size(s,1));
% pressure = zeros(size(s,1));
% MSD = zeros(size(s,1) - s(1,4));
% Debye = zeros(size(s,1) - s(1,4));
% Diff = zeros(size(s,1) - s(1,4));
% coh_e = zeros(size(s,1) - s(1,4));

for i=2:size(s,1)
total_energy(i-1)=s(i,1);
potential_energy(i-1)=s(i,2);
kinetic_energy(i-1)=s(i,3);
temp(i-1)=s(i,4);
pressure(i-1)=s(i,5);
%MSD(i-1)=s(i,6);
%Debye(i-1)=s(i,7);
%Diff(i-1)=s(i,8);
%coh_e(i-1)=s(i,9);
end

for i=2:(size(s,1)-s(1,4))
MSD(i-1)=s(i+s(1,4),6);
Debye(i-1)=s(i+s(1,4),7);
Diff(i-1)=s(i+s(1,4),8);
coh_e(i-1)=s(i+s(1,4),9);
end

%%subplot(3,3,1)%Plot Total energy
figure(1);
plot(time,total_energy);
title('Total Energy');
xlabel('Time [fs]')
ylabel('Energy [eV]')
grid on
%%subplot(3,3,2)%Plot Potential energy
figure(2);
plot(time,potential_energy);
title('Potential Energy');
xlabel('Time [fs]')
ylabel('Energy [eV]')
grid on
%%subplot(3,3,3)%Plot Kinetic energy
figure(3);
plot(time,kinetic_energy);
title('Kinetic Energy');
xlabel('Time [fs]')
ylabel('Energy [eV]')
grid on
%%subplot(3,3,4)%Plot Temperature
figure(4);
plot(time,temp);
title('Temperature');
xlabel('Time [fs]')
ylabel('Temp [K]')
grid on
%%subplot(3,3,5)%Plot Pressure
figure(5);
plot(time,pressure);
title('Pressure');
xlabel('Time [fs]')
ylabel('Pressure [eV][fs]^2[Å]^-2')
grid on
%%subplot(3,3,6)%Plot MSD
figure(6);
plot(time2,MSD);
title('MSD');
xlabel('Time [fs]')
ylabel('MSD [Å]^2')
grid on
%%subplot(3,3,7)%Plot Debye
figure(7);
plot(time2,Debye);
title('Debye Temp');
xlabel('Time [fs]')
ylabel('Temp [K]')
grid on
%%subplot(3,3,8)%Plot Diffusion coeff
figure(8);
plot(time2,Diff);
title('Diffusion coefficient');
xlabel('Time [fs]')
ylabel('Diffusion coeff')
grid on
%%subplot(3,3,9)%Plot Cohesive energy
figure(9);
plot(time2,coh_e);
title('Cohesive energy');
xlabel('Time [fs]')
ylabel('Cohesive energy per atom [eV]')
grid on


%Plot Total


