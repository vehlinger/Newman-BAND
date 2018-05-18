% Victoria Ehlinger
% UC Berkeley Department of Chemical Engineering
%
% test of band code based on Newman's approach
% multicomponent diffusion in a stefan tube
% Taylor & Krishna Example 2.1.1

clear all
close all
clc

n = 3;               % number of species, including the solvent
                     % 1 = acetone, 2 = methanol, 3 = air
nm = n;              % number of species  
n = 2*n-1;           % number of unknowns at each mesh point (including fluxes)
nj = 21;             % number of mesh points
C = zeros(n,nj);     % change variable (mol fraction)
L = 0.238;           % modeling length (m)

% initial conditions
C(1:3,1:nj-1) = repmat([0.319;0.528;0.153],1,nj-1);
C(1:3,nj) = [0,0,1];
C(4:5,:) = repmat([1e-5; 1e-5],1,nj);
% initial fluxes  = 0

P = 99.4; % kPa
T = 328.5; % K
R = 8.314; % (kPa L)/(K mol)

ctot = P/(R*T);      % mol/L
ctot = ctot*1000;    % mol/cm3
jcount = 0;          % current iteration
dC(1:3) = 1e-4;      % Delta C = small variation in value of C
dC(4:6) = 1e-10;     % Delta C = small variation in value of C
err = 1e-11;
CC = zeros(n,nj);
E = zeros(n,n+1,nj); % placeholder matrix
    
while jcount < 2
    jcount = jcount+1;   % update iteration
    CC = C;              % initialize CC
    j = 0;
    while j < nj
        j = j+1;         % move to next mesh point
        [C,E]=autoband(j,n,nj,C,dC,E);
    end
    for j = 1:nj
        for i = 1:n
            C(i,j) = CC(i,j)+C(i,j);
        end
    end
end

x = 0:L/(nj-1):L;

plot(x,C(1,:))
hold on
plot(x,C(2,:))
plot(x,C(3,:))
xlabel('z (cm)')
ylabel('x')
title('Solution of Maxwell-Stefan Equations')
legend('acetone','methanol','air','Location','Best')
