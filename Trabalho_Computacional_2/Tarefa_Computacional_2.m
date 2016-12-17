% Programa para analisar os esforços internos numa barra de um mecanismo,
% em função do ângulo theta no qual o mecanismo é posicionado. 

%% Preparação do MATLAB
close all
clear all
clc

%% Definição das constantes no SI
L1 = 2; % Tamanho da barra levantadora
L2 = 0.8; % Tamanho da barra suporte
m1 = 20; % Massa da barra levantadora
m2 = 11; % Massa da barra suporte
m = 1000; % Massa pontual levantada pelo sistema
b = 0.5; % Distância entre os suportes das barras
g = 9.81; % Gravidade

Pd = m*g; % Peso da massa pontual
Pbd = m1*g; % Peso da barra levantadora
Pac = m2*g; % Peso da barra suporte

%% Definição das forças aplicadas sobre a barra para ângulos de theta a 2pi
theta = linspace(0, pi/2, 6000);

alpha = asin(b .* sin(theta) ./ L2) + theta; % Define alfa em função de theta
L3 = L2 .* sin(alpha) ./ sin(theta); % Define a distância entre B e a manga

Rc = (L1 .* cos(theta) .* 1./(2.*L3)).*(2.*Pd + Pbd); % Valor da reação em C
Rbx = Rc .* sin(theta); % Valor da reação em B na direção x
Rby = (Rc .* cos(theta)) - Pd - Pbd; % Valor da reação em B na direção y
Rax = Rbx; % Valor da reação A na direção x
Ray = Pd + Pac + Pbd + Rby; % Valor da reação B na direção x
T = (Pac .* L2./2 .* cos(alpha)) + (Pd .* ((L1 .* cos(theta)) - b)) ...
    + (Pbd .* ((L1./2 .* cos(theta)) - b)) - (Rby .* b); % Valor do torque

%% Cálculo dos esforços internos
x = linspace(0, L1, 10000); % Cria um vetor de valores de x

% Inicialização das matrizes (reduz tempo do loop)
Nx = zeros(6000, 10000);
Vy = zeros(6000, 10000);
Mz = zeros(6000, 10000);

% Por meio de um for, criamos matrizes de Nx, Vy e Mz, em que cada linha
% apresenta todos os valores do esforço interno em função da posição para
% cada um dos ângulos

for angle = (1:1:6000)
    for pos = (1:1:10000)        
        % No caso em que a reação encontra-se acima da metade da barra
        if L3(angle) > L1/2
            % Seção 1 - Entre 0 e L1-L3
            if x(pos) <= L1 - L3(angle)
                Nx(angle, pos) = -Pd * sin(theta(angle));
                Vy(angle, pos) = -Pd * cos(theta(angle));
                Mz(angle, pos) = -Pd * cos(theta(angle)) * x(pos);
            % Seção 2 - Entre L1-L3 e L1/2  
            elseif x(pos) > L1 - L3(angle) && x(pos) <= L1/2
                Nx(angle, pos) = -Pd * sin(theta(angle));
                Vy(angle, pos) = Rc(angle) - (Pd * cos(theta(angle)));
                Mz(angle, pos) = (Rc(angle) * (x(pos) - ...
                   (L1 - L3(angle)))) - (Pd * cos(theta(angle)) * x(pos));
            % Seção 3 - Entre L1/2 e L1
            elseif x(pos) > L1/2 && x(pos) <= L1
                Nx(angle, pos) = - (Pd * sin(theta(angle))) - ... 
                   (Pbd * sin(theta(angle)));
                Vy(angle, pos) = Rc(angle) - (Pd * cos(theta(angle))) - ...
                   (Pbd * cos(theta(angle)));
                Mz(angle, pos) = (Rc(angle) * (x(pos) - ...
                   (L1 - L3(angle)))) - (Pd * cos(theta(angle)) * ...
                   x(pos)) - ((Pbd * cos(theta(angle))) * (x(pos) - L1/2));
            end
        % No caso em que a reação fica abaixo da metada da barra
        else
            % Seção 1 - Entre 0 e L1/2
            if x(pos) <= L1/2
                Nx(angle, pos) = -Pd * sin(theta(angle));
                Vy(angle, pos) = -Pd * cos(theta(angle));
                Mz(angle, pos) = -Pd * cos(theta(angle)) * x(pos);
            % Seção 2 - Entre L1/2 e L1-L3
            elseif x(pos) > L1/2 && x(pos) <= L1 - L3(angle)
                Nx(angle, pos) = - (Pd * sin(theta(angle))) - ... 
                   (Pbd * sin(theta(angle)));
                Vy(angle, pos) = - (Pd * cos(theta(angle))) - ...
                   (Pbd * cos(theta(angle)));
                Mz(angle, pos) = - (Pd * cos(theta(angle)) * x(pos)) - ...
                   ((Pbd * cos(theta(angle))) * (x(pos) - L1/2));
            % Seção 3 - Entre L1-L3 e L1
            elseif x(pos) > L1/2 && x(pos) <= L1
                Nx(angle, pos) = - (Pd * sin(theta(angle))) - ... 
                   (Pbd * sin(theta(angle)));
                Vy(angle, pos) = Rc(angle) - (Pd * cos(theta(angle))) - ...
                   (Pbd * cos(theta(angle)));
                Mz(angle, pos) = (Rc(angle) * (x(pos) - ...
                   (L1 - L3(angle)))) - (Pd * cos(theta(angle)) * ...
                   x(pos)) - ((Pbd * cos(theta(angle))) * (x(pos) - L1/2));
            end
        end
    end
end

%% Definição dos esforços para os ângulos indicados e plotagem
% Seleciona os valores exclusivos para o ângulo de pi/6
Nx_30 = Nx(2000, :)/1000;
Vy_30 = Vy(2000, :)/1000;
Mz_30 = Mz(2000, :)/1000;

% Seleciona os valores exclusivos para o ângulo de pi/3
Nx_60 = Nx(4000, :)/1000;
Vy_60 = Vy(4000, :)/1000;
Mz_60 = Mz(4000, :)/1000;

subplot(3, 2, 1);
plot(x, Nx_30, 'r');
title('Esforço Normal para \pi/6');
xlabel('Posição (m)');
ylabel('Nx (kN)');
grid on;

subplot(3, 2, 2);
plot(x, Nx_60, 'r');
title('Esforço Normal para \pi/3');
xlabel('Posição (m)');
ylabel('Nx (kN)');
grid on;

subplot(3, 2, 3);
plot(x, Vy_30, 'r');
title('Esforço Cortante para \pi/6');
xlabel('Posição (m)');
ylabel('Vy (kN)');
grid on;

subplot(3, 2, 4);
plot(x, Vy_60, 'r');
title('Esforço Cortante para \pi/3');
xlabel('Posição (m)');
ylabel('Vy (kN)');
grid on;

subplot(3, 2, 5);
plot(x, Mz_30, 'r');
title('Momento Fletor para \pi/6');
xlabel('Posição (m)');
ylabel('Mz (kN)');
grid on;

subplot(3, 2, 6);
plot(x, Mz_60, 'r');
title('Momento Fletor para \pi/3');
xlabel('Posição (m)');
ylabel('Mz (kN)');
grid on;