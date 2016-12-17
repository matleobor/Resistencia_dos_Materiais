% Programa para analisar as rea��es e o torque necess�rios para 
% manter equilibrado um mecanismo de um grau de liberdade, em 
% fun��o de um �ngulo theta. 

%% Prepara��o do MATLAB
close all
clear all
clc

%% Defini��o das constantes no SI
L1 = 2; % Tamanho da barra levantadora
L2 = 0.8; % Tamanho da barra suporte
m1 = 20; % Massa da barra levantadora
m2 = 11; % Massa da barra suporte
m = 1000; % Massa pontual levantada pelo sistema
b = 0.5; % Dist�ncia entre os suportes das barras
g = 9.81; % Gravidade

Pd = m*g; % Peso da massa pontual
Pbd = m1*g; % Peso da barra levantadora
Pac = m2*g; % Peso da barra suporte

%% Defini��o das rela��es trigonom�tricas e equa��es de equil�brio
syms theta;

alpha = asin(b * sin(theta) / L2) + theta; % Define alfa em fun��o de theta
L3 = L2 * sin(alpha) / sin(theta); % Define a dist�ncia entre B e a manga

Rc = (L1 * cos(theta) * 1/(2*L3))*(2*Pd + Pbd); % Valor da rea��o em C
Rbx = Rc * sin(theta); % Valor da rea��o em B na dire��o x
Rby = (Rc * cos(theta)) - Pd - Pbd; % Valor da rea��o em B na dire��o y
Rax = Rbx; % Valor da rea��o A na dire��o x
Ray = Pd + Pac + Pbd + Rby; % Valor da rea��o B na dire��o x
T = (Pac * L2/2 * cos(alpha)) + (Pd * ((L1 * cos(theta)) - b)) ...
    + (Pbd * ((L1/2 * cos(theta)) - b)) - (Rby * b); % Valor do torque T

%% Defini��o dos valores para 0 <= theta <= 2pi
theta = linspace(0, 2*pi, 10000); % Cria um vetor entre 0 e 2pi

% Cria, para todas as vari�veis, vetores de valores obtidos atrav�s da 
% substitui��o de theta nas equa��es de equil�brio (as divis�es por 1000 
% s�o realizadas para converter as unidades em kN)
Rc_ev = eval(Rc)/1000;
Rbx_ev = eval(Rbx)/1000;
Rby_ev = eval(Rby)/1000;
Rax_ev = eval(Rax)/1000;
Ray_ev = eval(Ray)/1000;
T_ev = eval(T)/1000;

%% Plotagem e formata��o dos gr�ficos
subplot(3, 2, 2);
plot(theta, Rc_ev, 'r');
xlim([0, 2*pi]);
title('Rea��o em C vs. Theta');
xlabel('Theta (rad)');
ylabel('Rea��o C (kN)');
grid on;

subplot(3, 2, 3);
plot(theta, Rbx_ev, 'r');
xlim([0, 2*pi]);
title('Rea��o Horizontal em B vs. Theta');
xlabel('Theta (rad)');
ylabel('Rea��o Horizontal B (kN)');
grid on;

subplot(3, 2, 4);
plot(theta, Rby_ev, 'r');
xlim([0, 2*pi]);
title('Rea��o Vertical em B vs. Theta');
xlabel('Theta (rad)');
ylabel('Rea��o Vertical B (kN)');
grid on;

subplot(3, 2, 5);
plot(theta, Rax_ev, 'r');
xlim([0, 2*pi]);
title('Rea��o Horizontal em A vs. Theta');
xlabel('Theta (rad)');
ylabel('Rea��o Horizontal A (kN)');
grid on;

subplot(3, 2, 6);
plot(theta, Ray_ev, 'r');
xlim([0, 2*pi]);
title('Rea��o Vertical em A vs. Theta');
xlabel('Theta (rad)');
ylabel('Rea��o Vertical A (kN)');
grid on;

subplot(3, 2, 1);
plot(theta, T_ev, 'r');
xlim([0, 2*pi]);
title('Torque vs. Theta');
xlabel('Theta (rad)');
ylabel('Torque (kN.m)');
grid on;

%% Torque absoluto m�ximo
fprintf('Torque m�ximo = %f kN.m\n', max(abs(T_ev)));