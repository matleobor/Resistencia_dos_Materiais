% Programa para:
% I) Calcular centroide, momento de inércia, gráficos de esforços internos
% e tensões em pontos específicos para uma viga de seção transversal 
% previamente definida
% II) Otimizar as dimensões da seção transversal de uma viga para suportar
% para valores máximos de tensão e deslocamento

%% Preparação do MATLAB
close all
clear all
clc

%% PARTE 1 - Seção transversal definida

% Definição das Constantes
L = 4.5; % Comprimento da viga (m)
E = 200e9; % Módulo de Elasticidade (Pa)
a = 2e-3; % Dimensões das seções (m) 
b = 10e-3;
c = 1e-3;
w = 10e-3;
x = 1e-3;
y = 3e-3;
q = 30e3; % Carga distribuída (N.m)
M = 45e3; % Momento pontual (N.m)
Rb = 41.25e3; % Reação B (N)
X = linspace(0, L, 1000); % Vetor posição na viga

% Cálculo do Centroide (em relação à base)
area = (b*a) + (2*x*y) + (w*c); % Define a área total
center_ba = b/2; % Define a posição dos centroides de cada parte
center_wc = b + (c/2);
center_xy = b + c - (y/2);
% Utiliza a somatória para cálculo de centroide
centroid = ((center_ba*b*a) + (2*center_xy*x*y) + (center_wc*w*c))/area;

% Cálculo do Momento de Inércia
I_ba = (a*b^3)/12; % Define os momentos de inércia de cada parte
I_wc = (w*c^3)/12;
I_xy = (x*y^3)/12;
d_ba = center_ba - centroid; % Define as distâncias entre os centroides das 
d_wc = center_wc - centroid; % partes e o centroide da seção total
d_xy = center_xy - centroid;
% Calcula o momento de inércia da seção total
Izz = (I_ba + (b*a*d_ba^2)) + (I_wc + (w*c*d_wc^2)) + ...
    2*(I_xy + (x*y*d_xy^2));

% Equações de Esforço Cortante e Momento Fletor
Vy = (-q*X) + (q.*(X - 1.5) + Rb).*(X > 1.5);
Mz = (-q/2*X.^2) + (q/2*(X - 1.5) + Rb).*(X - 1.5).*(X > 1.5) + M.*(X > 3);

% Cálculo das Tensões
M_max = max(abs(Mz)); % Definição do máximo momento absoluto
d_A = (b + c - y) - centroid; % Distâncias dos pontos à linha neutra
d_B = b - centroid;
sigma_A = (M_max*d_A*1e-6)/Izz; % Cálculo das tensões em MPa
sigma_B = (M_max*d_B*1e-6)/Izz;

% Resultados e Plotagem dos Gráficos
fprintf('---------------------------------------------------------\n');
fprintf('\t\t\t\t\t\tPARTE 1\n');
fprintf('---------------------------------------------------------\n');
fprintf('Altura do centroide em relação à base: %.2e m\n', centroid);
fprintf('Momento de inércia: %e m^4\n', Izz);
fprintf('Momento fletor máximo: %.2f kN.m\n', M_max*1e-3);
fprintf('Tensão de flexão em A: %.2f MPa\n', sigma_A);
fprintf('Tensão de flexão em B: %.2f MPa\n', sigma_B);
fprintf('\n');

% Os termos comentados são ajustes para melhor aspecto dos gráficos quando
% exportados como imagem
figure(1);
subplot(1, 2, 1);
plot(X, Vy*1e-3, '-r');
title('Esforço Cortante');
xlabel('Posição (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(1, 2, 2);
plot(X, Mz*1e-3, '-r');
title('Momento Fletor');
xlabel('Posição (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

% set(gcf, 'papersize', [17.6, 8]);
% set(gcf, 'paperposition', [0,0,17.6, 8]);
% print('graficos_4', '-dbmp');

%% PARTE 2 - Seção transversal para otimização

% Definição das Constantes
clear a; % Limpa os valores de a previamente definidos
rho = 7900; % Densidade Aço 1020 (kg/m^3)
Sy = 352; % Limite de Escoamento Aço 1020 (MPa)
g = 9.81; % Gravidade (m/s^2)
epsilon = 1e-5; % Precisão dos cálculos
a = epsilon; % Valores iniciais para entrada no loop
sigma_max = Sy;
v_max = L;

% Definição da Dimensão a
while (sigma_max > Sy/2) || (v_max > 0.002*L) % Enquanto os critérios não 
                                              % são satisfeitos
    b = 32*a; % Definição das dimensões em função de a
    c = 2*a;
    w = 3*a;
    x = 0.1*a;
    y = 15*a;
    
    % Cálculo do Centroide (em relação à base)
    area = (b*a) + (2*x*y) + (w*c); % Define a área total
    center_ba = b/2; % Define a posição dos centroides de cada parte
    center_wc = b + (c/2);
    center_xy = b + c - (y/2);
    % Utiliza a somatória para cálculo de centroide
    centroid = ((center_ba*b*a) + (2*center_xy*x*y) + (center_wc*w*c))/ ...
        area;
    
    % Cálculo do Momento de Inércia
    I_ba = (a*b^3)/12; % Define os momentos de inércia de cada parte
    I_wc = (w*c^3)/12;
    I_xy = (x*y^3)/12;
    d_ba = center_ba - centroid; % Define as distâncias entre os centroides
    d_wc = center_wc - centroid; % das partes e o centroide da seção total
    d_xy = center_xy - centroid;
    % Calcula o momento de inércia da seção total
    Izz = (I_ba + (b*a*d_ba^2)) + (I_wc + (w*c*d_wc^2)) + ...
        2*(I_xy + (x*y*d_xy^2));
    
    % Cálculo do Peso Distribuído, Reação B e Constantes (utilização do
    % método da equação diferencial)
    p = area*rho*g;
    Rb = (123.75e3 + 2.25*p*L)/3;
    c1 = (1/L)*((q/2*L^2) - (q/2*(L - 1.5)^2) - Rb*(L - 1.5) - M + ...
        (p/2*L^2));
    c3 = (1/(L - 1.5))*(- c1/6*((L^3) - (1.5)^3) + p/24*((L^4) - ...
        (1.5)^4) - M/2*((L - 3)^2) - Rb/6*((L - 1.5)^3) - ...
        q/24*((L - 1.5)^4) + q/24*((L^4) - (1.5)^4));
    c4 = 1/24*(1.5^4)*(q + p) - 1.5*c3;
    
    % Equações de Esforço Cortante, Momento Fletor e Deslocamentos
    % (utilização do método da equação diferencial)
    Vy = (-p*X) - (q*X) + (q.*(X - 1.5) + Rb).*(X > 1.5) + c1;
    Mz = (-p/2*X.^2) - (q/2*X.^2) + (q/2*(X - 1.5) + Rb).*(X - 1.5).* ...
        (X > 1.5) + M.*(X > 3) + c1*X;
    theta = (1/(E*Izz))*((-p/6*X.^3) - (q/6*X.^3) + ...
        (q/6*(X - 1.5) + Rb/2).*((X - 1.5).^2).*(X > 1.5) + M.*(X - 3).*...
        (X > 3) + c1/2*X.^2 + c3);
    v = (1/(E*Izz))*((-p/24*X.^4) - (q/24*X.^4) + ...
        (q/24*(X - 1.5) + Rb/6).*((X - 1.5).^3).*(X > 1.5) + ...
        M/2.*((X - 3).^2).*(X > 3) + c1/6*X.^3 + c3*X + c4);
    
    % Cálculo da Tensão e Deslocamento Máximos
    M_max = max(abs(Mz)); % Definição do máximo momento fletor absoluto
    v_max = max(abs(v)); % Definição do máximo deslocamento
    d_max = centroid; % Maior distância à linha neutra (maior tensão)
    sigma_max = (M_max*d_max*1e-6)/Izz; % Cálculo da tensão em MPa
        
    a = a + epsilon; % Aumenta a um épsilon previamente definido (precisão)
end % Ao final do loop, é obtido o valor de a que cumpre os requisitos 
    % previamente definidos

% Resultados e Plotagem dos Gráficos
fprintf('---------------------------------------------------------\n');
fprintf('\t\t\t\t\t\tPARTE 2\n');
fprintf('---------------------------------------------------------\n');
fprintf('Dimensões otimizadas (mm):\n');
fprintf('a = %.2f\tb = %.2f\t c = %.2f\t\n', (a - epsilon)*1e3, ...
    b*1e3, c*1e3);
fprintf('w = %.2f\tx = %.2f\t y = %.2f\t\n', w*1e3, x*1e3, y*1e3);
fprintf('---------------------------------------------------------\n');
fprintf('Massa final da viga: %.2f kg\n', L*area*rho);
fprintf('Altura do centroide em relação à base: %.2f m\n', centroid);
fprintf('Momento de inércia: %e m^4\n', Izz);
fprintf('Tensão máxima de flexão: %.2f MPa\n', sigma_max);
fprintf('Deslocamento vetical máximo: %.2e m\n', v_max);

% Os termos comentados são ajustes para melhor aspecto dos gráficos quando
% exportados como imagem
figure(2);
subplot(1, 2, 1);
plot(X, theta*1e3, '-r');
title('Deslocamento Angular');
xlabel('Posição (m)');
ylabel('rad (x10^{3})');
% set(gca, 'FontSize', 25);
grid on;

subplot(1, 2, 2);
plot(X, v*1e3, '-r');
title('Deslocamento Vertical');
xlabel('Posição (m)');
ylabel('mm');
% set(gca, 'FontSize', 25);
grid on;

% set(gcf, 'papersize', [17.6, 8]);
% set(gcf, 'paperposition', [0,0,17.6, 8]);
% print('graficos_5', '-dbmp');

% Caso se queira plotar os gráficos de esforços internos e tensão para um
% ponto, basta escolher um valor de d, a distância do ponto à linha neutra,
% em metros
d = 19e-2;  % Distância do ponto à linha neutra
sigma = (Mz*d*1e-6)/Izz; % Cria o vetor de tensões à distância d da linha 
                         % neutra

figure(3);
subplot(1, 3, 1);
plot(X, Vy*1e-3, '-r');
title('Esforço Cortante');
xlabel('Posição (m)');
ylabel('kN');
grid on;

subplot(1, 3, 2);
plot(X, Mz*1e-3, '-r');
title('Momento Fletor');
xlabel('Posição (m)');
ylabel('kN.m');
grid on;

subplot(1, 3, 3);
plot(X, sigma, '-r');
title('Tensão no ponto d');
xlabel('Posição (m)');
ylabel('MPa');
grid on;