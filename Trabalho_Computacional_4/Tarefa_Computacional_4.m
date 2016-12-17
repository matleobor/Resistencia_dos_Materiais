% Programa para:
% I) Calcular centroide, momento de in�rcia, gr�ficos de esfor�os internos
% e tens�es em pontos espec�ficos para uma viga de se��o transversal 
% previamente definida
% II) Otimizar as dimens�es da se��o transversal de uma viga para suportar
% para valores m�ximos de tens�o e deslocamento

%% Prepara��o do MATLAB
close all
clear all
clc

%% PARTE 1 - Se��o transversal definida

% Defini��o das Constantes
L = 4.5; % Comprimento da viga (m)
E = 200e9; % M�dulo de Elasticidade (Pa)
a = 2e-3; % Dimens�es das se��es (m) 
b = 10e-3;
c = 1e-3;
w = 10e-3;
x = 1e-3;
y = 3e-3;
q = 30e3; % Carga distribu�da (N.m)
M = 45e3; % Momento pontual (N.m)
Rb = 41.25e3; % Rea��o B (N)
X = linspace(0, L, 1000); % Vetor posi��o na viga

% C�lculo do Centroide (em rela��o � base)
area = (b*a) + (2*x*y) + (w*c); % Define a �rea total
center_ba = b/2; % Define a posi��o dos centroides de cada parte
center_wc = b + (c/2);
center_xy = b + c - (y/2);
% Utiliza a somat�ria para c�lculo de centroide
centroid = ((center_ba*b*a) + (2*center_xy*x*y) + (center_wc*w*c))/area;

% C�lculo do Momento de In�rcia
I_ba = (a*b^3)/12; % Define os momentos de in�rcia de cada parte
I_wc = (w*c^3)/12;
I_xy = (x*y^3)/12;
d_ba = center_ba - centroid; % Define as dist�ncias entre os centroides das 
d_wc = center_wc - centroid; % partes e o centroide da se��o total
d_xy = center_xy - centroid;
% Calcula o momento de in�rcia da se��o total
Izz = (I_ba + (b*a*d_ba^2)) + (I_wc + (w*c*d_wc^2)) + ...
    2*(I_xy + (x*y*d_xy^2));

% Equa��es de Esfor�o Cortante e Momento Fletor
Vy = (-q*X) + (q.*(X - 1.5) + Rb).*(X > 1.5);
Mz = (-q/2*X.^2) + (q/2*(X - 1.5) + Rb).*(X - 1.5).*(X > 1.5) + M.*(X > 3);

% C�lculo das Tens�es
M_max = max(abs(Mz)); % Defini��o do m�ximo momento absoluto
d_A = (b + c - y) - centroid; % Dist�ncias dos pontos � linha neutra
d_B = b - centroid;
sigma_A = (M_max*d_A*1e-6)/Izz; % C�lculo das tens�es em MPa
sigma_B = (M_max*d_B*1e-6)/Izz;

% Resultados e Plotagem dos Gr�ficos
fprintf('---------------------------------------------------------\n');
fprintf('\t\t\t\t\t\tPARTE 1\n');
fprintf('---------------------------------------------------------\n');
fprintf('Altura do centroide em rela��o � base: %.2e m\n', centroid);
fprintf('Momento de in�rcia: %e m^4\n', Izz);
fprintf('Momento fletor m�ximo: %.2f kN.m\n', M_max*1e-3);
fprintf('Tens�o de flex�o em A: %.2f MPa\n', sigma_A);
fprintf('Tens�o de flex�o em B: %.2f MPa\n', sigma_B);
fprintf('\n');

% Os termos comentados s�o ajustes para melhor aspecto dos gr�ficos quando
% exportados como imagem
figure(1);
subplot(1, 2, 1);
plot(X, Vy*1e-3, '-r');
title('Esfor�o Cortante');
xlabel('Posi��o (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(1, 2, 2);
plot(X, Mz*1e-3, '-r');
title('Momento Fletor');
xlabel('Posi��o (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

% set(gcf, 'papersize', [17.6, 8]);
% set(gcf, 'paperposition', [0,0,17.6, 8]);
% print('graficos_4', '-dbmp');

%% PARTE 2 - Se��o transversal para otimiza��o

% Defini��o das Constantes
clear a; % Limpa os valores de a previamente definidos
rho = 7900; % Densidade A�o 1020 (kg/m^3)
Sy = 352; % Limite de Escoamento A�o 1020 (MPa)
g = 9.81; % Gravidade (m/s^2)
epsilon = 1e-5; % Precis�o dos c�lculos
a = epsilon; % Valores iniciais para entrada no loop
sigma_max = Sy;
v_max = L;

% Defini��o da Dimens�o a
while (sigma_max > Sy/2) || (v_max > 0.002*L) % Enquanto os crit�rios n�o 
                                              % s�o satisfeitos
    b = 32*a; % Defini��o das dimens�es em fun��o de a
    c = 2*a;
    w = 3*a;
    x = 0.1*a;
    y = 15*a;
    
    % C�lculo do Centroide (em rela��o � base)
    area = (b*a) + (2*x*y) + (w*c); % Define a �rea total
    center_ba = b/2; % Define a posi��o dos centroides de cada parte
    center_wc = b + (c/2);
    center_xy = b + c - (y/2);
    % Utiliza a somat�ria para c�lculo de centroide
    centroid = ((center_ba*b*a) + (2*center_xy*x*y) + (center_wc*w*c))/ ...
        area;
    
    % C�lculo do Momento de In�rcia
    I_ba = (a*b^3)/12; % Define os momentos de in�rcia de cada parte
    I_wc = (w*c^3)/12;
    I_xy = (x*y^3)/12;
    d_ba = center_ba - centroid; % Define as dist�ncias entre os centroides
    d_wc = center_wc - centroid; % das partes e o centroide da se��o total
    d_xy = center_xy - centroid;
    % Calcula o momento de in�rcia da se��o total
    Izz = (I_ba + (b*a*d_ba^2)) + (I_wc + (w*c*d_wc^2)) + ...
        2*(I_xy + (x*y*d_xy^2));
    
    % C�lculo do Peso Distribu�do, Rea��o B e Constantes (utiliza��o do
    % m�todo da equa��o diferencial)
    p = area*rho*g;
    Rb = (123.75e3 + 2.25*p*L)/3;
    c1 = (1/L)*((q/2*L^2) - (q/2*(L - 1.5)^2) - Rb*(L - 1.5) - M + ...
        (p/2*L^2));
    c3 = (1/(L - 1.5))*(- c1/6*((L^3) - (1.5)^3) + p/24*((L^4) - ...
        (1.5)^4) - M/2*((L - 3)^2) - Rb/6*((L - 1.5)^3) - ...
        q/24*((L - 1.5)^4) + q/24*((L^4) - (1.5)^4));
    c4 = 1/24*(1.5^4)*(q + p) - 1.5*c3;
    
    % Equa��es de Esfor�o Cortante, Momento Fletor e Deslocamentos
    % (utiliza��o do m�todo da equa��o diferencial)
    Vy = (-p*X) - (q*X) + (q.*(X - 1.5) + Rb).*(X > 1.5) + c1;
    Mz = (-p/2*X.^2) - (q/2*X.^2) + (q/2*(X - 1.5) + Rb).*(X - 1.5).* ...
        (X > 1.5) + M.*(X > 3) + c1*X;
    theta = (1/(E*Izz))*((-p/6*X.^3) - (q/6*X.^3) + ...
        (q/6*(X - 1.5) + Rb/2).*((X - 1.5).^2).*(X > 1.5) + M.*(X - 3).*...
        (X > 3) + c1/2*X.^2 + c3);
    v = (1/(E*Izz))*((-p/24*X.^4) - (q/24*X.^4) + ...
        (q/24*(X - 1.5) + Rb/6).*((X - 1.5).^3).*(X > 1.5) + ...
        M/2.*((X - 3).^2).*(X > 3) + c1/6*X.^3 + c3*X + c4);
    
    % C�lculo da Tens�o e Deslocamento M�ximos
    M_max = max(abs(Mz)); % Defini��o do m�ximo momento fletor absoluto
    v_max = max(abs(v)); % Defini��o do m�ximo deslocamento
    d_max = centroid; % Maior dist�ncia � linha neutra (maior tens�o)
    sigma_max = (M_max*d_max*1e-6)/Izz; % C�lculo da tens�o em MPa
        
    a = a + epsilon; % Aumenta a um �psilon previamente definido (precis�o)
end % Ao final do loop, � obtido o valor de a que cumpre os requisitos 
    % previamente definidos

% Resultados e Plotagem dos Gr�ficos
fprintf('---------------------------------------------------------\n');
fprintf('\t\t\t\t\t\tPARTE 2\n');
fprintf('---------------------------------------------------------\n');
fprintf('Dimens�es otimizadas (mm):\n');
fprintf('a = %.2f\tb = %.2f\t c = %.2f\t\n', (a - epsilon)*1e3, ...
    b*1e3, c*1e3);
fprintf('w = %.2f\tx = %.2f\t y = %.2f\t\n', w*1e3, x*1e3, y*1e3);
fprintf('---------------------------------------------------------\n');
fprintf('Massa final da viga: %.2f kg\n', L*area*rho);
fprintf('Altura do centroide em rela��o � base: %.2f m\n', centroid);
fprintf('Momento de in�rcia: %e m^4\n', Izz);
fprintf('Tens�o m�xima de flex�o: %.2f MPa\n', sigma_max);
fprintf('Deslocamento vetical m�ximo: %.2e m\n', v_max);

% Os termos comentados s�o ajustes para melhor aspecto dos gr�ficos quando
% exportados como imagem
figure(2);
subplot(1, 2, 1);
plot(X, theta*1e3, '-r');
title('Deslocamento Angular');
xlabel('Posi��o (m)');
ylabel('rad (x10^{3})');
% set(gca, 'FontSize', 25);
grid on;

subplot(1, 2, 2);
plot(X, v*1e3, '-r');
title('Deslocamento Vertical');
xlabel('Posi��o (m)');
ylabel('mm');
% set(gca, 'FontSize', 25);
grid on;

% set(gcf, 'papersize', [17.6, 8]);
% set(gcf, 'paperposition', [0,0,17.6, 8]);
% print('graficos_5', '-dbmp');

% Caso se queira plotar os gr�ficos de esfor�os internos e tens�o para um
% ponto, basta escolher um valor de d, a dist�ncia do ponto � linha neutra,
% em metros
d = 19e-2;  % Dist�ncia do ponto � linha neutra
sigma = (Mz*d*1e-6)/Izz; % Cria o vetor de tens�es � dist�ncia d da linha 
                         % neutra

figure(3);
subplot(1, 3, 1);
plot(X, Vy*1e-3, '-r');
title('Esfor�o Cortante');
xlabel('Posi��o (m)');
ylabel('kN');
grid on;

subplot(1, 3, 2);
plot(X, Mz*1e-3, '-r');
title('Momento Fletor');
xlabel('Posi��o (m)');
ylabel('kN.m');
grid on;

subplot(1, 3, 3);
plot(X, sigma, '-r');
title('Tens�o no ponto d');
xlabel('Posi��o (m)');
ylabel('MPa');
grid on;