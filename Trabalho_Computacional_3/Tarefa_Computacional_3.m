% Programa para analisar a equivalência entre um carregamento distribuído e
% um momento pontual quando aplicados em uma viga bi-apoiada

%% Preparação do MATLAB
close all
clear all
clc

%% Definição das Constantes e Incialização das Matrizes
L = 2; % Comprimento da Viga
Vy_eq = zeros(1000); % Esforço cortante do modelo pontual
Mz_eq = zeros(1000); % Momento fletor do modelo pontual
Vy_real = zeros(1000); % Esforço cortante da carga distribuída
Mz_real = zeros(1000); % Momento fletor da carga distribuída
a = linspace(0.001, L/2, 1000); % Valores possíveis para a
x = linspace(0, L, 1000); % Valores de x ao longo da viga

%% Definição das Equações de Esforço Cortante e Momento Fletor
% Valores para o modelo de momento pontual (2 seções)
for i = 1:1000 % Para cada valor de A, realiza os cálculos de esforço
    M = (10*a(i));
    for j = 1:1000;
        if x(j) <= L/2
            Vy_eq(i, j) = M/L;
            Mz_eq(i, j) = (M/L)*x(j);
        else
            Vy_eq(i, j) = M/L;
            Mz_eq(i, j) = M*((x(j)/L) - 1);
        end
    end
end

% Valores para o modelo de carga distribuída (4 seções)
for i = 1:1000 % Para cada valor de a, realiza os cálculos de esforço
    q0 = 10/a(i);
    for j = 1:1000
        if x(j) <= (L/2 - a(i))
            Vy_real(i, j) = (q0*a(i)^2)/L;
            Mz_real(i, j) = (x(j)*q0*a(i)^2)/L;
        elseif x(j) > (L/2 - a(i)) && x(j) <= L/2
            Vy_real(i, j) = q0*((a(i)^2)/L + L/2 - a(i) - x(j));
            Mz_real(i, j) = (x(j)*q0*a(i)^2)/L - ...
                (q0/2)*(x(j) - L/2 + a(i))^2;
        elseif x(j) > L/2 && x(j) <= (L/2 + a(i))
            Vy_real(i, j) = q0*((a(i)^2)/L - L/2 - a(i) + x(j));
            Mz_real(i, j) = (x(j)*q0*a(i)^2)/L - (q0*a(i))*(x(j) - L/2 ...
                + a(i)/2) + (q0/2)*(x(j) - L/2)^2;
        else
            Vy_real(i, j) = (q0*a(i)^2)/L;
            Mz_real(i, j) = (q0*a(i)^2)*((x(j)/L) - 1);
        end
    end
end

%% Plotagem dos Gráficos
% As linhas comentadas são ajsutes no tamanho dos gráficos, utilizados para
% salvar as imagens no tamanho certo
subplot(4, 2, 1);
plot(x, Vy_real(500, :), '-r');
hold on;
plot(x, Vy_eq(500, :), '-b');
title({'Esforço Cortante', 'a = 0,5m'});
xlabel('Posição (m)');
ylabel('kN');
legenda = legend('Real', 'Equivalente');
% set(legenda, 'FontSize', 8);
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 2);
plot(x, Mz_real(500, :), '-r');
hold on;
plot(x, Mz_eq(500, :), '-b');
title({'Momento Fletor', 'a = 0,5m'});
xlabel('Posição (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 3);
plot(x, Vy_real(100, :), '-r');
hold on;
plot(x, Vy_eq(100, :), '-b');
title('a = 0,1m');
xlabel('Posição (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 4);
plot(x, Mz_real(100, :), '-r');
hold on;
plot(x, Mz_eq(100, :), '-b');
title('a = 0,1m');
xlabel('Posição (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 5);
plot(x, Vy_real(50, :), '-r');
hold on;
plot(x, Vy_eq(50, :), '-b');
title('a = 0,05m');
xlabel('Posição (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 6);
plot(x, Mz_real(50, :), '-r');
hold on;
plot(x, Mz_eq(50, :), '-b');
title('a = 0,05m');
xlabel('Posição (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 7);
plot(x, Vy_real(1, :), '-r');
hold on;
plot(x, Vy_eq(1, :), '-b');
title('a = 0,001m');
xlabel('Posição (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 8);
plot(x, Mz_real(1, :), '-r');
hold on;
plot(x, Mz_eq(1, :), '-b');
title('a = 0,001m');
xlabel('Posição (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

% set(gcf, 'papersize', [18.46,25.16]);
% set(gcf, 'paperposition', [0,0,18.46,25.16]);
% print('graficos', '-dbmp');