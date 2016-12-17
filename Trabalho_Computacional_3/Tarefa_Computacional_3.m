% Programa para analisar a equival�ncia entre um carregamento distribu�do e
% um momento pontual quando aplicados em uma viga bi-apoiada

%% Prepara��o do MATLAB
close all
clear all
clc

%% Defini��o das Constantes e Incializa��o das Matrizes
L = 2; % Comprimento da Viga
Vy_eq = zeros(1000); % Esfor�o cortante do modelo pontual
Mz_eq = zeros(1000); % Momento fletor do modelo pontual
Vy_real = zeros(1000); % Esfor�o cortante da carga distribu�da
Mz_real = zeros(1000); % Momento fletor da carga distribu�da
a = linspace(0.001, L/2, 1000); % Valores poss�veis para a
x = linspace(0, L, 1000); % Valores de x ao longo da viga

%% Defini��o das Equa��es de Esfor�o Cortante e Momento Fletor
% Valores para o modelo de momento pontual (2 se��es)
for i = 1:1000 % Para cada valor de A, realiza os c�lculos de esfor�o
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

% Valores para o modelo de carga distribu�da (4 se��es)
for i = 1:1000 % Para cada valor de a, realiza os c�lculos de esfor�o
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

%% Plotagem dos Gr�ficos
% As linhas comentadas s�o ajsutes no tamanho dos gr�ficos, utilizados para
% salvar as imagens no tamanho certo
subplot(4, 2, 1);
plot(x, Vy_real(500, :), '-r');
hold on;
plot(x, Vy_eq(500, :), '-b');
title({'Esfor�o Cortante', 'a = 0,5m'});
xlabel('Posi��o (m)');
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
xlabel('Posi��o (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 3);
plot(x, Vy_real(100, :), '-r');
hold on;
plot(x, Vy_eq(100, :), '-b');
title('a = 0,1m');
xlabel('Posi��o (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 4);
plot(x, Mz_real(100, :), '-r');
hold on;
plot(x, Mz_eq(100, :), '-b');
title('a = 0,1m');
xlabel('Posi��o (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 5);
plot(x, Vy_real(50, :), '-r');
hold on;
plot(x, Vy_eq(50, :), '-b');
title('a = 0,05m');
xlabel('Posi��o (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 6);
plot(x, Mz_real(50, :), '-r');
hold on;
plot(x, Mz_eq(50, :), '-b');
title('a = 0,05m');
xlabel('Posi��o (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 7);
plot(x, Vy_real(1, :), '-r');
hold on;
plot(x, Vy_eq(1, :), '-b');
title('a = 0,001m');
xlabel('Posi��o (m)');
ylabel('kN');
% set(gca, 'FontSize', 25);
grid on;

subplot(4, 2, 8);
plot(x, Mz_real(1, :), '-r');
hold on;
plot(x, Mz_eq(1, :), '-b');
title('a = 0,001m');
xlabel('Posi��o (m)');
ylabel('kN.m');
% set(gca, 'FontSize', 25);
grid on;

% set(gcf, 'papersize', [18.46,25.16]);
% set(gcf, 'paperposition', [0,0,18.46,25.16]);
% print('graficos', '-dbmp');