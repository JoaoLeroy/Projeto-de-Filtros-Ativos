% Joao Pedro Mendes Leroy - 220950008
clc;
clear;
close all;

% ========================================== Entradas interativas ===================================
Amax = input('Insira o valor de Amax [dB] (ex: 1): ');
Amin = input('Insira o valor de Amin [dB] (ex: 40): ');

% Verifica Fp < Fs
while true
    Fp = input('Insira a frequencia de passagem Fp [Hz] (ex: 1000): ');
    Fs = input('Insira a frequencia de rejeicao Fs [Hz] (ex: 2000): ');
    if Fp < Fs
        break;
    else
        disp('Erro: Fp deve ser menor que Fs. Tente novamente.');
    end
end

K = input('Insira o ganho K (ex: 10): ');

Wp = 2 * pi * Fp; % Conversão para rad/s da frequencia de passagem
Ws = 2 * pi * Fs; % Conversão para rad/s da frequencia de rejeicao

fprintf('\nTipos de Filtro:\n 1 - Butterworth\n 2 - Chebyshev\n'); %escolhe o tipo de topologia 
entrada1 = input('Escolha o tipo de filtro: ');
fprintf('\nTopologias:\n 1 - Sallen-Key\n 2 - MFB\n');
entrada2 = input('Escolha a topologia: ');

C1v = input('Insira o valor de C1 [F] (ex: 10e-9): ');  %escolhe o capacitor dos circuitos
fprintf('______________________________________________________________________________________________________________________\n');

% =================================== Calculo dos polos =====================================
if entrada1 == 1  % Butterworth
    e = sqrt(10^(Amax/10) - 1); % É o "ripple" q p Amax gera (quase sempre zero)
    N = ceil((log10((10^(Amin/10) - 1)/(e^2))) / (2 * log10(Ws/Wp)));  % Ordem do filtro
    W0 = Wp * (1/e)^(1/N); %frequencia natural 
    polos = [];  %for pra calcular todos os polos da função de transferência 
    for i = 1:N     
        theta = pi * (2*i - 1)/(2*N);
        p = W0 * (-sin(theta) + 1j * cos(theta));
        polos = [polos; p];
    end
else   % Chebyshev
    e = sqrt(10^(Amax/10) - 1);  % É o "ripple" q p Amax gera (quase sempre zero)
    N = ceil(acosh(sqrt((10^(Amin/10) - 1)/e^2)) / acosh(Ws/Wp));  % Ordem do filtro
    polos = [];
    for i = 1:N  %for pra calcular todos os polos da função de transferência
        alpha = sinh(1/N * asinh(1/e));
        beta  = cosh(1/N * asinh(1/e));
        theta = pi * (2*i - 1)/(2*N);
        p = -Wp * sin(theta) * alpha + 1j * Wp * cos(theta) * beta;
        polos = [polos; p];
    end
    W0 = mean(abs(polos));  % Estimativa da frequência natural média dos polos pra calcular os parâmetros dos circuitos
end

% ================================== Exibir polos =============================================
figure('Name', 'Polos do Filtro', 'NumberTitle', 'off');

% --- Subplot 1: aq tá no plano S
subplot(1,2,1)
plot(real(polos), imag(polos), 'rx'); hold on;
sgrid;    % insere linhas de zeta e wn
grid off;
xlabel('Parte Real');
ylabel('Parte Imaginária');
title('Polos no plano S');

% --- Subplot 2: Apenas os polos sem estar no plano S
subplot(1,2,2)
plot(real(polos), imag(polos), 'rx'); hold on;
grid on;
xlabel('Parte Real');
ylabel('Parte Imaginária');
title('Polos com Grade Normal');


disp('Polos do filtro:');  % Exibe também os valores numéricos na command window
disp(polos);

fprintf('______________________________________________________________________________________________________________________\n');
% ============================= Funcao de Transferencia Completa ===============================

if entrada1 == 1
    Tnum = K * W0^N; % Numerador da função de transferência (Butterworth)
else
    Tnum = (K * Wp^N)/(e*2^(N-1)); % Numerador para Chebyshev %adicio
end
Tden = poly(polos.'); % Faz o denominador usando os polos
T = tf(Tnum, real(Tden)); %junta o numerador com denominador pra poder plotar 
disp('Funcao de Transferencia Completa:');
display(T);

fprintf('______________________________________________________________________________________________________________________\n');
% ============================= Decomposição em blocos de 1ª e 2ª ordem =========================
z = [];  % zeros (não usados aqui)
p = polos;
k = Tnum;

blocos = {};  %Cria a lista dos blocos individuais
usado = false(size(p));  % marca os polos já usados pra não repetir

fprintf('\n--- Funções de Transferência em Blocos ---\n');
bloco_idx = 1;

for i = 1:length(p)
    if usado(i)
        continue;
    end
    pi1 = p(i);
    encontrado = false;

    % Tenta encontrar conjugado, ai vira um bloco de segunda ordem
    for j = i+1:length(p)
        if ~usado(j)
            pi2 = p(j);
            if abs(real(pi1) - real(pi2)) < 1e-2 && abs(imag(pi1) + imag(pi2)) < 1e-2  %margem de erro entre os polos pra poder agrupar na msm função
                den = poly([pi1, pi2]);
                bloco = tf([k], real(den));
                wn = sqrt(real(den(3)));  % Frequência natural do bloco
                zeta = real(den(2)) / (2 * wn);   % Amortecimento 
                Q = 1 / (2 * zeta);  % Fator de qualidade
                fprintf('__________________________________________\n');
                fprintf('Bloco %d (2ª ordem):\n', bloco_idx);
                display(bloco)
                fprintf('  → Wn = %.2f rad/s\n', wn);
                fprintf('  → Q  = %.3f\n', Q);
                blocos{end+1} = bloco;
                usado(i) = true;
                usado(j) = true;
                encontrado = true;
                break;
            end
        end
    end

    if ~encontrado
        % Quando é de primeira ordem
        den = poly(pi1);
        bloco = tf([k], real(den));
        fprintf('__________________________________________\n');
        fprintf('Bloco %d (1ª ordem):\n', bloco_idx);
        display(bloco)
        fprintf('  → Wn = %.2f rad/s (pólo real)\n', abs(real(pi1)));
        blocos{end+1} = bloco;
        usado(i) = true;
    end

    bloco_idx = bloco_idx + 1;
end

fprintf('______________________________________________________________________________________________________________________\n');

% =========================================== Diagrama de Bode ========================================
x = logspace(log10(Wp/100), log10(100*Wp), 1000);    % Frequências em log
H = abs(polyval(Tnum, 1j*x) ./ polyval(Tden, 1j*x)); % Resposta em frequência
Hdb = 20 * log10(H); % Converte para dB
fase = zeros(size(x));
for i = 1:length(polos)
    wp = abs(polos(i));
    fase = fase - atan(x/wp);
end
fase_deg = rad2deg(fase);

figure;
subplot(2,1,1)
semilogx(x, Hdb, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequencia [rad/s]');
ylabel('Magnitude [dB]');
title('Diagrama de Bode - Magnitude');

subplot(2,1,2)
semilogx(x, fase_deg, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Frequencia [rad/s]');
ylabel('Fase [°]');
title('Diagrama de Bode - Fase');

sgtitle('Diagrama de Bode Completo');
drawnow;

% ===================================== Cálculo dos componentes ======================================
fprintf('Calculo dos componentes (valores reais e comerciais - serie E24):\n');
C1_real = []; C2_real = []; R1_real = []; R2_real = []; R3_real = [];
C1_com = []; C2_com = []; R1_com = []; R2_com = []; R3_com = [];

if entrada2 == 1  % Sallen-Key !!!!!!!!!!!!
    for i = 1:floor(length(polos)/2)
        polo = polos(i);   %calculo dos conpoenetes sem usar os blocos, aq foi separado pra usar os números crus e fieis
        zeta = -real(polo)/abs(polo);
        Q = 1 / (2 * zeta);
        n = 1.1 * 4 * Q^2;  %delta foi fixado como 1.1
        m = (-(2 - (n/Q^2)) + sqrt((2 - (n/Q^2))^2 - 4)) / 2;
        C2v = n * C1v;
        R2v = 1 / (W0 * C1v * sqrt(m * n));
        R1v = m * R2v;
        C1_real(end+1) = C1v; 
        C2_real(end+1) = C2v;
        R1_real(end+1) = R1v;
        R2_real(end+1) = R2v;
        C1_com(end+1) = nearestE24(C1v);
        C2_com(end+1) = nearestE24(C2v);
        R1_com(end+1) = nearestE24(R1v);
        R2_com(end+1) = nearestE24(R2v);
    end
    T = table(C1_real.', C1_com.', C2_real.', C2_com.', R1_real.', R1_com.', R2_real.', R2_com.', ...
        'VariableNames', {'C1_real','C1_com','C2_real','C2_com','R1_real','R1_com','R2_real','R2_com'});
    disp(T);
else % MFB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    for i = 1:floor(length(polos)/2)
        polo = polos(i);   %calculo dos conponentes sem usar os blocos, aq foi separado pra usar os números crus e fieis
        zeta = -real(polo)/abs(polo);
        Q = 1 / (2 * zeta);
        g = 1 + K;  %alfa do caderno 
        n = 1.1 * g * 4 * Q^2;  %delta foi fixado como 1.1
        a = (2/g) - (n / (g^2 * Q^2));
        m = (-a + sqrt(a^2 - 4/(g^2)))/2; %4º formula do caderno de baixo pra cima 
        C2v = n * C1v;
        R2v = 1 / (W0 * C1v * sqrt(m * n));
        R3v = m * R2v;
        R1v = R3v / K;
        C1_real(end+1) = C1v; 
        C2_real(end+1) = C2v;
        R1_real(end+1) = R1v; 
        R2_real(end+1) = R2v; R3_real(end+1) = R3v;
        C1_com(end+1) = nearestE24(C1v);
        C2_com(end+1) = nearestE24(C2v);
        R1_com(end+1) = nearestE24(R1v);
        R2_com(end+1) = nearestE24(R2v);
        R3_com(end+1) = nearestE24(R3v);
    end
    componentes = table(C1_real.', C1_com.', C2_real.', C2_com.', R1_real.', R1_com.', ...
              R2_real.', R2_com.', R3_real.', R3_com.','VariableNames', {'C1_real','C1_com','C2_real','C2_com','R1_real','R1_com','R2_real','R2_com','R3_real','R3_com'});
    disp(componentes);
end

fprintf('______________________________________________________________________________________________________________________\n');

% ================================== Mostrar imagem da topologia escolhida =================================

% Caminhos das imagens (ajuste os nomes se necessário)
img_SK = 'sallen_key.png';  % arquivo com png do Sallen-Key
img_MFB = 'mfb.png';        % arquivo com png do MFB

% Verifica qual topologia foi escolhida e exibe a imagem
figure('Name', 'Topologia do Filtro', 'NumberTitle', 'off');
if entrada2 == 1
    imshow(img_SK);
    title('Topologia: Sallen-Key');
else
    imshow(img_MFB);
    title('Topologia: Multiple Feedback (MFB)');
end
% ================================ Funcao auxiliar: aproximacao para serie E24 ==================================
function valor = nearestE24(x)   %so criou uma lista coms os valores comerciais da série E24
    decades = 10.^(-12:12);
    e24 = [1.0 1.1 1.2 1.3 1.5 1.6 1.8 2.0 2.2 2.4 2.7 3.0 ...
           3.3 3.6 3.9 4.3 4.7 5.1 5.6 6.2 6.8 7.5 8.2 9.1];
    baseList = sort(unique(e24' * decades));  %aq ele analisa independente do 10 elevado
    [~, idx] = min(abs(baseList - x));    %acha o valor mais próximo de x
    valor = baseList(idx);  %retorna o valor comercial
end
