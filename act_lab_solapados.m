clear; clc;

% Parámetros generales
f0 = 1;             % Frecuencia de banda base (Hz)
alphas = [0, 0.25, 0.75, 1];   

% Tiempo para respuesta al impulso 
t = linspace(0, 5/f0, 1000);     

% Gráfica de la respuesta al impulso en un solo gráfico
figure;
hold on;
colors = ['b', 'r', 'g', 'k']; 

for i = 1:length(alphas)
    alpha = alphas(i);
    fDelta = alpha * f0;

    sinc_term = sinc(2*f0*t); 
    denom = 1 - (4*fDelta*t).^2;
    denom(denom==0) = eps; 

    % Respuesta al impulso
    he_t = 2*f0 .* sinc_term .* cos(2*pi*fDelta*t) ./ denom;

    plot(t, he_t, 'LineWidth', 1.5, 'Color', colors(i));
end

title('Respuesta al impulso para diferentes \alpha');
xlabel('t (s)');
ylabel('h_e(t)');
legend('\alpha = 0', '\alpha = 0.25', '\alpha = 0.75', '\alpha = 1');
grid on;
hold off;

% Frecuencia para respuesta en frecuencia
f = linspace(-2*f0*(1+1), 2*f0*(1+1), 1000); % desde -2B hasta 2B

% Gráfica de la respuesta en frecuencia en un solo gráfico
figure;
hold on;

for i = 1:length(alphas)
    alpha = alphas(i);
    f1 = f0 * (1 - alpha);   % Límite inferior de transición
    f2 = f0 * (1 + alpha);   % Límite superior de transición

    He_f = zeros(size(f));
    abs_f = abs(f);

    % Construcción de la respuesta en frecuencia
    for k = 1:length(f)
        if abs_f(k) <= f1
            He_f(k) = 1;
        elseif abs_f(k) > f1 && abs_f(k) <= f2
            He_f(k) = 0.5 * (1 + cos(pi * (abs_f(k) - f1) / (2 * alpha * f0)));
        else
            He_f(k) = 0;
        end
    end

    plot(f, He_f, 'LineWidth', 1.5, 'Color', colors(i));
end

title('Respuesta en frecuencia para diferentes \alpha');
xlabel('f (Hz)');
ylabel('H_e(f)');
legend('\alpha = 0', '\alpha = 0.25', '\alpha = 0.75', '\alpha = 1');
grid on;
hold off;
