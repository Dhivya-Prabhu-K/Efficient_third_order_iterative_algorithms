clear all
clc;
n=input('value of n=');
    % Get data from each method
    [x1, y1] = r_error_fom_2(n);          % First: Red
    [x2, y2] = r_error_som_2(n);          % Second: Black
    %[x3 ,y3] = hermite_errors_asy(n); % Third: Green (returns only y)
    [x4, y4] = r_error_tom_2(n);          % Fourth: Blue

    % Define index for y3 if x not available
   

    % Plot each curve with specific color and marker
    figure;
    semilogy(x1, y1, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM');
    hold on;
    semilogy(x2, y2, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'SOM');
    %semilogy(x3,  y3, '^-', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'ASY');
    semilogy(x4, y4, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'TOM');

    % Labels and grid
    grid on;
    xlabel('zeros of H_{10000}(x)');
    ylabel('Relative error');
    % title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
    legend('Location', 'best');
    set(gca, 'FontSize', 12);
%%

figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-H');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM-H');  % Black
%plot(x3, y3, '-', 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'ASY');  % Green
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue

% Labels and grid
grid on;
xlabel('zeros of H_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
legend('Location', 'best');
set(gca, 'FontSize', 12);
%%
 s1=y1;
s2=y2;
% s3=time_method3;
s4=y4;
a=zeros(1, length(s1));
for k=1:length(s1)
    a(k) = min([s1(k), s2(k), s4(k)]);
end
m=1.2;
for i=1:length(s1)
    t1(i) = s1(i)/a(i);
    t2(i) = s2(i)/a(i);
    
    t4(i) = s4(i)/a(i);
end

tow = 0.9:0.00005:m;

n1 = zeros(1, length(tow));
n2 = zeros(1, length(tow));

n4 = zeros(1, length(tow));

for i = 1:length(tow)
    n1(i) = sum(t1 < tow(i));
    n2(i) = sum(t2 < tow(i));
   
    n4(i) = sum(t4 < tow(i));
end

% Normalize
n1 = n1 / length(s1);
n2 = n2 / length(s1);

n4 = n4 / length(s1);

% Plot
figure;
plot(tow, n1, 'r', 'LineWidth', 2); hold on;
plot(tow, n2, 'k', 'LineWidth', 2);
% plot(tow, n3, 'g', 'LineWidth', 2);
plot(tow, n4, 'b', 'LineWidth', 2);
hold off;

xlabel('\tau');
ylabel('\rho_s(\tau)');
legend({'FOM-H','SOM-H','TOM-H'}, 'Location','northwest');
title('Comparison Plot');
grid on;
