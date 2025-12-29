clear all
clc;
n=input('value of n=');

% Storing values
[x1, y1] = compare_leg_fourth_1(n);         
[x2, y2] = compare_leg_second_1(n);         
[x3, y3] = compare_leg_third_1(n);           


figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM'); hold on;  % Red
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM');  % Black
plot(x3, y3, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue


grid on;
xlabel('zeros of P_{10000}(x)');
ylabel('Relative error');
legend('Location', 'best');
set(gca, 'FontSize', 12);




%%
s1=y1;
s2=y2;
% s3=time_method3;
s3=y3;
a=zeros(1, length(s1));
for k=1:length(s1)
    a(k) = min([s1(k), s2(k), s3(k)]);
end
m=1.1;
for i=1:length(s1)
    t1(i) = s1(i)/a(i);
    t2(i) = s2(i)/a(i);
    
    t3(i) = s3(i)/a(i);
end

tow = 0.98:0.00005:m;

n1 = zeros(1, length(tow));
n2 = zeros(1, length(tow));

n3 = zeros(1, length(tow));

for i = 1:length(tow)
    n1(i) = sum(t1 < tow(i));
    n2(i) = sum(t2 < tow(i));
   
    n3(i) = sum(t3 < tow(i));
end

% Normalize
n1 = n1 / length(s1);
n2 = n2 / length(s1);

n3 = n3 / length(s1);

% Plot
figure;
plot(tow, n1, 'r', 'LineWidth', 2); hold on;
plot(tow, n2, 'k', 'LineWidth', 2);
% plot(tow, n3, 'g', 'LineWidth', 2);
plot(tow, n3, 'b', 'LineWidth', 2);
hold off;

xlabel('\tau');
ylabel('\rho_s(\tau)');
legend({'FOM-L','SOM-L','TOM-L'}, 'Location','northwest');
title('Comparison Plot');
grid on;
