clear all
clc;
mu=input('value of mu=');
alpha=input('value of alpha=');

    % Get data from each method
    [x1, y1] = error_cyl_fom(mu, alpha); 
    [x2, y2] = error_cyl_som(mu, alpha);          % Second: Black
    [x3 ,y3] = error_cyl_mnm(mu, alpha); % Third: Green (returns only y)
    [x4, y4] = error_cyl_tom(mu, alpha);          % Fourth: Blue

    % Define index for y3 if x not available
   

    % Plot each curve with specific color and marker
    % figure;
    % % semilogy(x1, y1, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM');
    % 
    % semilogy(x2, y2, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'SOM');
    %  hold on;
    % semilogy(x3,  y3, '^-', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'MNM');
    % semilogy(x4, y4, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'TOM');
    % 
    % % Labels and grid
    % grid on;
    % xlabel('zeros of C_{10000}(x)');
    % ylabel('Relative error');
    % % title(sprintf('Relative Error Comparison of cylinder (mu = %d)', mu));
    % legend('Location', 'best');
    % set(gca, 'FontSize', 12);
%%

figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM-C');  % Black

plot(x3, y3, '-', 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'MNM-C');  % Green
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue

% Labels and grid
grid on;
xlabel('zeros of C_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of cylinder (mu = %d)', mu));
legend('Location', 'best');
set(gca, 'FontSize', 12);


% s1=y1;
% s2=y2;
% s3=y3;
% s4=y4;
% a=zeros(1, length(s1));
% for k=1:length(s1)
%     a(k) = min([s1(k), s2(k), s3(k), s4(k)]);
% end
% m=1.5;
% for i=1:length(s1)
%     t1(i) = s1(i)/a(i);
%     t2(i) = s2(i)/a(i);
%     t3(i) = s3(i)/a(i);
%     t4(i) = s4(i)/a(i);
% end
% 
% tow = 1:0.00005:m;
% 
% n1 = zeros(1, length(tow));
% n2 = zeros(1, length(tow));
% n3 = zeros(1, length(tow));
% n4 = zeros(1, length(tow));
% 
% 
% for i = 1:length(tow)
%     n1(i) = sum(t1 < tow(i));
%     n2(i) = sum(t2 < tow(i));
%     n3(i) = sum(t3 < tow(i));
%     n4(i) = sum(t4 < tow(i));
% end
% 
% % Normalize
% n1 = n1 / length(s1);
% n2 = n2 / length(s1);
% n3 = n3 / length(s1);
% n4 = n4 / length(s1);
% 
% % Plot
% figure;
% plot(tow, n1, 'r', 'LineWidth', 2); hold on;
% plot(tow, n2, 'k', 'LineWidth', 2);
% plot(tow, n3, 'g', 'LineWidth', 2);
% plot(tow, n4, 'b', 'LineWidth', 2);
% hold off;
% 
% xlabel('\tau');
% ylabel('\rho_s(\tau)');
% legend({'FOM-C','SOM-C', 'MNM-C', 'TOM-C'}, 'Location','northwest');
% title('Comparison Plot');
% grid on;
