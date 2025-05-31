%Stock Market Volatility Analysis
clear;clc;close all;
data = readtable('AAPL.csv');

if iscell(data.Price)
    prices = str2double(strrep(data.Price, ',', ''));
else
    prices = data.Price;
end

prices = flip(prices);
returns = diff(log(prices));

S0 = prices(end);
mu = mean(returns);
sigma = std(returns);
T = 1;
N = 252;
dt = T/N;
nSimulations = 100;

simulated_prices = zeros(N+1, nSimulations);
simulated_prices(1, :) = S0;

for i = 1:nSimulations
    for t = 2:N+1
        dW = sqrt(dt) * randn;
        simulated_prices(t, i) = simulated_prices(t-1, i) * exp((mu - 0.5*sigma^2)*dt + sigma*dW);
    end
end

figure;
plot(simulated_prices(:, 1:10), 'LineWidth', 1);  % Plot only 10 simulations
title('Monte Carlo Simulation of Stock Price Paths');
xlabel('Days');
ylabel('Price');
grid on;

[h_ttest, p_ttest] = ttest(returns);
fprintf('\nT-Test Result: H = %d, P = %.4f\n', h_ttest, p_ttest);

[h_ks, p_ks] = kstest(returns);
fprintf('K-S Test Result: H = %d, P = %.4f\n', h_ks, p_ks);

figure;
histogram(returns, 'Normalization', 'pdf');
hold on;
x = linspace(min(returns), max(returns), 100);
pdf_fit = normpdf(x, mu, sigma);
plot(x, pdf_fit, 'r', 'LineWidth', 2);
title('Histogram of Returns with Normal Distribution Fit');
xlabel('Log Return');
ylabel('Probability Density');
legend('Returns', 'Normal Fit');
grid on;

fprintf('\nMean Return: %.6f\n', mu);
fprintf('Standard Deviation: %.6f\n', sigma);