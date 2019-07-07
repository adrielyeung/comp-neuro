%% 1. Modelling the inputs
% 1.1
N = 50;  % no. of neurons to simulate
theta_test = linspace(-pi/2, pi/2, N);  % 1.2
theta_0 = 0;
c_test = 3;
epsilon_test = 0.1;
input = h_input(theta_0, theta_test, c_test, epsilon_test);  % 1.2

% 1.2 Plot result
figure(1);
plot(theta_test/pi, input, 'b-')
xlabel("\theta/\pi")
ylabel("\it{h_i^{ext}}")
grid();
title("Stimuli orientation \theta_0 " + sprintf("= %u", theta_0))

% 1.4
T_test = 0;
beta_test = 0.1;
h_test = linspace(-15, 15, 1000);
g_test = g(h_test, beta_test, T_test);  % Filter function

figure(2);
plot(h_test, g_test, 'b-')
xlabel("Input \it{h}")
ylabel("Filter \it{g}")
ylim([-0.1, 1.1]);
grid();
title("Filter function \it{g(h)}")

%% 2. Modelling the neurons
% 2.1 Test
N_i = 30;  % no. of iterations
m = zeros(1, N_i+1);

for t = 2:N_i+1  % iterate
    m(t) = one_neuron(input(1), m(t-1), 5);
end

figure(3);
plot(0:N_i, m)
title('Activity of neuron vs time')
xlabel('Time (ms)')
ylabel('Voltage')
grid();

% 2.2 Inputs: m_0, theta_0, N, N_i, epsilon, c
act_testc15 = neurons(0, 0, 50, 30, 0.9, 1.5) * 100;  % Amplify activity by 100 times

figure(4);
image([0, 30], [1, 50], act_testc15) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity over time, c = 1.5')

% 2.3 Inputs: m_0, theta_0, N, N_i, epsilon, c
act_testc12 = neurons(0, 0, 50, 30, 0.9, 1.2) * 100;  % Amplify activity by 100 times
act_testc40 = neurons(0, 0, 50, 30, 0.9, 4) * 100;  % Amplify activity by 100 times

figure(5);
image([0, 30], [1, 50], act_testc12) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity over time, c = 1.2')

figure(6);
image([0, 30], [1, 50], act_testc40) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity over time, c = 4.0')

%% 3. Modelling the network
% 3.1
J0 = 86;
J2 = 112;

J_ij = J(J0, J2, theta_test);
figure(7);
image(J_ij, 'CDataMapping', 'scaled')
xlabel('Neuron number')
ylabel('Neuron number')
title('\it{J_{ij}}')
colorbar

% 3.2 Inputs: m_0, theta_0, N, N_i, epsilon, c, J_0, J_2
act_con15 = neurons_con(0, 0, 50, 30, 0.1, 1.5, J0, J2, true) * 100;  % Amplify activity by 100 times
act_con12 = neurons_con(0, 0, 50, 30, 0.1, 1.2, J0, J2, true) * 100;  % Amplify activity by 100 times
act_con40 = neurons_con(0, 0, 50, 30, 0.1, 4, J0, J2, true) * 100;  % Amplify activity by 100 times

figure(8);
image([0, 30], [1, 50], act_con15) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with connections, c = 1.5')

figure(9);
image([0, 30], [1, 50], act_con12) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with connections, c = 1.2')

figure(10);
image([0, 30], [1, 50], act_con40) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with connections, c = 4.0')

%% 4. Apply different stimuli
% 4.1 Inputs: m_0, theta_0, N, N_i, epsilon, c, J_0, J_2
act_41_start = neurons_con(0, 0, 50, 30, 0.8, 100, J0, J2, true);  % Amplify activity by 100 times
act_41_after = neurons_con(act_41_start(:, end), 2*pi/3, 50, 500, 0.8, 100, J0, J2, true);

figure(11);
hold on
subplot(1, 2, 1);
image([0, 30], [1, 50], act_41_start*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with stimulus \theta_0 = 0')
subplot(1, 2, 2);
image([30, 530], [1, 50], act_41_after*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with stimulus \theta_0 = 2\pi/3')
hold off

% No connectivity
act_41_start_ncon = neurons_con(0, 0, 50, 30, 0.8, 100, 0, 0, true);  % Amplify activity by 100 times
act_41_after_ncon = neurons_con(act_41_start_ncon(:, end), 2*pi/3, 50, 500, 0.8, 8, 0, 0, true);

figure(12);
hold on
subplot(1, 2, 1);
image([0, 30], [1, 50], act_41_start_ncon*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with stimulus \theta_0 = 0')
subplot(1, 2, 2);
image([30, 530], [1, 50], act_41_after_ncon*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with stimulus \theta_0 = 2\pi/3')
hold off

% 4.2 Inputs: m_0, theta_0, N, N_i, epsilon, c, J_0, J_2
act_42_start = neurons_con(0, 0, 50, 30, 0.8, 100, J0, J2, true);
act_42_after = neurons_con(act_42_start(:, end), 0, 50, 30, 0.8, 100, J0, J2, false);

figure(13);
hold on
subplot(1, 2, 1);
image([0, 30], [1, 50], act_42_start*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with external input')
subplot(1, 2, 2);
image([30, 60], [1, 50], act_42_after*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity without external input')
hold off

% No connectivity after removing stimulus
act_42_start_ncon = neurons_con(0, 0, 50, 30, 0.1, 1.2, J0, J2, true);  % answer was plotted with connections initially
act_42_after_ncon = neurons_con(act_42_start_ncon(:, end), 0, 50, 30, 0.1, 1.2, 0, 0, false);

figure(14);
hold on
subplot(1, 2, 1);
image([0, 30], [1, 50], act_42_start_ncon*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity with external input')
subplot(1, 2, 2);
image([30, 60], [1, 50], act_42_after_ncon*100) %, 'CDataMapping', 'scaled')
colorbar
xlabel('Time (ms)')
ylabel('Neuron number')
title('Neuron activity without external input')
hold off

%% Function definitions
% 1.1
function out = h_input(theta0, theta, c, epsilon)
out = c*((1 - epsilon) + epsilon * cos(2*(theta - theta0)));
end

% 1.3
function g_out = g(h, beta, T)
g_out = zeros(1, length(h));
hltT = h <= T;
hltTbeta = h <= (T + 1/beta);
for i = 1:length(h)
    if hltT(i)
        g_out(i) = 0;
    elseif hltTbeta(i)
        g_out(i) = beta*(h(i) - T);
    else
        g_out(i) = 1;
    end
end
end

% 2.1
function m_after = one_neuron(h, m, tau)
m_after = (tau - 1)/tau*m + g(h, 0.1, 0)/tau;
end

% 2.2
function activity = neurons(m_0, theta_0, N, N_i, epsilon, c)
activity = zeros(N, N_i+1);  % Each row represents a neuron
activity(1, :) = m_0;  % Set initial activity
h_in = h_input(theta_0, linspace(-pi/2, pi/2, N), c, epsilon);
for neu = 1:N
    for t = 2:N_i+1
        activity(neu, t) = one_neuron(h_in(neu), activity(neu, t-1), 5);
    end
end
end

% 3.1
function J_mat = J(J_0, J_2, theta)
N_neurons = length(theta);
J_mat = zeros(N_neurons);  % N_neurons * N_neurons
for i = 1:N_neurons
    for j = i:N_neurons
        J_mat(i, j) = -J_0 + J_2*cos(2*(theta(i) - theta(j)));
        if i ~= j
            J_mat(j, i) = -J_0 + J_2*cos(2*(theta(j) - theta(i)));
            % because cos(a) = cos(-a), matrix should be symmteric
        end
    end
end
end

% 3.2 and later
function h_in_con = h_connect(theta0, theta, c, epsilon, m, J_0, J_2, ext)
% ext = whether to apply external input
% m = vector of current activity rate across all neurons
N_neurons = length(theta);
sum_con = zeros(1, N_neurons);
J_mat = J(J_0, J_2, theta);
for i = 1:N_neurons  % current neuron
    for j = 1:N_neurons  % connected neurons (including self)
        sum_con(i) = sum_con(i) + J_mat(i, j)*m(j);
    end
end
if ext
    h_in_con = sum_con + h_input(theta0, theta, c, epsilon);
else
    h_in_con = sum_con;
end
end

function act_con = neurons_con(m_0, theta_0, N, N_i, epsilon, c, J_0, J_2, ext)
% ext = whether to apply external input
act_con = zeros(N, N_i+1);  % Each row represents a neuron
act_con(:, 1) = m_0;  % Set initial activity
for t = 2:N_i+1
    h_in = h_connect(theta_0, linspace(-pi/2, pi/2, N), c, epsilon, act_con(:, t-1), J_0, J_2, ext);
    for neu = 1:N
        act_con(neu, t) = one_neuron(h_in(neu), act_con(neu, t-1), 5);
    end
end
end
