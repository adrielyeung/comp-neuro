%% 1. BCM Rule
%
simtime = 10/1e-3;  % in ms
eta = 1e-6;  % in ms^{-1}
y0 = 10;
tau = 50;  % in ms

x = [20, 0; 0, 20];  % all possible inputs

w = zeros(2, simtime+1);
w(:, 1) = [0.5, 0.5];  % initial weight = 1 (equal) - no competition
theta = zeros(1, simtime+1);  % initial condition theta = 0
y = zeros(1, simtime+1);

for t = 1:simtime
    % Pick new input
    pick = rand;  % uniform rand no. between 0 and 1
    if pick < 0.5
        xpick = x(1,:);
    else
        xpick = x(2,:);
    end
    
    % New output
    y(t) = dot(w(:,t), xpick);
    % Update weight and threshold
    %theta(t+1) = mean(y(1:t))^2/y0;
    theta(t+1) = theta(t) - theta(t)/tau + y(t)^2/(y0*tau);
    delta_w = eta*xpick*y(t)*(y(t)-theta(t));
    w(:,t+1) = w(:,t) + delta_w.';
    for wi = 1:length(w(:,1))
        if w(wi,t+1) < 0
            w(wi,t+1) = 0;
        end
    end
end
%}

% Plot
figure(1)
subplot(3,1,1);
hold on
for wi = 1:length(w(:,1))
    plot(0:simtime, w(wi,:), 'DisplayName', sprintf('Neuron %u',wi))
end
xlim([0, simtime])
xlabel("Time (ms)")
ylabel("Synaptic weight \it{w}")
hold off
legend()
%}

subplot(3,1,2);
plot(0:simtime, theta)
xlabel('Time (ms)')
ylabel('\theta')
xlim([0, simtime])

subplot(3,1,3);
plot(0:simtime, y)
xlabel('Time (ms)')
ylabel('\it{y}')

%% 2. STDP
%
period = 1000; % time between pairing, in ms
pair = 60; % no. of pairings
A_pl = 1;
A_mn = 1;
tau_pl = 10; % in ms
tau_mn = 20; % in ms

t_lag = -50:5:50;  % in ms

w_change = zeros(1, length(t_lag));

% Simulate
for tl_ind = 1:length(t_lag)
    simtime_2 = period*(pair - 1) + abs(t_lag(tl_ind)) + 1; % +1 because simulation from t=1 to t=T+1
    t_pre = 1; % Set initial presynaptic time, in ms
    t_post = t_pre - t_lag(tl_ind);  % Find initial postsynaptic time, in ms
    if t_post < 1
        t_post = t_post + period;  % Require postsynaptic time to be positive
    end

    x_pre = zeros(1, simtime_2);
    y_post = zeros(1, simtime_2);

    for t = 1:simtime_2
        x_pre(t+1) = x_pre(t) + 1/tau_pl*(-x_pre(t) + (t - t_pre == 0));
        y_post(t+1) = y_post(t) + 1/tau_mn*(-y_post(t) + (t - t_post == 0));
        delta_w = A_pl*x_pre(t)*(t - t_post == 0) - A_mn*y_post(t)*(t - t_pre == 0);
        w_change(tl_ind) = w_change(tl_ind) + delta_w;
        if t - t_pre == 0
            t_pre = t_pre + period;  % Find next presynaptic time after 1 period
        end
        if t - t_post == 0
            t_post = t_post + period;  % Find next postsynaptic time after 1 period
        end
    end
end

% Plot
figure(2)
plot(t_lag, w_change, 'b-')
xlabel('Time lag (pre - post) (ms)')
ylabel('Weight change \Delta \it{w_{ij}}')
