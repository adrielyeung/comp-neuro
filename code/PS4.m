% Discretise the integrate-and-fire model (forward Euler method) equation
% and use it to calculate the value at each timestep:
sim_time = 100; % simulation time


delta_t = 1; % time step, unit = ms
tau_v = 10; % time constant, unit = ms
v_th = 10; % threshold voltage
R = 1; % resistance
I_arr = [9, 11, 15]; % test cases of current

t_plot = 0:delta_t:sim_time;

for i = 1:3
    I = I_arr(i);
    v = zeros(sim_time/delta_t + 1,1); % first timestep when t = 0, v = 0
    for t = 1:sim_time
       v(t+1) = v(t) + delta_t/tau_v*(-v(t) + R*I);
       if v(t+1) > v_th
          v(t+1) = 0; 
       end
    end
    subplot(3,1,i);
    plot(t_plot, v)
    hold on
    plot([0, 100], [v_th, v_th], 'b--')
    title(['Voltage of neuron vs time for current I = ', num2str(I)])
    xlabel('Time (ms)')
    ylabel('Voltage')
    ylim([0, 11])
    hold off
end

