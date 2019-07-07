% Defining network model parameters
vt = 1;                     % Spiking threshold
tau_m = 10;                 % Membrane time constant [ms]
g_m = 1;                    % Neuron conductance
Nsig = 0.25;                % Variance amplitude of current
Nmean = 0.75;               % Mean current to neurons
tau_I = 10;                 % Time constant to filter the synaptic inputs
N_all = [500, 5000, 50000, 500000]; % Number of neurons in total
dt = 1;                     % Simulation time bin [ms]
T = 300/dt;                 % Simulation length

for N_ind = 1:length(N_all)
    N = N_all(N_ind);               % Set total number of neurons
    NE = 0.5*N;                 % Number of excitatory neurons
    NI = 0.5*N;                 % Number of inhibitory neurons
    W = 2/sqrt(N);                  % Connectivity strength

    % Initialization
    rng(100);                   % RNG seed = 100
    raster = [];                % save spike times for plotting
    raster_spk = [];            % save spike times for plotting (with perturbation)
    v_init = rand(N,1)*vt;      % save initial membrane potential
    noise = randn(N,T);         % save noise input for all timesteps
    I_oneneu = zeros(T,1);      % Step 5: synaptic current for example (excitatory) neuron

    % loop over with and without perturbation
    for i = 1:2
        v = v_init;                 % reset to initial membrane potential
        vv = zeros(N,1);            % variable that notes if v crosses the threshold
        Iback = zeros(N,1);         % building up the external current
        SP = 0;                     % recording spike times
        Ichem = zeros(N,1);         % current coming from synaptic inputs
        Iext = zeros(N,1);          % external current

        % loop over the time
        for t = 1:T
            Iback = Iback + dt/tau_I*(-Iback + noise(:,t));          % generate a colored noise for the current
            Iext = Iback/sqrt(1/(2*(tau_I/dt)))*Nsig+Nmean;         % rescaling the noise current to have the correct mean and variance

            Ichem(1:NE) = Ichem(1:NE) + dt/tau_I*(-Ichem(1:NE) + W*(sum(vv(1:NE))-vv(1:NE))-W*(sum(vv(NE+1:end)))); % current to excitatory neurons coming from the synaptic inputs
            Ichem(NE+1:end) = Ichem(NE+1:end) + dt/tau_I*(-Ichem(NE+1:end) -W*(sum(vv(NE+1:end))-vv(NE+1:end))+W*(sum(vv(1:NE)))); % current to inhibitory neurons coming from the synaptic inputs
            Itot = Iext+Ichem;
            if i == 2  % no extra spike
                I_oneneu(t) = Ichem(1);  % choosing neuron 1 (excitatory) as example
            end
            %%%%%%%%%%% To insert integrate-and-fire model here  %%%%%%%%%%%%%
            v = v + dt/tau_m*(-v + 1/g_m*Itot);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% Step 4: Extra spike %%%%%%%%%%%%%%
            if t == 200 && i == 1
                spk_ind = round(linspace(1, N, 10));  % change last element in linspace to change how many neurons to spike
                v(spk_ind) = vt;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            vv =(v>=vt);                                        % spike if voltage crosses the threshold    

            v = (1-vv).*v;                                      % reset after spike
            SP = find(vv);                                      % find the spike times
            if i == 1
                raster_spk = [raster_spk;t*ones(length(SP),1),SP];
            else
                raster=[raster;t*ones(length(SP),1),SP];            % save spike times for plotting
            end
        end
    end


    % Plot the raster output
    h = figure; hold on;
    plot(raster(:,1)*dt, raster(:,2),'.b')
    plot(raster_spk(:,1)*dt, raster_spk(:,2), '.r')
    xlim([100 300])
    xlabel('time [ms]','fontsize',20)
    ylabel('neuron index','fontsize',20)
    title(['Mean = ', num2str(Nmean), ', Sigma = ', num2str(Nsig)])
    set(gca,'fontsize',20);
    set(gca,'YDir','normal')
    hold off;

    % Step 5: Plot synaptic current
    h1 = figure(99); hold on;
    plot(1:T, I_oneneu(:,1), '-', 'DisplayName', ['N = ', num2str(N)])
    xlabel('time [ms]')
    ylabel('Synaptic current \it{I}_{chem}')
    legend()
    title('Synaptic current of example neuron, weight \propto 1/sqrt(N)')
    hold off
end