%% 1.1 Logic gates
b = 1;
alpha = 1;
simtime = 10; % max no. of iterations

% 1.1 OR gate
x = [0, 0; 0, 1; 1, 0; 1, 1];  % p*n matrix
yt_or = [0, 1, 1, 1];  % Length p

or_err = perc(x, yt_or, b, alpha, simtime);

% 1.1 XOR gate
yt_xor = [0, 1, 1, 0];

xor_err = perc(x, yt_xor, b, alpha, simtime);

figure(1)
hold on
plot(1:simtime, or_err, 'b-')  % line diff from answer due to diff in Heaviside
plot(1:simtime, xor_err, 'r-')
xlim([0,simtime])
xlabel("Number of iterations")
ylabel("Squared error")
legend('OR', 'XOR')
hold off
%}

%% 1.2 Random dataset
%
simtime_rnd = 50;
p_try = [50, 100, 200];  % no. of patterns to try
p_try_str = ["p = 50", "p = 100", "p = 200"];

figure(2)
hold on

% Generate random x (input patterns) and y_target (labels)
for p = p_try
    x_rnd = randi(2, p, 100) - 1; % -1 so that the values are 0 or 1
    yt_rnd = randi(2, 1, p) - 1;

    rnd_err = perc(x_rnd, yt_rnd, b, alpha, simtime_rnd);

    plot(1:simtime_rnd, rnd_err, '-')
end
xlim([0,simtime_rnd])
xlabel("Number of iterations")
ylabel("Squared error")
legend(p_try_str)
hold off
%}

%% The Perceptron
%

function E = perc(x_inp, y_tar, bb, alal, simt)
pn = size(x_inp);
w = zeros(1,pn(2));

y = zeros(1,pn(1));
E = zeros(1,simt);

for t = 1:simt
    for mu = 1:pn(1)
        y(mu) = ((w * x_inp(mu,:)' - bb) >= 0);
        % note diff in def of Heavside cf answers
        w = w + alal*(y_tar(mu) - y(mu))*x_inp(mu,:);
        E(t) = E(t) + ((y_tar(mu) - y(mu)).^2);
    end
    if E(t) < 1e-5  % tolerance
        break
    end
end
end