%% Temporal Difference Learning
%cue (stimulus) at second 5, reward at second 20

Trials=100; %number of trials
Time=20;    %total time
rewTime=20; %reward time
cueTime=5; %start cue
endCueTime=Time; %end cue
n=endCueTime-cueTime+1; %cue duration (5s - 20s, so 16 elements)

X= eye(n);  % n*n identity (the stimulus go from cueTime to endCueTime so "slides" across like identity matrix
X=[zeros(n,cueTime-1), X, zeros(n,Time-endCueTime)];

V=zeros(Time,Trials);  % prediction of total future reward expected from time t at each trial
w = zeros(n,1); %weights: sum up everything until reward
r = zeros(Time,Trials); %reward
norew_tr = 50;  % omission of reward at trial norew_tr
r(rewTime,1:norew_tr)=1;  % Reward at rewTime
delta = zeros(Time, Trials); %prediction error


gamma= 1;
alpha= 0.6;

%t=time, i=trial
for i=1:Trials
    V(:,i)= w.'*X(:,:); %value function
    V_tplus1 = [V(2:end,:); zeros(1,Trials)];  % value func at next timestep
    delta(:,i)= r(:,i) + gamma*V_tplus1(:,i) - V(:,i);%prediction error
    w= w + alpha*(X(:,:)*delta(:,i)); %weights
end


%% Plot 

%Plot prediction error
figure
surf(delta')
ylabel('trials')
xlabel('time')
zlabel('prediction error')

%Plot value function 
figure
surf(V)
xlabel('trials')
ylabel('time')
zlabel('V')