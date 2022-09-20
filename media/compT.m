function CfirstT = compT(I)

% Plot T x winp using computational method

%% Set Parameters

tau = 20;    % membrane time constant (ms)
vrest = -70;  % resting potential (mV), original= -70
vth = -54;    % threshold voltage (mV)


Nplc = 100;   % number of place fields/cells
winp = linspace(10.2,20,Nplc)';     % p cell VOLTAGE input (mV), original=2,20,Nplc
tinp = 0;     % time of next p cell VOLTAGE input spike
simT = 300;   % total sim time
dt = 1;       % timestep (ms)

%% Run the simulation

% Constants
Nstps = ceil(simT/dt);    % total num of steps

% Monitor cells
mon = [];
vmon = vrest + zeros(length(mon), Nstps);    % row of voltage (mV) in time-steps
Nmon = length(vmon);

ee = ones(Nplc,1);
v = vrest*ee;   % column, v at t, each cell

SPK = zeros(Nplc, 1);
Ctimes = zeros(Nplc, 1);   % cell x time

j = 0;   % only necessary for monitoring

t = 0;


while t <= simT
    
    j = j + 1;
    
    % Update all cells' voltages
    num = 2 * vrest * dt + v * (2*tau - 1);
    den = 2*tau + 1;
    Cvnext = num ./ den;  % Cvnext is volts of t
    
    % External input, update v
    if t >= tinp
        Cvnext = Cvnext + winp;
        tinp = tinp + I;
    end
   
    
    if any(Cvnext >= vth)
        S=find(Cvnext>=vth);
        indS = find(Ctimes(S) == 0);
        SPK(S) = 1;
        if t==0
            Ctimes(S(indS)) = 1;
        else Ctimes(S(indS)) = t;
        end
    end
    
    if Nmon
        vmon(:,j) = Cvnext(mon);
    end
    v = Cvnext;
    v(S) = vth+5;
    t = t + dt;
end  % t<=simT

% plot(winp, Ctimes)
% title(['Computation-First spike time as a function of Winp, I = ', int2str(I)],'fontsize',20)
% xlabel('Winp (mV)','fontsize',20)
% ylabel('Time of First Spike (ms)','fontsize',20)

if Nmon
    plot(1:length(vmon),vmon)
end
CfirstT = Ctimes;
return
