function w = compW(I)

% Outputs lowest w for input for I

%% Set Parameters

tau = 20;    % membrane time constant (ms)
vrest = -70;  % resting potential (mV), original= -70
vth = -54;    % threshold voltage (mV)


Nplc = 100;   % number of place fields/cells
winp = linspace(2,20,100)';     % p cell VOLTAGE input (mV), original=10
tinp = 0;     % time of next p cell VOLTAGE input spikesw
simT = 300;   % total sim time
dt = 1;       % timestep (ms)

%% Run the simulation

% Constants
Nstps = ceil(simT/dt);    % total num of steps

% Monitor cells
mon = [];
vmon = vrest + zeros(length(mon), Nstps);      % row of voltage (mV) in time-steps
Nmon = length(vmon);

ee = ones(Nplc,1);
v = vrest*ee;   % v at t

SPK = zeros(Nplc, 1);

j = 0;   % only necessary for monitoring

t = 0;

while t <= simT
    t = t + dt;
    j = j + 1;
    
    % Update all cells' voltages
    num = 2 * vrest * dt + v * (2*tau - 1); % v is volts of t-1
    den = 2*tau + 1;
    vnext = num ./ den;  % vnext is volts of t
    
    % External input, update v
    if t >= tinp
        vnext = vnext + winp;
        tinp = tinp + I;
    end
    
    if any(vnext >= vth)
        S=find(vnext>=vth);
        SPK(S) = 1;
    end
    
    if Nmon
        vmon(:,j) = vnext(mon);
    end
    v = vnext;
    v(S) = vth+5;
    
end  % for t<=simT

if any(SPK),
    cell = find(SPK>0,1);
    w = winp(cell);
else w = 0;
end

if Nmon
    plot(1:length(vmon),vmon)
    title(['Voltage (mV) of cell ',int2str(mon),', I = ',int2str(I),' ms'])
    xlabel('Time (ms)')
    ylabel('mV')
end

return


