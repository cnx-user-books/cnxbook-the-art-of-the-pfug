function AnpeakV(I)

% Plot peak voltages

%% Set Parameters

tau = 20;    % membrane time constant (ms)
vrest = -70;  % resting potential (mV), original= -70
vth = -52;    % threshold voltage (mV)


Nplc = 100;   % number of place fields/cells
w = compW(I);
winp = linspace(w,20,Nplc)';     % p cell VOLTAGE input (mV), original=10
tinp = 0;     % time of next p cell VOLTAGE input spike
simT = 300;   % total sim time
dt = 1;       % timestep (ms)

%% Run the simulation

% Constants
Nstps = ceil(simT/dt);    % total num of steps

% Monitor cells
mon = [10];
vmon = vrest + zeros(length(mon), Nstps);    % row of voltage (mV) in time-steps
Nmon = length(vmon);

ee = ones(Nplc,1);
v = vrest*ee;   % column, v at t, each cell
Avnext = zeros(Nplc,1);   % column, v at t+1, each cell

SPK = zeros(Nplc, 1);
Atimes = zeros(Nplc, 1);   % cell x time
vspk = zeros(ceil(Nstps/I)+1,1);
tspk = zeros(ceil(Nstps/I)+1,1);
j = 0;   % only necessary for monitoring

t = 0;

Ninp = 0;   % num of inputs

while t <= simT
    
    j = j + 1;    
    
    % Update voltages
    if t < tinp   % between input spikes
        Avnext = zeros(Nplc,1);
        
        for m=0:Ninp-1
            Avnext = Avnext + winp*exp(-(t-m*I)/tau);
        end
        Avnext = Avnext + vrest;
    end
    
    % External inputs at n*I, update Avnext
    if t == 0
        Ninp = Ninp + 1;   % first inp
        Avnext = vrest + winp;
        tspk(Ninp) = 1;
        tinp = tinp + I;
        vspk(Ninp) = Avnext(mon);
    elseif t>0 && t >= tinp
        Ninp = Ninp + 1;
        Avnext = zeros(Nplc,1);
        tspk(Ninp) = tinp+1;
        for k=0:Ninp-1
            Avnext = Avnext + winp*exp(-I/tau)^k;
        end
        
        Avnext = Avnext + vrest;
        vspk(Ninp) = Avnext(mon);
        tinp = tinp + I;
    end
    
    if any(Avnext >= vth)
        S=find(Avnext>=vth);  % gets ind of Avnext (spkd cells)
        indS = find(Atimes(S) == 0);   % gets ind of S (cells where this is 1st spk)
        SPK(S) = 1;  % gives spk to spkd cells
        if t==0
            Atimes(S(indS)) = 1;  % sets spktime to t of cells where this is 1st spk
        else Atimes(S(indS)) = t;
        end
        %Avnext(S) = vth+5;   % plots spk of all spkd cells
    end
    
    if Nmon
        vmon(:,j) = Avnext(mon);
    end
    
    t = t + dt;

end  % t<=simT


%% Plot vmon
if Nmon
    figure
    plot(1:length(vmon),vmon)
    title(['Voltage, (I=', int2str(I),' ms, winp=', int2str(winp(mon)),' mV)'])
    xlabel('Time (ms)')
    xlim([0,simT])
    ylim([vrest,vth+1])
    ylabel('mV')
    hold all
    plot(tspk,vspk,':*m')   % plot peak volts
    plot(1:length(vmon),vth*ones(length(vmon),1),'--r')   % plot threshold
    legend('v(t)','v(kI)','Threshold','location','best')
    hold off

end
AfirstT = Atimes;
return
