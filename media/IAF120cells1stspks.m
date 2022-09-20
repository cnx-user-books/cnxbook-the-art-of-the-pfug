function IAF120cells1stspks

% (Same program as IAF120cellsSTDP, different plots)

% Plot first spike times as a function of time
% (Can also plot conductance, voltage)

%% Set Parameters

taum = 20;    % membrane time constant (ms)
tauE = 5;     % exc. time constant (ms)
tref = 5;     % refractory period (ms)
vrest = -70;  % resting potential (mV)
vreset = -60; % reset potential (mV)
vth = -54;    % threshold voltage (mV)
vE = 0;       % exc reversal potential (mV)

Nplc = 120;   % number of place fields
winp = 10;     % synaptic input
ISI = 20;     % interspike interval for external input (ms)
tinp = 0;     % time of next input (conductance) spike
tplc = 100;   % time spent in each place field (ms)
dt = 1;       % timestep (ms)

% Form place cell weight matrix
winit = 0.5;    % initial input weight
ee = winit * ones(Nplc-1,1);
W = diag(ee,1) + diag(ee,-1);
W(1,end) = winit;  W(end,1) = winit;

wmax = 5;
Aplus = .08 * wmax;
Aminus = .084 * wmax;
tauplus = 20;
tauminus = 20;

%% Run the simulation

% Constants
gc1 = (2*tauE - dt) / (2*tauE + dt);
gc2 = 2/(2*tauE + dt);

nlaps = 20;
laptime = Nplc * tplc;  % sim time of each lap (ms)
Nstps = ceil(nlaps*laptime/dt);    % total num of steps

mon = [];
gEmon = zeros( length(mon), Nstps );  % matrix of exc conductance in time-steps
vmon = vrest + gEmon;      % matrix of voltage (mV) in time-steps

wmon = zeros( length(mon), Nstps );

gE = zeros(Nplc,1);  v = vrest+gE;   % conductances and voltages of each cell

Tspk = -1000 * ones(Nplc,1);  % time of last volt spikes of cell 1, cell 2 (ms)
spkcells = [];

spkmon = 2;   % cell to monitor first spk of each lap
firstspk = zeros(1, nlaps) - 10;   % time of first spk of cell spkmon at each lap

t = 0;  
j = 0;   % only necessary for monitoring

for lap=1:nlaps    % laps
	pfield = 0;   % reset to 0 instead of mod
    nspkmon = 0;   % reset num of spks of spkmon
    
while pfield < Nplc   % place fields of one lap (120)
    
    pfield = pfield+1;
    tfld = 0;
    
    while tfld < tplc   % time steps of each place field
        tfld = tfld + dt;
        t = t + dt;
        j = j + 1;
        
        % Update conductances
        gEnext = gc1 * gE;   % at every time-step
        
        % Inputs from external (grid cells)
        if t >= tinp      % when there is an input spike, reset gE
            gEnext(pfield) = gEnext(pfield) + gc2 * winp;
            tinp = tinp + ISI;
        end
        
        % Inputs from internal (place cells)
        if ~isempty(spkcells)
            gEnext = gEnext + gc2 * sum(W(:,spkcells),2);  % 2 sums along rows
        end
      
        % Update voltages
        
        % 1) update volt's (ODE) at every time-step
        num = (2*taum/dt - 1 - gE) .* v + 2*vrest + (gEnext+gE)*vE;
        den = 2*taum/dt + 1 + gEnext;
        vnext = num ./ den;
        
        % 2) find: refrac cells, reset volt's
        refcells = find(-Tspk+t < tref);     % = cell 1, cell 2, both, or empty
        
        if ~isempty(refcells)     % if there are ref'y cells, reset their volt's
            vnext(refcells) = vreset;    % reset voltages
        end
        
        % 3) find: spiked cells, make spike, reset gE(2), adjust weights
        spkcells = find(vnext >= vth);

        if ~isempty(spkcells)	% if there are spiking cells, reset gE
            vnext(spkcells) = vreset;
            
            % LTP, spkcells = post
            post = spkcells; pre = [spkcells - 1, mod(spkcells, Nplc) + 1];
            if any(pre == 0)
                pre( find(pre==0) ) = Nplc;
            end
            deltat = t - Tspk(pre);
            winc = Aplus * exp( -deltat / tauplus );
            W(post, pre) = min( W(post, pre) + winc', wmax );
            
            % LTD, spkcells = pre
            pre = spkcells; post = [spkcells - 1, mod(spkcells, Nplc) + 1];
            if any(post == 0)
                post( find(post==0) ) = Nplc;
            end
            deltat = t - Tspk(post);
            wdec = Aminus * exp( -deltat / tauminus );
            W(post, pre) = max( W(post, pre) - wdec, 0 );
            
            if any(spkcells==spkmon)
                nspkmon = nspkmon+1;
                if nspkmon==1
                    firstspk(lap) = mod(t, Nplc*tplc);
                    firstspk(lap)
                end
            end
            
            Tspk(spkcells) = t;
        end
        
%         gEmon(:,j) = gEnext(mon);
%         vmon(:,j) = vnext(mon);
        wmon(1, j) = W(1,2);
        wmon(2, j) = W(2,1);
        
        gE = gEnext;
        v = vnext;
     
    end  % while tfld < tplc
    
end  % while pfield < Nplc
        
end  % while lap <= nlaps


%% Plots
% Move 'return' to include desired plots

% Plot first spikes of cell spkmon
figure(4)
plot(1:nlaps, firstspk, '*m')
title(['The time at which cell ', int2str(spkmon), ' fires each lap'])
ylabel('Time (ms)')
xlabel('Lap Number')
axis tight

return

% Plot the conductance of first cell in mod
figure(1)
subplot(2,1,1)

plot(tvec,gEmon(1,:))
title(['Cell ', int2str(min(mon)), ', ', int2str(nlaps), ' laps'])
ylabel('Conductance')
axis tight
 xlim([(min(mon)-1)*tplc (max(mon)+.5)*tplc])

% Plot the voltage of first cell in mod
subplot(2,1,2)
plot(tvec,vmon(1,:))
% xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis tight
 xlim([(min(mon)-1)*tplc (max(mon)+.5)*tplc])

% Plot the threshold line
hold on
plot(tvec,-54*ones(size(tvec)),'--r')
hold off



% Plot conductance of second cell in mod
figure(2)
subplot(2,1,1)
plot(tvec,gEmon(2,:))
title(['Cell ', int2str(max(mon)), ', ', int2str(nlaps), ' laps'])
ylabel('Conductance')
axis tight
 xlim([(min(mon)-1)*tplc (max(mon)+.5)*tplc])

% Plot the voltage of second cell in mod
subplot(2,1,2)
plot(tvec,vmon(2,:))
xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis tight
 xlim([(min(mon)-1)*tplc (max(mon)+.5)*tplc])

% Plot the threshold line
hold on
plot(tvec,-54*ones(size(tvec)),'--r')
hold off