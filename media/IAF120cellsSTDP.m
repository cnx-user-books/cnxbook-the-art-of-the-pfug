function IAF120cellsSTDP

% (Same program as IAF120cells1tspks, different plot)

% Plot weights of synapses as a function of time

%% Set Parameters

taum = 20;    % membrane time constant (ms)
tauE = 5;     % exc. time constant (ms)
tref = 5;     % refractory period (ms)
vrest = -70;  % resting potential (mV)
vreset = -60; % reset potential (mV)
vth = -54;    % threshold voltage (mV)
vE = 0;       % exc reversal potential (mV)

Nplc = 120;
winp = 10;     % synaptic input
ISI = 20;     % interspike interval for inputs (ms)
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

mon = [1 2];   % cells to monitor
gEmon = zeros( length(mon), Nstps );  % matrix of exc conductance in time-steps
vmon = vrest + gEmon;      % matrix of voltage (mV) in time-steps

wmon = zeros( length(mon), Nstps );

% conductances and voltages of each cell at t
gE = zeros(Nplc,1);
v = vrest+gE;

Tspk = -1000 * ones(Nplc,1);  % time of last volt spikes of cell 1, cell 2 (ms)
spkcells = [];

t = 0;  
j = 0;   % only necessary for monitoring

for lap=1:nlaps    % laps
	lap
	pfield = 0;   % reset to 0 instead of mod

while pfield < Nplc   % place fields of one lap (120)
    
    pfield = pfield+1;
    tfld = 0;
    
    while tfld < tplc   % time steps of each place field
        tfld = tfld + dt;
        t = t + dt;
        j = j + 1;
        
        % Update conductances
        gEnext = gc1 * gE;   % at every time-step
        
        % External input (grid cells)
        if t >= tinp      % when there is an input spike, reset gE
            gEnext(pfield) = gEnext(pfield) + gc2 * winp;
            tinp = tinp + ISI;
        end
        
        % Internal input (place cells)
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
        
        % 3) find: spiked cells, make spike, adjust weights
        spkcells = find(vnext >= vth);

        if ~isempty(spkcells)	% if there are spiking cells, reset vnext
            vnext(spkcells) = vreset;
            
            % STDP
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
            
            Tspk(spkcells) = t;
        end
        
        gEmon(:,j) = gEnext(mon);
        vmon(:,j) = vnext(mon);
        wmon(1, j) = W(1,2);
        wmon(2, j) = W(2,1);
        
        gE = gEnext;
        v = vnext;
     
    end  % while tfld < tplc
    
end  % while pfield < Nplc
        
end  % while lap <= nlaps

tvec = (1:size(wmon,2))*dt;

%% Plot weight
figure(3)
plot(tvec,wmon(:,:))
%title('Weight')
xlabel('Time (ms)')
ylabel('Weight (mV)')
legend('w21','w12','location','best')  % change labels as needed
axis tight

return
