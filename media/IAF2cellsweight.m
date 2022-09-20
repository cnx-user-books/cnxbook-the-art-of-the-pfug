function IAF2cellsweight

% make it for 120 cells
% diag

%% Set Parameters

taum = 20;    % membrane time constant (ms)
tauE = 5;     % exc. time constant (ms)
tref = 5;     % refractory period (ms)
vrest = -70;  % resting potential (mV)
vreset = -60; % reset potential (mV)
vth = -54;    % threshold voltage (mV)
vE = 0;       % exc reversal potential (mV)

winp = 10;     % synaptic input
ISI = 20;     % interspike interval for inputs (ms)
tinp = 0;     % time of next input (conductance) spike
tplc = 100;   % time spent in each place field (ms)
w21 = 0.5;    % input from cell 1 to cell 2
dt = .1;       % timestep (ms)

%% Run the simulation - through place fields of cells 1 and 2

% Constants
gc1 = (2*tauE - dt) / (2*tauE + dt);
gc2 = 2/(2*tauE + dt);

tfin = 2 * tplc;  % total sim time (ms)
Nstps = ceil(tfin/dt);    % total num of steps
% pfield = floor(tfin/tplc);   % which place field you are in

gE = zeros(2,Nstps);  % matrix of exc conductance in time-steps
v = vrest + gE;      % matrix of voltage (mV) in time-steps

Tspk = [-1000 -1000];  % time of last volt spikes of cell 1, cell 2 (ms)

t = 0;

W = [0 .5; .5 0];  % w12 = w21 = 0.5

for j = 1:Nstps
    ispk = 0;
   
    if t >= tinp
	if t < tplc, ispk = 1;             % input spike place field 1
	elseif t < 2 * tplc, ispk = 2;     % input spike place field 2
	end
    end

    % Update conductances
    gE(:,j+1) = gc1 * gE(:,j);   % at every time-step

    if ispk      % when there is an input spike, resets gE
	gE(ispk,j+1) = gE(ispk,j+1) + gc2 * winp;
        tinp = tinp + ISI;
    end

        
    % Update voltages

    % 1) update volt's (ODE) at every time-step
    num = (2*taum/dt - 1 - gE(:,j)) .* v(:,j) + 2*vrest + (gE(:,j+1)+gE(:,j))*vE;
    den = 2*taum/dt + 1 + gE(:,j+1);
    v(:,j+1) = num ./ den;
      
    % 2) find: refrac cells, reset volt's
    refcells = find(-Tspk+t < tref);     % = cell 1, cell 2, both, or empty
    % isempty(refcells);                    % if no ref'y cells, this is 1
	
    if ~isempty(refcells)     % if there are ref'y cells, reset their volt's
	v(refcells,j+1) = vreset;    % reset voltages
    end

    % 3) find: spiked cells, make spike, reset gE(2)
    spkcells = find(v(:,j+1) >= vth);
    % isempty(spkcells);	% if no spiking cells, this is 1

    if ~isempty(spkcells)	% if there are spiking cells, make spike, reset gE
	v(spkcells,j) = vth + 5;
	v(spkcells,j+1) = vreset;
	Tspk(spkcells) = t;
        gE(:,j+1) = gE(:,j) + gc2 * sum(W(:,spkcells),2);  % 2 sums along rows
    end
	
    t = t+dt;

end


% Plot the conductance
figure(1)
subplot(2,1,1)
tvec = (1:Nstps+1)*dt;
plot(tvec,gE(1,:),'b')
hold on
plot(tvec,gE(2,:),'g')
hold off
ylabel('Conductance')
axis tight

% Plot the voltage
subplot(2,1,2)
plot(tvec,v(1,:),'b')
hold on
plot(tvec,v(2,:),'g')
plot(tvec,-54*ones(size(tvec)),'--r')
legend('Cell 1','Cell 2','Threshold','location','NorthOutside')
hold off

xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis tight

return
