function App

% Plot   w12 as a function of n1*I
% and    n2 as a function of w12 

%% Set Parameters
I = 20;
tau = 20;    % membrane time constant (ms)
vth = -54;    % threshold voltage (mV)
vrest = -70;   % resting voltage (mV)
wext = linspace(10.2,20,100)';     % external input weights (mV)

%% Calculate n1
x = exp(-I/tau);
a = 1-( ((vth-vrest)./wext)*(1-x) );
n1vec = ceil(-tau/I * log(a));   % vector of min num of wext inp spks necessary for Cell 1 to fire
w12vec = (vth-vrest)*(1-exp(-n1vec*I./tau));

% Calculate min w12 for Cell 2 activity for chosen n1
n1 = 2;   % Chosen value for n1
MINw12 = (vth-vrest)*(1-exp(-n1*I./tau));   % min w12 for chosen n1
w12 = linspace(MINw12,20,100);

% Plot w12 as a function of n1*I
figure
plot(n1vec*I,w12vec)
hold on
ind=find(n1vec==n1);
plot(n1vec(ind)*I,w12vec(ind),'sk','MarkerFaceColor','k')
hold off
xlabel('n1*I (ms)','fontsize',16)
ylabel('w12 (mV)','fontsize',16)

%% Calculate n2
b = 1-( ((vth-vrest)./w12)*(1-x) );
n2 = ceil(-tau/I * log(b));   % min num of w12 inp spks necessary for Cell 2 to fire

% Plot n2 as a function of w12
figure
plot(w12,n2)
hold on
ind2=find(n2==1)
plot(w12(ind2(1)),n2(ind2(1)),'sk','MarkerFaceColor','k')
hold off
xlabel('w12 (mV)','fontsize',16)
ylabel('n2', 'fontsize', 16)
c = max(n2);
xlim([MINw12 20])
ylim([0 c+1])
return