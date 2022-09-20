function AnalysisT(I)

% Plot num of spks needed for activity as a function of winp from
% computation and analysis

%% Set Parameters

tau = 20;    % membrane time constant (ms)
vth = -54;    % threshold voltage (mV)
vrest = -70;   % resting voltage (mV)
Nplc = 100;   % num of cells for winp
winp = linspace(10.2,20,Nplc)';     % p cell VOLTAGE input (mV), original=10

x = exp(-I/tau);
a = 1-( ((vth-vrest)./winp)*(1-x) );
Aninp = ceil(-tau/I * log(a));   % min num of inp spks necessary for activity

%% Plot Cninp x winp from computation
figure
CfirstT = compT(I)
Cninp= floor(CfirstT./I)+1;
plot(winp, Cninp, 'k')
hold on
%% Plot Aninp x winp from analytics
plot(winp,Aninp,'--g')
% title(['First spike times from computation and analysis x Winp, I = ', int2str(I),' (ms)'],'fontsize',20)
legend('Computation ','Analysis')
%set = (leg,'FontSize',16)
xlabel('Winp (mV)','fontsize',16)
ylabel('Num of input spikes','fontsize',16)
hold off
return
