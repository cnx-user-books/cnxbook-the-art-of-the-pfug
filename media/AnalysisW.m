function AnalysisW(I)

% Plots comparison of winp x I using computation and analysis

Ivec = linspace(2,30,15);
tau = 20; % ms
vth = -54; % mV
vrest = -70; 
%% Make plot with winp from computation
for k = 1:length(Ivec)
    I = Ivec(k);
    w(k) = compW(I);  % w is smallest weight where v reaches vth
end

plot(Ivec,w,'k')
hold on

%% Make plots with winp from analysis
Wan = (vth-vrest)*(1 - exp(-Ivec/tau));
plot(Ivec,Wan,'--r')
%% Note chosen value of winp
a = find(Ivec==20);
plot(Ivec(a),w(a),'sk','MarkerFaceColor','k')

xlabel('I (ms)', 'FontSize', 20)
ylabel('Winp (mV)', 'FontSize', 20)
%title('Winp from computation and analysis x I', 'FontSize', 26)
legend('Computation', 'Analysis','location','southeast')
hold off