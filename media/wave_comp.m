t= 0:.01:1;

amps = [1 -4];
freqs = [4 1];

[tg fg] = meshgrid(t,freqs);

y = zeros(size(t));

for i = 1:numel(amps)
  y = y + amps(i)*sin(freqs(i)*2*pi*t);
end

%amps'*sin((fg.*tg)*2*pi);

plot(t,y)