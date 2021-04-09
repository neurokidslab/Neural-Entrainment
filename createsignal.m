% n : length of the signal
% fs : sampling rate
% frq : frequency of the cosine
% t0 : position of the signal
function signal = createsignal(n, fs, frq, t0)

% generating peak
signal = zeros (1, n);
t = 1/fs*(1:round(1/frq/2*fs));
y = cos(2*pi*t*frq-pi/2);
t0 = t0-floor(length(y)/2);
signal(t0:t0+length(y)-1) = y;

end
