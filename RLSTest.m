clear all;
close all;

N  = 32;

x  = rand(1,500);
b  = fir1(N-1,0.5);     % FIR system to be identified
n  = 0.01 * rand(1,500); % Observation noise signal
d  = filter(b,1,x)+n;  % Desired signal
lam = 0.9;            % RLS forgetting factor
w  = zeros(N, 1);

hMyRLS = fnRLSCreate(N, lam);

for ii = 1:1:length(x)
    if ii < N + 1
        y(ii) = [zeros(1,N-ii),x(1:ii)] * w;
        e = d(ii) - y(ii);
        [deltaW, hMyRLS] = fnRLS( hMyRLS, [zeros(1,N-ii),x(1:ii)].', e, 0 );
        w = w + deltaW;
    else
        y(ii) = x(ii-N+1:ii) * w;
        e = d(ii) - y(ii);
        [deltaW, hRLS] = fnRLS( hMyRLS, x(ii-N+1:ii).', e, 0 );
        w = w + deltaW;
    end
end

figure;
stem(b);
hold on;
plot(w,'*-');
