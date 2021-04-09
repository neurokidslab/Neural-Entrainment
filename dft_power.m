function [data_power, freq] = dft_power(Y, srate, cndidx, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.DataField         = 'data';
P.fmin              = 0.1;
P.fmax              = 25;
P.dohilbert         = 0; 

% Optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dft_power: Non recognized inputs')
end



%% ------------------------------------------------------------------------
%% Extract signal envelopp with hilbert filtering
if P.dohilbert==1
    for i = 1:size(Y, 3)
        Y(:, :, i) = (envelope(Y(:, :, i)',30, 'peak'))';
    end
end

%% ------------------------------------------------------------------------
%% compute DFT
fprintf('DFT\n')
Ne = size(Y,1);
Ns = size(Y,2);
Y = Y - repmat(mean(Y,2), [1 Ns 1]);
f_res = srate / Ns;             % resolution (Hz)
freq  = (0:floor(Ns/2))*f_res;  % frequencies
idx_freq = logical(freq>=P.fmin & freq<=P.fmax);
freq = freq(idx_freq);
Ydft = fft(Y,Ns,2);
Ydft = Ydft(:,1:floor(Ns/2)+1,:);
Ydft = Ydft(:,idx_freq,:);


%% ------------------------------------------------------------------------
%% Power Spectrum
theCNDs = unique(cndidx);
Ncnd = length(unique(cndidx));
Nfrq = length(freq);
data_power = nan(Ne,Nfrq,Ncnd);
for i=1:Ncnd
    cndi = theCNDs(i);
    data_power(:,:,i) = 2 * abs(nanmean( Ydft(:,:,cndidx==cndi), 3)).^2 / (Ns.^2); % Power Spectrum
end

end

