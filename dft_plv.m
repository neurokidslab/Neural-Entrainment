function [data_plv, freq] = dft_plv(Y, srate, cndidx, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.DataField         = 'data';
P.fmin              = 0.1;
P.fmax              = 25;

% Optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dft_plv: Non recognized inputs')
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
%% PLV
Ncnd = length(unique(cndidx));
Nfrq = length(freq);
data_plv = nan(Ne,Nfrq,Ncnd);
for i=1:Ncnd
    data_plv(:,:,i) = abs( nanmean( Ydft(:,:,cndidx==i) ./ abs(Ydft(:,:,cndidx==i)) , 3 ) );
end

end
