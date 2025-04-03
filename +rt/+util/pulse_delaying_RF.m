function [s] = pulse_delaying_RF(sig, delay_time, Fs)
    Nb = size(sig, 1);
    nb_pts_fft = 2 ^ nextpow2(Nb);

    fv = [(0:nb_pts_fft / 2 - 1), (-nb_pts_fft / 2:-1)] * Fs / nb_pts_fft;
    fv = fv.';
    TF_sig = fft(sig, nb_pts_fft);
    TF_sig_dec = TF_sig .* exp(-1i * 2 * pi .* fv .* delay_time);
    sig_delayed = real(ifft(TF_sig_dec));
    s = sig_delayed(1:Nb, :);
end
