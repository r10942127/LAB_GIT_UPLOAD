function x = Generate_signal_time_domain(Fs,s_s,s_i,ss_freq,si_freq,angle,M,noise_power)

    sz_s = length(s_s);
    iter = floor(sz_s/Fs);
    K = 2;
    sensor_out = [];
    for it=1:iter
        x0 = zeros(M,Fs);
        for k=1:3
            sstr = zeros(M,1);
            for m=1:M
                sstr(m) = exp(-1j*2*pi*(m-1)*(ss_freq(k)/Fs)*sin(angle(1)));
            end
            x0 = x0 + sstr*s_s((it-1)*Fs+1:it*Fs);
        end

        for k=1:3
            istr = zeros(M,1);
            for m=1:M
                istr(m) = exp(-1j*2*pi*(m-1)*(si_freq(k)/Fs)*sin(angle(2)));
            end
            x0 = x0 + istr*s_i((it-1)*Fs+1:it*Fs);
        end
        sensor_out = [(sensor_out) x0];
    end
    x = awgn((sensor_out),noise_power);
    %sound(x(1,:))

end

