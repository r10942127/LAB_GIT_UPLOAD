function x = Generate_signal_freq_domain(Fs,s_s,s_i,max_freq,angle,M,noise_power,d)
    sz_s = length(s_s);
    iter = floor(sz_s/Fs);
    K = 2;
    sensor_out = [];
    for it=1:iter
        s1 = s_s(((it-1)*Fs+1):(it*Fs));
        s2 = s_i(((it-1)*Fs+1):(it*Fs));
        S1 = fft(s1);
        S2 = fft(s2);
        sw = [S1;S2];
        Aw = zeros(M,K,Fs);
        for i=1:Fs
            for m=1:M
                for k=1:K
                    if d==1
                        Aw(m,k,i)=exp(-1j*2*pi*((i-1)/Fs)*(m-1)*sin(angle(k)));
                    elseif d==2
                        Aw(m,k,i)=exp(-1j*pi*((i-1)/max_freq)*(m-1)*sin(angle(k)));
                    end
                end
            end
        end
        xw = zeros(M,Fs);
        for i=1:Fs
            xw(:,i) = squeeze(Aw(:,:,i))*sw(:,i);
        end
        x0 = zeros(M,Fs);
        for m=1:M
            x0(m,:) = ifft(xw(m,:));
        end
        sensor_out = [(sensor_out) x0];
    end
    x = awgn((sensor_out),noise_power);
    %sound(x(1,:))

end

