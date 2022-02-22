function [y_out,w_out] = LCMV(Fs,x,source_freq,interfere_freq,angle,M,d)
    
    J = input('幾個Tap? ');
    max_freq = max(max(source_freq),max(interfere_freq));
    a_source = zeros(J*M,length(source_freq));
    a_interfere = zeros(J*M,length(interfere_freq));
    Cs=[];
    for k=1:length(source_freq)
        for j=1:J
            for m=1:M
                if d==1
                    a_source(M*(j-1)+m,k) = exp(-1j*2*pi*(source_freq(k)/Fs)*((m-1)*sin(angle(1))+(j-1)));
                elseif d==2
                    a_source(M*(j-1)+m,k) = exp(-1j*pi*(source_freq(k)/max_freq)*((m-1)*sin(angle(1))+(2*max_freq*(j-1)/Fs)));
                end
            end
        end
        Cs=[Cs a_source(:,k)];
    end
    Ci=[];
    for k=1:length(interfere_freq)
        for j=1:J
            for m=1:M
                if d==1
                    a_interfere(M*(j-1)+m,k) = exp(-1j*2*pi*(interfere_freq(k)/Fs)*((m-1)*sin(angle(2))+(j-1)));
                elseif d==2
                    a_interfere(M*(j-1)+m,k) = exp(-1j*pi*(source_freq(k)/max_freq)*((m-1)*sin(angle(2))+(2*max_freq*(j-1)/Fs)));
                end
            end
        end
        Ci=[Ci a_interfere(:,k)];
    end
    gs = ones(length(source_freq),1);
    gi = zeros(length(interfere_freq),1);
    ifFrost = input('用Frost? (Y/N) ','s');
    if ifFrost=='Y'
        [y1,w1] = Frost(Cs,Ci,gs,gi,x,M,J,Fs,source_freq,interfere_freq,d);
    end
    ifOpt = input('Optimization? (Y/N) ','s');
    if ifOpt=='Y'
        C = [Cs Ci];
        g = [gs;gi];
        freq = [source_freq interfere_freq];
        [y2,w2] = LCMV_opt(C,g,x,M,J,Fs,freq,d);
    end
    
    if ifFrost && ifOpt
        choose = input('要使用哪一個algorithm的output與weights? 1.Frost 2.LCMV optimization ');
        if choose==1
            y_out = y1;
            w_out = w1;
        else
            y_out = y2;
            w_out = w2;
        end
    end
end

