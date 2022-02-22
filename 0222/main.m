clc;
close all;
clear all;

Fs = 1200 ;
t = linspace(0,1,Fs);
source_freq = [262 294 330]; 
interfere_freq = [349 392 440]; 
source_angle = input('Source的角度: ');
interfere_angle = input('Interfere signal的角度: ');
angle = [source_angle*pi/180 interfere_angle*pi/180];
source = 2*[exp(1j*2*pi*source_freq(1)*t) exp(1j*2*pi*source_freq(2)*t) exp(1j*2*pi*source_freq(3)*t)]  ;        % source
interfere = 10*[exp(1j*2*pi*interfere_freq(1)*t) exp(1j*2*pi*interfere_freq(2)*t) exp(1j*2*pi*interfere_freq(3)*t)] ;    % interfere
% source = 2*[cos(2*pi*source_freq(1)*t) cos(2*pi*source_freq(2)*t) cos(2*pi*source_freq(3)*t)]  ;    
% interfere = 4*[cos(2*pi*interfere_freq(1)*t) cos(2*pi*interfere_freq(2)*t) cos(2*pi*interfere_freq(3)*t)] ;   
M = input('幾個sensor? ');
noise = input('AWGN的功率: ');
choose_d = input('鄰近sensor之間的距離d (1.sample frequency的波長 2.訊號最高頻率的波長/2): ');
x = Generate_signal_freq_domain(Fs,source,interfere,max(max(source_freq),max(interfere_freq)),angle,M,noise,choose_d);  % 此function只允許1個source與1個interfere
%x = Generate_signal_time_domain(Fs,source,interfere,source_freq,interfere_freq,angle,M,0.01);
test_freq = input('是否進行頻率分析? (Y/N) ','s');
if test_freq=='Y'
    figure()
    plot(linspace(0,Fs,length(x(1,:))),abs(fft(x(1,:))))
    title('sensor輸出data的frequency response')
end
%%
test_DOA = input('是否進行DOA測試? (Y/N) ','s');
if test_DOA=='Y'
    [DOA_S,DOA_I] = find_DOA(x,180);
end
%% 
choose = input('1.LCMV 2.Response Variation: ');
if choose==1
    [y,w] = LCMV(Fs,x,source_freq,interfere_freq,angle,M,choose_d);
elseif choose==2
    [y,w] = RV(Fs,x,source_freq,interfere_freq,angle,M,choose_d);
end
%%
display = input('播出聲音? (Y/N) ','s');
if display=='Y'
    sound(real(y),Fs)
end