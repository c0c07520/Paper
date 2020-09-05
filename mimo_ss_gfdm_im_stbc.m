%-------------STBC-MIMO-GFDM---------------%
clear all;
%****************************************************************
fprintf( 'STBC-MIMO-GFDM 2x2 \n\n') ;
%****************************************************************
global Num_Sub_Per_Group
global Num_Bit_Per_Group
global Num_Space_Bit_Per_Group
global Walsh_Metric
global index
%---------IM-SS setting------------------%
Num_Sub_Per_Group = 4;
index = 2;
Num_Space_Bit_Per_Group = floor(log2(Num_Sub_Per_Group));
Num_Bit_Per_Group = Num_Space_Bit_Per_Group + index;
%---------Spreading Code-----------------%
%Walsh_Metric = hadamard(Num_Sub_Per_Group);
Chirp = exp(1i*pi/Num_Sub_Per_Group*[0:Num_Sub_Per_Group-1].^2);
Walsh_Metric = zeros(Num_Sub_Per_Group);
for ii = 0:Num_Sub_Per_Group-1
    Walsh_Metric(:,ii+1) = circshift(Chirp.',ii).';
end

Nt=2;
Nr=2;   
Nsc=128;  %Number subcarrier
Num_Group = Nsc / Num_Sub_Per_Group;
Ncp=20;
%frame_len = Nsc+Ncp;  % 一個frame的長度  Nsc + Ncp = 128 + 32 = 160
framelen=96;   % 96
Num_Path=6;  %Num path
num_bits=Num_Bit_Per_Group*Num_Group*framelen; % total bits SC * frame  128x96 = 12288
num_symbol=Num_Sub_Per_Group*Num_Group*framelen/Nt; % QPSK  

nc_symbol=num_symbol/Nsc;  % 3072 / 128 = 24 Num carrier symbol
Group=2;  % number groups
num_STBC_Groups = 1;
K_Nsc = Group*Nsc; % block * Nsc
%sub_symbol_len = K*nc_symbol*frame_len;
%%% GFDM setting%%%
p.pulse='rc';
p.a=0.2;
p.K =128;
p.L=2;
p.Kon=p.K;
p.M=nc_symbol*2;
p.mu=2;
p.oQAM = 0;

t = linspace(-p.M/2, p.M/2, p.M*p.K+1); t = t(1:end-1); t = t';
g1 = (sinc(t) .* cos(pi*p.a*t) ./ (1-4*p.a*p.a*t.*t));
g = fftshift(g1);
g(p.K+1:p.K:end) = 0;
g2 = g / sqrt(sum(g.*g));

g = g2(round(1:(p.K/p.L):length(g2)));
G = fft(g);

K = p.K; M = p.M;
N = p.M*p.K;
L = length(G) / M;
F = 1/sqrt(Nsc)*exp(-j*2*pi/Nsc*[0:Nsc-1]'*[0:Nsc-1]); % FFT matrix
SNR_dB=[0:2:20];
iter=10;

BER1=[];
for i=1:length(SNR_dB)
    tic;
    sigma2=10^(-SNR_dB(i)/10);
    BER1_num=0;
    %--Monte carlo
    for t=1:iter
        bits_source = randi([0,1], 1,num_bits);  % 1 x 12288
        vblast_out = reshape(bits_source,Nt,length(bits_source)/Nt);  % 2 x 6144 (Nt =2)
        for mo=1:Nt
            Random_Bit(:,:,mo) = reshape(vblast_out(mo,:),Num_Bit_Per_Group,length(vblast_out(mo,:))/Num_Bit_Per_Group);
            Trans_Sym(:,:,mo) = Generate_Trans_Sym_GFDM(Random_Bit(:,:,mo));
            symbols(mo,:) = reshape(Trans_Sym(:,:,mo),1,num_symbol);
            %symbols_mo(mo,:) = modulation(vblast_out(mo,:),index); %mapping 2 x 3072
        end
        for vt=1:num_STBC_Groups
            STBC_out(:,:,vt) = symbols(2*vt-1:2*vt,:);  % 2 x 3072
        end
        %%%STBC %%%
        symbol_Tx_out = [];
        for va=1:num_STBC_Groups
            STBC1 = STBC_out(1,:,va); %[ x1 x1 .... ]  1x3072
            STBC2 = STBC_out(2,:,va); %[ x2 x2 .... ]  1x3072
            temp1 = conj(STBC1); % [ x1* x1* .... ]
            temp2 = -conj(STBC2); % [ -x2* -x2* ....]
            STBC1_out = [STBC1; temp2]; %[ x1 ... ; -x2*...]
            STBC2_out = [STBC2; temp1]; %[ x2 ... ;  x1...*]
            STBC1_symbol_out(1,:,va) = STBC1_out(:).' ; % P/S  串轉並 1 x 6144  [ x1 -x2* .... ]
            STBC2_symbol_out(1,:,va) = STBC2_out(:).' ; % P/S  串轉並 1 x 6144  [ x2  x1* .... ]
            symbol_sp1 = reshape(STBC1_symbol_out(1,:,va), K_Nsc,length(STBC1_symbol_out(1,:,va))/K_Nsc);
            % S/P  256 * 24   (24 = number carrier symbol )
            symbol_sp2 = reshape(STBC2_symbol_out(1,:,va), K_Nsc,length(STBC2_symbol_out(1,:,va))/K_Nsc);
            %---------S/P  IFFT -------%
            for k=1:Nsc
                temp1 = symbol_sp1(2*k-1:2*k,:);  %interleaver
                symbol_ifft1 = temp1(:).' ;
                symbol_IFFT1(k,:) = symbol_ifft1;  % 128 x 24
                
                temp2 = symbol_sp2(2*k-1:2*k,:);
                symbol_ifft2 = temp2(:).' ;
                symbol_IFFT2(k,:) = symbol_ifft2.';
            end
            %d1 = reshape(symbol_IFFT1,p.K,p.M);
            x1 = do_modulate(p,symbol_IFFT1,'F');
            xcp1 =[x1(N - Ncp +1 :N) ;x1];
            
            %d2 = reshape(symbol_IFFT2,p.K,p.M);
            x2 = do_modulate(p,symbol_IFFT2,'F');
            xcp2 =[x2(N - Ncp +1 :N) ;x2];
            
            %symbol_IFFT1_out = F'*symbol_IFFT1; % IFFT 64 x 28
            %symbol_IFFT1_out = ifft(symbol_IFFT1);
            %symbol_IFFT2_out = F'*symbol_IFFT2;
            %symbol_IFFT2_out = ifft(symbol_IFFT2);
            
            %%% add cp
            %symbol_scp1_out = [symbol_IFFT1_out(Nsc -Ncp + 1:end,:); symbol_IFFT1_out]; % CP  160 x 48
            %symbol_scp2_out = [symbol_IFFT2_out(Nsc -Ncp +1:end,:); symbol_IFFT2_out];
            symbol_txx1 = reshape(xcp1,1,N+Ncp); % 1 x 3079
            symbol_txx2 = reshape(xcp2,1,N+Ncp); % 1 x 3079
            %symbol_tx1 = sqrt( Nsc/(Nsc+Ncp) )*reshape(symbol_scp1_out,1,sub_symbol_len);%P/S  1 x 7680
            %symbol_tx2 = sqrt( Nsc/(Nsc+Ncp) )*reshape(symbol_scp2_out,1,sub_symbol_len);
            
            symbol_tx_out = [symbol_txx1; symbol_txx2]; % 2 x 3079
            %symbol_Tx_out = [symbol_tx1; symbol_tx2];
        end
        %symbols_tx=symbol_Tx_out;
        symbols_Tx=symbol_tx_out;
        %---------------- generate channel tap------------------------
        h=(randn(Nt*Nr,Num_Path)+i*randn(Nt*Nr,Num_Path))/sqrt(2); % channel
        %h=h.*(ones(Nt*Nr,1)*(exp(-0.5).^[0:L-1]));
        h=h.*(ones(Nt*Nr,1)*(exp(-0.5).^[1:Num_Path])); % channel gain
        h=h./(sqrt(sum(abs(h).^2,2))*ones(1,Num_Path)); % channle impluse normalize
        
        %CL=size(h,2);
        Symbls_rx = zeros(Nr,N+Ncp);
        %symbols_rx=zeros(Nr,sub_symbol_len); % 2 x 3840
        
        for nr=1:Nr
            for nt=1:Nt
                symbols_rx_cp=conv(symbols_Tx(nt,:),h((nr-1)*Nt+nt,:));
                Symbls_rx(nr,:)=Symbls_rx(nr,:)+symbols_rx_cp(1:end-Num_Path+1);
            end
            %--------noise--------%
            noise_power=(sum(abs(Symbls_rx(nr,:)).^2)/length(Symbls_rx(nr,:)))*sigma2;
            noise = sqrt(noise_power)*(randn(size(Symbls_rx(nr,:)))+j*randn(size(Symbls_rx(nr,:))));
            %noise = sqrt(sigma2)*(randn(size(symbols_rx(nr,:)))+j*randn(size(symbols_rx(nr,:))))/sqrt(2);
            %--------receive signal --------%
            Symbls_rx(nr,:)=Symbls_rx(nr,:)+noise;  % 2 x 3079
            %--------remove CP | FFT--------%
            %symbols_cp = reshape(Symbls_rx(nr,:),Nsc+Ncp,sub_symbol_len/(Nsc+Ncp));%S/P  80 x 48
            symbols_ifft=Symbls_rx(nr,Ncp+1:N + Ncp).'; %remove CP  1x3072
            %symobls_re_cp = reshape(symbols_ifft,Nsc,24); %128x24
            %symbols_freq(:,:,nr)=fft(symbols_ifft); % FFT
            symbols_freq(:,:,nr)=do_demodulate(p,symbols_ifft,'MF'); %128x24x2
        end
        
        
        
        h_temp=[h zeros(Nt*Nr,Nsc-Num_Path)]; % padding 補0  4 x 128
        
        for a = 1:Nsc
            H_c = [];
            r_c = [];
            for nr=1:Nr
                Hh=[];
                for nt=1:Nt
                    h_f2 = fft(h_temp(nt+(nr-1)*Nt,:)).';
                    h_f = sqrt(Nsc)*F(a,:)*h_temp(nt+(nr-1)*Nt,:).';  % H frequency domain
                    Hh = [Hh h_f];
                end
                H_c = [H_c; Hh]; %[ h11 h12; h21 h22];
                r_c = [r_c; symbols_freq(a,:,nr)]; %[ y11 y12 ; y21 y22];
            end
            H = H_c;
            r = r_c;
            
            
            H_qr = [];
            for nR=1:Nr
                Hm = [];
                for nT=1:Nt/2
                    temp_qr = [conj(H(nR,2*nT)) -1*conj(H(nR,2*nT-1))];
                    Hm = [Hm temp_qr];
                end
                temp1_qr = [H(nR,:); Hm];
                H_qr = [H_qr; temp1_qr];
            end
            y_qr = [];
            for s_qr=1:nc_symbol
                r_qr = [];
                for it=1:Nr
                    rm(:,:,it) = [r(it,2*s_qr-1); conj(r(it,2*s_qr))];
                    r_qr = [r_qr; rm(:,:,it)];
                end
                y_qr = [y_qr r_qr];
            end
            %-------------detection algorithm-------------------%
            G1=inv(H_qr'*H_qr+sigma2*eye(Nt))*H_qr'; %MMSE
            S1_out = G1*y_qr;
            for n=1:Nt
                S1_est(a,:,n) = S1_out(n,:);
            end
        end
        bits_vblast_est1 = [];
        for est=1:Nt
            symbol_stbc = reshape(S1_est(:,:,est),Num_Sub_Per_Group, num_symbol/Num_Sub_Per_Group);
            bit_est1 = ML_detection_GFDM(symbol_stbc);
            %symbol_stbc1 = reshape(S1_est(:,:,est), 1, num_symbol); % P/S
            %bit_est1 = demodulation(symbol_stbc1,index);
            bit_est2 = reshape(bit_est1,1,num_bits/Nt);
            bits_vblast_est1 = [bits_vblast_est1; bit_est2];
        end
        bits_est1 = bits_vblast_est1(:).'; % P/s
        BER1_num = BER1_num + sum( bits_est1 ~= bits_source );
    end
    BER1 = [BER1 BER1_num/(iter*num_bits)];
    %disp(['SNR=', num2str(SNR_dB(i)) '; BER1 = ',num2str(BER1(i))]);
    toc;    
end
figure (1)
semilogy(SNR_dB,BER1,'b-s','LineWidth',1.5)
hold on
%axis([0 20 10^-5 10^0])
grid on
xlabel('SNR (dB)')
ylabel('BER')

% semilogy(SNR_dB,BER1,'k-s','LineWidth',1.5); %
% xlabel('SNR dB');   ylabel('BER'); %title('Nt=4,Nr=2');
% grid on;
% hold on;
% time_end=datestr(now);
% disp(time_begin);
% disp(time_end);