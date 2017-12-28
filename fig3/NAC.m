%
% **********************基于  NAC 的频谱感知  的频谱感知 仿真**************
% 
% 输入 信号，噪声
% 输出 ROC曲线

clc
clear
disp(['The simulation of MaxCorrentyopy vs 文献 Hindawi Publishing Corporation ...'])
disp(['Mobile Information Systems ...'])
disp(['Volume 2016, Article ID 6753830, 6 pages  ...'])
disp(['http://dx.doi.org/10.1155/2016/6753830  ...'])
disp(['Publication Year 2016, ...'])
disp(['Spectrum Sensing Based on Nonparametric Autocorrelation in ...'])
disp(['Wireless Communication Systems under Alpha Stable Noise'])
disp([' alpha噪声下， NAC , Eq.(7)'])
disp([])

tic
%% 
%********************* 参数 **************************
delta = sqrt(0.5);
% Gsnrdb = [3.026 4.275 6.04 ];%3.026对应GSNR=0dB；4.275对应GSNR=3dB；6.04对应GSNR=6dB；
Gsnrdb = [-10 : 5:  0];

epsilon = 0;

%********************* 样本数 ******************************
t = 0:0.25:320;
nSample = length(t);         % samples in signal
N = 5;

%****************** alpha noise 参数  *******************
alpha = 1.5;
gama = 1;
beta = 0;
a=0;

%***************** thresh ******************************
% thresh_max = 1/(sqrt(2*pi)) ;
% thresh = (0:0.001: thresh_max )   ;

%%
%% ********* calculate  NAC 的频谱感知   *********************
L = numel(Gsnrdb);

nloop = 1000;
hWait = waitbar(0,'please wait...');
for ii = 1:L
    d = 0;
    d_2 = 0 ;
    d_3 = 0 ;
    d_4 = 0 ;
%     for jj = 1: nloop
                
        %**********  transmit signal   *******************************************
%         infoSignal = exp(-t.^2)' ; 
        infoSignal = (cos(2 * pi * t /10) + j * sin(2 * pi * t /10))' ;
        txSignal = infoSignal ;
        
        %***************  根据 Gsnrdb 计算信号的幅度  ***********************
        attn_2 = sqrt( 10.^(Gsnrdb(ii)/10) * nSample / sum(abs(txSignal).^2));  
        
        %********** 经过信道后 receive signal   *******************************************
        rxSignal = txSignal  *  attn_2;       % H1的情况
        rxSignal_2 = 0 ;            % H0的情况
        
%         %***************  calculate 噪声的大小 *********************
%         spow = sum(txSignal.^2);%/nSample;
%         attn = spow * 10 .^( -snrdb(ii)/10);
%         attn = sqrt(attn);                      % 噪声电压
        
       %%   噪声  和  统计量   

        %*************  生成噪声 ************************
        noise = alpha_channel(nSample,alpha,gama,beta,a)';  % alpha稳定分布噪声 
        x_H0 = noise ; %H0：接收信号
        x_H1 = rxSignal + noise ; %H1：接收信号
        point_delay = [2 3];
       %%  延迟为 2
        x_H1_dalay_1 = [zeros(1,point_delay(1)), x_H1(1:nSample- point_delay(1))' ];
        R_1 = 1/ nSample * sum ( x_H1 .* conj( x_H1_dalay_1 )' );

       %%  延迟为 3     
        x_H1_dalay_2 = [zeros(1,point_delay(2)), x_H1(1:nSample- point_delay(2))' ];        
        R_2 = 1/ nSample * sum ( x_H1 .* conj( x_H1_dalay_2 )' );     
       %% 计算 lamda
       clear jj;
       jj = 0;
       var_1 = 1:1:nSample';
       for s = [-0.5 : 0.2 :0.5]
           jj = jj + 1;
           var_2 = exp (-2 * j * pi * s *  var_1 ./nSample )';
         %%  延迟为 2  
           var_3(jj) = sum(( x_H1 .* conj( x_H1_dalay_1 )')' .* var_2');       
           var_4(jj) = sum(( conj(x_H1) .* x_H1_dalay_1')' .* var_2');          %Eq.(7) 
         %%  延迟为 3  
           var_5(jj) = sum(( x_H1 .* conj( x_H1_dalay_2 )')' .* var_2');         
           var_6(jj) = sum(( conj(x_H1) .* x_H1_dalay_2')' .* var_2');         %Eq.(7)             
       end % for 
      %%  延迟为 2
       S_2f =  sum( var_3.^2)/(nSample * length(s) );
       S_2f_conj =  sum ( var_3 .* var_4 )  /(nSample * length(s) );      %Eq.(7)    
      %%  延迟为 3 
       S_2f_2 = sum( var_5.^2)/(nSample * length(s) );
       S_2f_2_conj =  sum ( var_5 .* var_6 )/(nSample * length(s) );      %Eq.(7)    
            
       var_7 = real(R_1).^2 /( S_2f_conj + S_2f) + imag(R_1).^2 / ( S_2f_conj - S_2f ) ;             
       var_8 = real(R_1).^2 /( S_2f_2_conj + S_2f_2) + imag(R_1).^2 / ( S_2f_2_conj - S_2f_2 ) ; 
       lamda = real( var_7 + var_8 );
           
       %%  门限  
        pf = [0.01:0.01:1 ];
        Thresh_R_0 = chi2inv(pf,4)  ;
        Thresh_R = sort(Thresh_R_0,'descend'); 

        %% 统计判决 （蒙特卡洛）
        for i = 1: length(pf)
            clear val;
            val = find (lamda > Thresh_R(i));
            Pd(i) = length(val) / nloop ;
        end 
        
        %%  检测概率 理论值
        lambda =  4 + 2 * nSample * lamda ;
        Pdtrue = Qchipr2(4,lambda,Thresh_R,epsilon);      
                   
%     end %end of for  jj = 1: nloop

    waitbar(ii/L,hWait);   
end   %end of for  ii = 1:L

close(hWait);
toc

%% 
%***************  output result  ***************
figure(1);
plot(pf,Pdtrue,'s-r',...
     pf,Pd,'+-r','LineWidth',1.5);
xlabel('False alarm probability,{\it P_f}');
ylabel('Detection probability,{\it P_d }');

set(leg,'Location','SouthEast');
set(gcf, 'position', [400 400 400 300]);
grid on


hold on

%**********************  end of file  ***************************
 



