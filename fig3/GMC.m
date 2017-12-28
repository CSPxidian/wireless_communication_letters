
% *********************** GMC ******************  
clc
clear
disp(['The simulation of MaxCorrentyopy  in alpha noise '])
disp([' alpha£¬GMC'])
%% 
disp([])

%
%delta is  the kernel size
%Gsnrdb  is GSNR
%epsilon is a constant
%nSample is the samples number
%N is the  number of antennas 
%alpha is a characteristic exponent
%gama  is a dispersion coefficient
%beta is a symmetrical parameter
%txSignal is 
%attn_2 is the 
%rxSignal is the 
%noise is the alpha noise
%RX() is the received signals expect noise in H1
%RX_2() is the received signals expect noise in H0
%kxi is the sparse vecto
%d is the search direction
%g is the Eq.(11)
%J1 is the Eq.(8) in H1
%J2 is the Eq.(8) in H0
%% 
%*********************  **************************
delta = sqrt(0.5);  %% 
Gsnrdb = [-20 : 5: -10];  
epsilon = 0;
t = 0:0.25:64;  
nSample = length(t);         % samples in signal
N = 5;

%****************** alpha noise  *******************
alpha = 1.5 ;
gama = 1;
beta = 0;
a=0;

%***************** thresh ******************************
thresh_max = 1/(sqrt(2*pi)) ;
thresh = (0:0.0001: thresh_max) ;
%%
%******************  calculation  *********************************
L = numel(Gsnrdb);
nloop = 1000;
hWait = waitbar(0,'please wait...');
for ii = 1:L
    d = 0;
    d_2 = 0 ;
    for jj = 1: nloop
                
        %**********  transmit signal   *******************************************
        infoSignal = cos(2 * pi * t /10)' ;
        txSignal = infoSignal ;
        
        %***************  attn ***********************
        attn_2 = sqrt( 10.^(Gsnrdb(ii)/10) * nSample / sum(abs(txSignal).^2)); 
                
        %**********  receive signal   *******************************************
        rxSignal = txSignal  *  attn_2;      % H1
        rxSignal_2 = 0 ;            % H0
        
        %***************** calculate noise ************************************
        RX =  zeros(nSample,N); 
        RX_2 =  zeros(nSample,N);
        for k = 1:1:N
            %*************   noise  ************************
            noise = alpha_channel(nSample,alpha,gama,beta,a) ;  % alpha 
            %*************  Signal + noise ***********************
            RX(:,k) = rxSignal + noise' ;  % H1
            RX_2(:,k) = rxSignal_2 + noise' ;  % H0
        end %for k

        g = zeros(5,1000);
        d = zeros(5,1000); 
        kxi = ones(1,5);  
        g_2 = zeros(5,1000);
        d_2 = zeros(5,1000);
        kxi_2 = ones(1,5);
       %% 
       %********* calculate correntropy in H1 ******************************
        for  k = 1:1:1000
            % ****************** calculate   d ************************
            e = txSignal -  RX * kxi';            %Eq.(4)
            h =  lambda * u * zeta .*  ( exp(-lambda .* ( abs(e) ).^u) .* ( abs(e) ).^(u-1) .* sign(e) );  %Eq.(5)
            g(:,k) =  - RX' * h ; 
            if  sum(kxi==0) == 5 
                d(:,k) = -min ( g(:,k) , 0);      %Eq.(10)
            elseif  sum(kxi>0)  > 1
                if k == 1                         %Eq.(11)
                    d(:,k) = -g(:,k); 
                elseif k > 1                      %Eq.(11)
                    a1 = g(:,k)' *  g(:,k-1) / ( g(:,k-1)' * g(:,k-1)  ) * g(:,k-1);
                    y = g(:,k) - a1;
                    beta = g(:,k)' *  y / ( g(:,k-1)' * g(:,k-1) );
                    d(:,k) = -g(:,k)  + beta * d(:,k-1) ;      
                else
                end
            else
            end 
            
            % ****************** judgment of d ************************                    
            if  abs( g(:,k)' * d(:,k) ) <= 0.01
                break ; 
            else 
                % ****************** calculate  step size   ************************
                for jjj = 0:1:50
                    StepSize = 0.5.^jjj;                   % step size in Eq.(9)
                    val_1 = kxi' + StepSize * d(:,k) ;     % Eq.(9)
                    val_2 = - ( sum( exp(-lambda * ( abs( txSignal -  RX * val_1 ) ).^u)) - sum( exp(-lambda * ( abs(txSignal -  RX * kxi') ).^u)) ) ;   %Eq.(12)
                    if ( sum(val_1 >= 0) == 5 ) && ( val_2 <= 0 )
                        break;
                    else
                    end
                end   
            end
            % ****************** calculate  kxi  ************************
            kxi =  kxi + StepSize * d(:,k)';               %Eq.(9) 
        end % for  k = 0 :1:100
                                     
        J1 = zeta * sum( exp ( -lambda .* abs( txSignal - RX * kxi' ).^u ) ) ;     %Eq.(8)
        % *************  Comparison of J1 and threshold  **************************   
        Num_1 = Num_1 + ( J1 > thresh );         %  H1
        
       %% 
       %********* calculate correntropy in H0 *************************************
        for  i = 1:1:1000    
            % ****************** calculate  d  ************************
            e = txSignal - RX_2 * kxi_2' ;                   %Eq.(4)
            h =  lambda * zeta *  u .* exp(-lambda .* ( abs(e) ).^u) .* ( abs( e ) ).^(u-1) .* sign(e) ;
            g_2(:,i) =  -( RX_2' * h );              
            if  sum(kxi_2==0) == 5 
                d_2(:,i) = -min ( g_2(:,i) , 0);             %Eq.(10)
            elseif sum(kxi_2>0)  > 1   
                if i == 1                                    %Eq.(11)
                    d_2(:,i) = -g_2(:,i); 
                elseif i > 1                                 %Eq.(11)
                    a2 =  g_2(:,k)' *  g_2(:,k-1)  / ( g_2(:,k-1)' * g_2(:,k-1))* g_2(:,k-1);
                    y_2 = g_2(:,k) - a2;
                    beta = g_2(:,i)' * y_2 / ( g_2(:,i-1)' * g_2(:,i-1) );
                    d_2(:,i) = -g_2(:,i)  + beta * d_2(:,i-1) ;
                else
                end
            else
            end 
                       
            % ****************** judgment of d ************************                
            if  abs( g_2(:,i)' * d_2(:,i) ) < 0.01
                break ;
            else           
                % ****************** calculate  step size  ************************
                for jjj = 0:1:20
                    StepSize = 0.5.^jjj;
                    val_1 = kxi_2' + StepSize * d_2(:,i) ;    %Eq.(9)
                    val_2 = -( ( exp(-lambda * ( abs( txSignal -  RX_2 * val_1) ).^u))) - ( exp(-lambda * ( abs( txSignal -  RX_2 * kxi_2') ).^u)) + 0.1 * alpha.^2 * d_2(:,i)' * d_2(:,i) ;

                    if ( sum(val_1 >= 0) == 64 ) && (sum(val_2 <= 0)  == 64)    
                        break;
                    else
                    end
                end
            end
            % ****************** calculate  kxi  ************************
            kxi_2 =  kxi_2 + StepSize * d_2(:,i)';            %Eq.(9)
           
        end % 
        J2 = zeta * sum( exp ( -lambda .* abs( txSignal -  RX_2 * kxi_2' ).^u ) );   %Eq.(8)
        % ************* Comparison of J2 and threshold **************************   
        Num_2 = Num_2 + (  J2 > thresh );   %  H0    
            
    end 
    pd(:,ii) = Num_1;  
    pf(:,ii) = Num_2;                          
    waitbar(ii/L,hWait);   
end   %end of for  ii = 1:L
pd = pd/nloop;   
pf = pf/nloop;   
close(hWait);

%% 
%***************  output result  ***************
figure(1);
plot(pf(:,1),pd(:,1),'s-r', ...
     pf(:,2),pd(:,2),'+-r', ...
     pf(1:2:length(thresh),3),pd(1:2:length(thresh),3),'o-r', 'LineWidth',1.5); 

%
xlabel('False alarm probability,{\it P_f}');
ylabel('Detection probability,{\it P_d }');
    
leg = legend('GSNR = -20dB,{\it N} = 256','GSNR = -15dB,{\it N} = 256','GSNR = -10dB,{\it N} = 256') ;         
set(leg,'Location','SouthEast');
set(gcf, 'position', [400 400 400 300]);
grid on
hold on



%**********************  end of file  ***************************
 



