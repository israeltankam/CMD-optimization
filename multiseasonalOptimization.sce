function Optim(varargin) 
    clc;
    // argument in in order : contactType,contactRate,inference,aggregation
    // contactType  = "massAction"; for a beta*U*S mass action infection function
    //                "frequency"; for frequency-dependant beta*U*S/(U+V) function
    //                "inference"; for vector inferences in infection function
    // contactRate = real number; for the contact rate beta_p = beta_v
    // If the following variable exist, then beta_p = beta_v = k_v = k_b give the inoculation and acquisition rates
    // inference = real number; coefficient of vector inference
    // aggregation = coefficient of vector aggregation
    if (varargin(1) == "massAction" || varargin(1) == "frequency") then
        if length(varargin) < 2 then
            disp("Arguments missing in function Optim");abort;
        else
            beta_p = varargin(2)
            beta_v = varargin(2)
            cost_vector = varargin(3)
            numberOfSeasons = varargin(4)
        end        
    end
    if varargin(1) == "inference"  then
        if length(varargin) < 4 then
            disp("Arguments missing in function Optim");abort;
        else
            kp = varargin(2)
            kv = varargin(2)
            mp = varargin(3)
            mv = varargin(4)
            cost_vector = varargin(5)
            numberOfSeasons = varargin(6)
        end
    end
    if numberOfSeasons<2 then
        disp("There must be at least 2 seasons");
    end
    rho = 0.2
    m = 20
    mu = 1/20
    omega = 0.1882
    alpha = [0.01 0.007 0.003 0.0008]
    gammaq = [0.0006988 0.0006314 0.5899016 0.8241926 0.8621587]
    K = 10000 // Density of plants per hectare
    T = 300
    step = 0.01
    t0 = 0
    Area = 1 //hectare
    selling_price = 380000 // XOF per ton
    maintenance_price = 300000 // XOF per hectare
    aH = 0.95
    aL = 0.05
    theta = [0 0.3586 0.6414  0]
    // State variables are in plants/hectares and whiteflies/hectare
    
    //******************************************** Initial conditions ****************************************
    U0 = 6000; V0=100;
    //********************************** Hyperparameter lambda ***********************************************
    lambda_values = logspace(-5,6,12)
    //lambda_values = [10^5]
    
    //********************************************************************************************************
    function ag = aggregation(x,a) // If the aggregation of vectors is used, this function might be useful
    if x>0 then
        ag = x^a
    else
        ag = 0
    end
    endfunction
    //****************************************************************************************************
    function vert = vertical(a,b) // In case of vertical transmission
    if a+b> 0 then
        vert = a/b
    else
        vert = 0
    end
    endfunction
    //****************************************************************************************************
    //###### Mass Action
    if varargin(1) == "massAction" then
        function yprim= season(t,y)
        yprim=[-beta_p*y(30)*y(1)
               beta_p*y(30)*y(1) - mu*y(2)
               mu*y(2)-alpha(1)*y(3) - gammaq(1)*y(3)
               alpha(1)*(y(3)-y(4)) - gammaq(2)*y(4)
               alpha(1)*(y(4)-y(5)) - gammaq(3)*y(5)
               alpha(1)*(y(5)-y(6)) - gammaq(4)*y(6)
               alpha(1)*y(6) - gammaq(5)*y(7)              // y(7) = I^{(5)} of HS
               -beta_p*y(30)*y(8)
               beta_p*y(30)*y(8) - mu*y(9)
               mu*y(9)-alpha(2)*y(10) - gammaq(1)*y(10)
               alpha(2)*(y(10)-y(11)) - gammaq(2)*y(11)
               alpha(2)*(y(11)-y(12)) - gammaq(3)*y(12)
               alpha(2)*(y(12)-y(13)) - gammaq(4)*y(13)
               alpha(2)*y(13) - gammaq(5)*y(14)              // y(14) = I^{(5)} of S
               -beta_p*y(30)*y(15)
               beta_p*y(30)*y(15) - mu*y(16)
               mu*y(16)-alpha(3)*y(17) - gammaq(1)*y(17)
               alpha(3)*(y(17)-y(18)) - gammaq(2)*y(18)
               alpha(3)*(y(18)-y(19)) - gammaq(3)*y(19)
               alpha(3)*(y(19)-y(20)) - gammaq(4)*y(20)
               alpha(3)*y(20) - gammaq(5)*y(21)              // y(21) = I^{(5)} of R
               -beta_p*y(30)*y(22)
               beta_p*y(30)*y(22) - mu*y(23)
               mu*y(23)-alpha(4)*y(24) - gammaq(1)*y(24)
               alpha(4)*(y(24)-y(25)) - gammaq(2)*y(25)
               alpha(4)*(y(25)-y(26)) - gammaq(3)*y(26)
               alpha(4)*(y(26)-y(27)) - gammaq(4)*y(27)
               alpha(4)*y(27) - gammaq(5)*y(28)              // y(28) = I^{(5)} of HR
               rho*vertical(y(1)+y(2)+y(8)+y(9)+y(15)+y(16)+y(22)+y(23),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30))*(1-(y(29)+y(30))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30)))) - beta_v*y(29)*(y(3)+y(4)+y(5)+y(6)+y(7)+y(10)+y(11)+y(12)+y(13)+y(14)+y(17)+y(18)+y(19)+y(20)+y(21)+y(24)+y(25)+y(26)+y(27)+y(28)) - omega*y(29)
               rho*vertical(y(3)+y(4)+y(5)+y(6)+y(7)+y(10)+y(11)+y(12)+y(13)+y(14)+y(17)+y(18)+y(19)+y(20)+y(21)+y(24)+y(25)+y(26)+y(27)+y(28),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30))*(1-(y(29)+y(30))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30)))) + beta_v*y(29)*(y(3)+y(4)+y(5)+y(6)+y(7)+y(10)+y(11)+y(12)+y(13)+y(14)+y(17)+y(18)+y(19)+y(20)+y(21)+y(24)+y(25)+y(26)+y(27)+y(28)) - omega*y(30)
              ]
         endfunction
    end
    //##### Frequency
    if varargin(1) == "frequency" then
        function yprim= season(t,y)
        yprim=[-beta_p*y(30)*y(1)/(y(29)+y(30))
               beta_p*y(30)*y(1)/(y(29)+y(30)) - mu*y(2)
               mu*y(2)-alpha(1)*y(3) - gammaq(1)*y(3)
               alpha(1)*(y(3)-y(4)) - gammaq(2)*y(4)
               alpha(1)*(y(4)-y(5)) - gammaq(3)*y(5)
               alpha(1)*(y(5)-y(6)) - gammaq(4)*y(6)
               alpha(1)*y(6) - gammaq(5)*y(7)              // y(7) = I^{(5)} of HS
               -beta_p*y(30)*y(8)/(y(29)+y(30))
               beta_p*y(30)*y(8)/(y(29)+y(30)) - mu*y(9)
               mu*y(9)-alpha(2)*y(10) - gammaq(1)*y(10)
               alpha(2)*(y(10)-y(11)) - gammaq(2)*y(11)
               alpha(2)*(y(11)-y(12)) - gammaq(3)*y(12)
               alpha(2)*(y(12)-y(13)) - gammaq(4)*y(13)
               alpha(2)*y(13) - gammaq(5)*y(14)              // y(14) = I^{(5)} of S
               -beta_p*y(30)*y(15)/(y(29)+y(30))
               beta_p*y(30)*y(15)/(y(29)+y(30)) - mu*y(16)
               mu*y(16)-alpha(3)*y(17) - gammaq(1)*y(17)
               alpha(3)*(y(17)-y(18)) - gammaq(2)*y(18)
               alpha(3)*(y(18)-y(19)) - gammaq(3)*y(19)
               alpha(3)*(y(19)-y(20)) - gammaq(4)*y(20)
               alpha(3)*y(20) - gammaq(5)*y(21)              // y(21) = I^{(5)} of R
               -beta_p*y(30)*y(22)/(y(29)+y(30))
               beta_p*y(30)*y(22)/(y(29)+y(30)) - mu*y(23)
               mu*y(23)-alpha(4)*y(24) - gammaq(1)*y(24)
               alpha(4)*(y(24)-y(25)) - gammaq(2)*y(25)
               alpha(4)*(y(25)-y(26)) - gammaq(3)*y(26)
               alpha(4)*(y(26)-y(27)) - gammaq(4)*y(27)
               alpha(4)*y(27) - gammaq(5)*y(28)              // y(28) = I^{(5)} of HR
               rho*vertical(y(1)+y(2)+y(8)+y(9)+y(15)+y(16)+y(22)+y(23),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30))*(1-(y(29)+y(30))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30)))) - beta_v*y(29)*(y(3)+y(4)+y(5)+y(6)+y(7)+y(10)+y(11)+y(12)+y(13)+y(14)+y(17)+y(18)+y(19)+y(20)+y(21)+y(24)+y(25)+y(26)+y(27)+y(28))/sum(y(1:28)) - omega*y(29)
               rho*vertical(y(3)+y(4)+y(5)+y(6)+y(7)+y(10)+y(11)+y(12)+y(13)+y(14)+y(17)+y(18)+y(19)+y(20)+y(21)+y(24)+y(25)+y(26)+y(27)+y(28),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30))*(1-(y(29)+y(30))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12)+y(13)+y(14)+y(15)+y(16)+y(17)+y(18)+y(19)+y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(27)+y(28))*(y(29)+y(30)))) + beta_v*y(29)*(y(3)+y(4)+y(5)+y(6)+y(7)+y(10)+y(11)+y(12)+y(13)+y(14)+y(17)+y(18)+y(19)+y(20)+y(21)+y(24)+y(25)+y(26)+y(27)+y(28))/sum(y(1:28)) - omega*y(30)
              ]
         endfunction
    end
    //************************************************************************************************************************************
    function varargout = deal(a)
    varargout = list()
    for i = 1:length(a)
        varargout(i) = a(i)
    end
    endfunction
    //************************************************************************************************************************************
    function yy = tuber_yield(finalValues)
        ss = sum(finalValues(1:28));
        index = 0;
        for i=1:length(finalValues)-2
            if pmodulo(i,7)==0 then
                k = 5
            elseif pmodulo(i,7)<=3 then
                k=1
            else
                k = pmodulo(i,7)-2
            end
             index = index+k*finalValues(i)
        end
        mean_index_of_severity = index/ss
        //disp(ss)
        //disp(mean_index_of_severity)
        yy = (-4.03*mean_index_of_severity + 31.63)*ss/K //To calibrate the yield to the number of plants remaining in the field 
    end
    //************************************************************************************************************************************
    function p = profit_season(yield,st) // st = O in strategy 0 and 1 in strategy1
        cost = K*theta*cost_vector'
        //p = yield
        p = Area*selling_price*yield - cost*st - Area*maintenance_price
    endfunction
    function strategy = choose_strategy(n,lambda) // n = length of strategy
        strategy = []
        solution = []
        H1_0 = aH*theta(1)*K
        H2_0 = aH*theta(2)*K
        H3_0 = aH*theta(3)*K
        H4_0 = aH*theta(4)*K
        L1_0 = aL*theta(1)*K
        L2_0 = aL*theta(2)*K
        L3_0 = aL*theta(3)*K
        L4_0 = aL*theta(4)*K
        I11_0 = (1-aL-aH)*theta(1)*K
        I21_0 = (1-aL-aH)*theta(2)*K
        I31_0 = (1-aL-aH)*theta(3)*K
        I41_0 = (1-aL-aH)*theta(4)*K
        I12_0 = 0
        I22_0 = 0
        I32_0 = 0
        I42_0 = 0
        I13_0 = 0
        I23_0 = 0
        I33_0 = 0
        I43_0 = 0
        I14_0 = 0
        I24_0 = 0
        I34_0 = 0
        I44_0 = 0
        I15_0 = 0
        I25_0 = 0
        I35_0 = 0
        I45_0 = 0
        t_init = t0
        t = t_init+step:step:T
        Sol = ode([H1_0;L1_0;I11_0;I12_0;I13_0;I14_0;I15_0;H2_0;L2_0;I21_0;I22_0;I23_0;I24_0;I25_0;H3_0;L3_0;I31_0;I32_0;I33_0;I34_0;I35_0;H4_0;L4_0;I41_0;I42_0;I43_0;I44_0;I45_0;U0;V0],t_init,t,season);
            for k=1:n
                     t_init = (k-1)*T
                     t = t_init+step:step:k*T
                    // strategy = 0 : ratoon
                    ll = length(Sol(1,:));
                    finalValues = round(Sol(:,ll))
                    ss = sum(finalValues(1:28))
                    [H1_0,L1_0,I11_0,I12_0,I13_0,I14_0,I15_0,H2_0,L2_0,I21_0,I22_0,I23_0,I24_0,I25_0,H3_0,L3_0,I31_0,I32_0,I33_0,I34_0,I35_0,H4_0,L4_0,I41_0,I42_0,I43_0,I44_0,I45_0] = deal((K/ss)*finalValues(1:28))
                    sol0 = ode([H1_0;L1_0;I11_0;I12_0;I13_0;I14_0;I15_0;H2_0;L2_0;I21_0;I22_0;I23_0;I24_0;I25_0;H3_0;L3_0;I31_0;I32_0;I33_0;I34_0;I35_0;H4_0;L4_0;I41_0;I42_0;I43_0;I44_0;I45_0;U0;V0],t_init+T,t+T,season);
                    yield0 = tuber_yield(round(sol0(:,length(sol0(1,:)))))
                    regulator0 = sum(sol0(2:7,length(sol0(1,:)))+sol0(9:14,length(sol0(1,:)))+sol0(16:21,length(sol0(1,:)))+sol0(23:28,length(sol0(1,:))))
                   // strategy = 1 : from market
                    H1_0 = aH*theta(1)*K
                    H2_0 = aH*theta(2)*K
                    H3_0 = aH*theta(3)*K
                    H4_0 = aH*theta(4)*K
                    L1_0 = aL*theta(1)*K
                    L2_0 = aL*theta(2)*K
                    L3_0 = aL*theta(3)*K

                    L4_0 = aL*theta(4)*K
                    I11_0 = (1-aL-aH)*theta(1)*K
                    I21_0 = (1-aL-aH)*theta(2)*K
                    I31_0 = (1-aL-aH)*theta(3)*K
                    I41_0 = (1-aL-aH)*theta(4)*K
                    I12_0 = 0
                    I22_0 = 0
                    I32_0 = 0
                    I42_0 = 0
                    I13_0 = 0
                    I23_0 = 0
                    I33_0 = 0
                    I43_0 = 0
                    I14_0 = 0
                    I24_0 = 0
                    I34_0 = 0
                    I44_0 = 0
                    I15_0 = 0
                    I25_0 = 0
                    I35_0 = 0
                    I45_0 = 0
                    sol1 = ode([H1_0;L1_0;I11_0;I12_0;I13_0;I14_0;I15_0;H2_0;L2_0;I21_0;I22_0;I23_0;I24_0;I25_0;H3_0;L3_0;I31_0;I32_0;I33_0;I34_0;I35_0;H4_0;L4_0;I41_0;I42_0;I43_0;I44_0;I45_0;U0;V0],t_init+T,t+T,season);
                    yield1 = tuber_yield(round(sol1(:,length(sol1(1,:)))))
                    regulator1 = sum(sol1(2:7,length(sol1(1,:)))+sol1(9:14,length(sol1(1,:)))+sol1(16:21,length(sol1(1,:)))+sol1(23:28,length(sol1(1,:)))) 
                    if profit_season(yield0,0) - lambda*regulator0 > profit_season(yield1,1) - lambda*regulator1 then
                        strategy = [strategy 0]
                        Sol = sol0
                    else
                        strategy = [strategy 1]
                        Sol = sol1
                    end
                    
        end
    endfunction
//************************************************************************************************************************************
    function yy = tuber_yield(finalValues)
        ss = sum(finalValues(1:28));
        index = 0;
        for i=1:length(finalValues)-2
            if pmodulo(i,7)==0 then
                k = 5
            elseif pmodulo(i,7)<=3 then
                k=1
            else
                k = pmodulo(i,7)-2
            end
             index = index+k*finalValues(i)
        end
        mean_index_of_severity = index/ss
        //disp(ss)
        //disp(mean_index_of_severity)
        yy = (-4.03*mean_index_of_severity + 31.63)*ss/K //To calibrate the yield to the number of plants remaining in the field 
    end
    //************************************************************************************************************************************
    function [YY,yield] = Modele(strategy)
        yield = 0
        t_init = t0
        t_tot =[]
        solution = []
        H1_0 = aH*theta(1)*K
        H2_0 = aH*theta(2)*K
        H3_0 = aH*theta(3)*K
        H4_0 = aH*theta(4)*K
        L1_0 = aL*theta(1)*K
        L2_0 = aL*theta(2)*K
        L3_0 = aL*theta(3)*K
        L4_0 = aL*theta(4)*K
        I11_0 = (1-aL-aH)*theta(1)*K
        I21_0 = (1-aL-aH)*theta(2)*K
        I31_0 = (1-aL-aH)*theta(3)*K
        I41_0 = (1-aL-aH)*theta(4)*K
        I12_0 = 0
        I22_0 = 0
        I32_0 = 0
        I42_0 = 0
        I13_0 = 0
        I23_0 = 0
        I33_0 = 0
        I43_0 = 0
        I14_0 = 0
        I24_0 = 0
        I34_0 = 0
        I44_0 = 0
        I15_0 = 0
        I25_0 = 0
        I35_0 = 0
        I45_0 = 0
        for k=1:length(strategy)+1
            t = t_init+step:step:k*T
            Sol = ode([H1_0;L1_0;I11_0;I12_0;I13_0;I14_0;I15_0;H2_0;L2_0;I21_0;I22_0;I23_0;I24_0;I25_0;H3_0;L3_0;I31_0;I32_0;I33_0;I34_0;I35_0;H4_0;L4_0;I41_0;I42_0;I43_0;I44_0;I45_0;U0;V0],t_init,t,season);
            yield = yield + tuber_yield(round(Sol(:,length(Sol(1,:)))))
            solution = [solution Sol]
            t_tot =[t_tot t]
            t_init = k*T
            if k<=length(strategy)then
                if strategy(k)==0 then //from same field
                    ll = length(Sol(1,:));
                    finalValues = round(Sol(:,ll))
                    ss = sum(finalValues(1:28))
                    [H1_0,L1_0,I11_0,I12_0,I13_0,I14_0,I15_0,H2_0,L2_0,I21_0,I22_0,I23_0,I24_0,I25_0,H3_0,L3_0,I31_0,I32_0,I33_0,I34_0,I35_0,H4_0,L4_0,I41_0,I42_0,I43_0,I44_0,I45_0] = deal((K/ss)*finalValues(1:28))
                else // strategy = 1 : from market
                    H1_0 = aH*theta(1)*K
                    H2_0 = aH*theta(2)*K
                    H3_0 = aH*theta(3)*K
                    H4_0 = aH*theta(4)*K
                    L1_0 = aL*theta(1)*K
                    L2_0 = aL*theta(2)*K
                    L3_0 = aL*theta(3)*K

                    L4_0 = aL*theta(4)*K
                    I11_0 = (1-aL-aH)*theta(1)*K
                    I21_0 = (1-aL-aH)*theta(2)*K
                    I31_0 = (1-aL-aH)*theta(3)*K
                    I41_0 = (1-aL-aH)*theta(4)*K
                    I12_0 = 0
                    I22_0 = 0
                    I32_0 = 0
                    I42_0 = 0
                    I13_0 = 0
                    I23_0 = 0
                    I33_0 = 0
                    I43_0 = 0
                    I14_0 = 0
                    I24_0 = 0
                    I34_0 = 0
                    I44_0 = 0
                    I15_0 = 0
                    I25_0 = 0
                    I35_0 = 0
                    I45_0 = 0
                end 
            end
        end
        nn = length(solution(1,:))
        YY = round(solution(:,nn))
    endfunction
    //**************************************************************************************************************************************************
    function p = profit(strategy)
        cost = K*theta*cost_vector'
        [y,yield] = Modele(strategy)
        //p = yield
        p = Area*selling_price*yield - (cost + sum(cost*strategy)) - maintenance_price*(length(strategy)+1)*Area // cost in parentheses and +1 account for the first season, for which cuttings are assumed to come from the market
    endfunction
    best_strategy = []
    best_lambda = 0
    best_profit = -1e30
    for i=1:length(lambda_values)
        lambda = lambda_values(i)
        disp(lambda)
        strategy = choose_strategy(numberOfSeasons-1,lambda)
        disp(strategy)
        current_profit = profit(strategy)
        disp(current_profit)
        if current_profit > best_profit then
            best_profit = current_profit
            best_lambda = lambda
            best_strategy = strategy
        end
    end
    disp("The best strategy is: ")
    disp(best_strategy)
    disp("It yields the profit: $"+string(best_profit)+" with the hyperparameter lambda="+string(best_lambda))
endfunction
Optim("frequency",0.0066,[0 25 400 650],6) 
