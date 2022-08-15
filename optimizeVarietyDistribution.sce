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
            disp("Arguments missing in function Optim")
        else
            beta_p = varargin(2)
            beta_v = varargin(2)
            cost_vector = varargin(3)
        end        
    end
    if varargin(1) == "inference"  then
        if length(varargin) < 4 then
            disp("Arguments missing in function Optim")
        else
            kp = varargin(2)
            kv = varargin(2)
            mp = varargin(3)
            mv = varargin(4)
            cost_vector = varargin(5)
        end
    end
    rho = 0.2
    m = 20
    mu = 1/20
    omega = 0.1882
    alpha = [0.01 0.007 0.003 0.0008]
    gammaq = [0.000012  0.0010911 0.173672 0.329857 0.9128083]
    K = 10000 // Density of plants per hectare
    T = 300
    t0 = 0
    step = 0.01
    Area = 1 //hectare
    selling_price = 380000 // XOF per ton
    maintenance_price = 500000 // XOF per hectare
    aH = 0.95
    aL = 0.05
    // State variables are in plants/hectares and whiteflies/hectare
    
    //******************************************** Initial conditions ****************************************
    U0 = 6000; V0=200;
    //****************************************************************************************************
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
    function YY = Modele(theta) //Return all values at T
        t = t0+step:step:T
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
        Sol = ode([H1_0;L1_0;I11_0;0;0;0;0;H2_0;L2_0;I21_0;0;0;0;0;H3_0;L3_0;I31_0;0;0;0;0;H4_0;L4_0;I41_0;0;0;0;0;U0;V0],t0,t,season); 
        nn = length(Sol(1,:));
        YY = round(Sol(:,nn))
    endfunction
    //***********************************************************************************************
    function yy = tuber_yield(theta)
        final_values = Modele(theta)
        ss = sum(final_values(1:28));
        index = 0;
        for i=1:length(final_values)-2
            if pmodulo(i,7)==0 then
                k = 5
            elseif pmodulo(i,7)<=3 then
                k=1
            else
                k = pmodulo(i,7)-2
            end
             index = index+k*final_values(i)
        end
        mean_index_of_severity = index/ss
        //disp(ss)
        //disp(mean_index_of_severity)
        yy = (-4.03*mean_index_of_severity + 31.63)*ss/K //To calibrate the yield to the number of plants remaining in the field 
    end
    //**************************************************************************************************
    function pp = Profit(theta)
        pp = Area*(selling_price*tuber_yield(theta) - maintenance_price -cost_vector*theta'*K)
    endfunction
    //disp(Profit([0.25 0.25 0.25 0.25]))
    //disp(Profit([0.8 0.2 0 0]))
    //disp(Profit([0 0 0.2 0.8]))
    //Now let's write the algorithm to optimize this profit : Adaptive Random Search on Simlplexes (Tankam-Chedjou et al. (2021a))
    
    
   //************************************* Simplex ARS useful functions **********************************************************
   function D = randomDirection(n)
        M = [eye(n-1,n-1); -1*ones(1,n-1)];
        M = orth(M);
        //disp("M= "); disp(M)
        for i=1:n-1
            M(:,i)=M(:,i)/norm(M(:,i));
        end
        //disp("M= "); disp(M);
        coord = rand(1,n-1);
        D = zeros(n,1);
        for i=1:n-1
            D=D+coord(i)*M(:,i);
        end
        D = D/norm(D);
    endfunction
    //******************************************************************************************************************************
    function bool = positive(M)
        bool = %t;
        i=1;
        while i<= length(M)
            if M(i)<0 then
                bool = %f
            end
            i = i+1;
        end
    endfunction
    //*******************************************************************************************************************************
    function lambda = randomSimplexe(theta, sd) //Looking around point theta on the simplex, with standard deviation sd
        n = length(theta);
        D = (randomDirection (n))';
        rand ("normal");
        prog = sd*rand(1);
        thetaplus = theta+prog*D;
        if positive(thetaplus) then
            lambda = thetaplus;
        else
            tab_neg = [];
            tab_thetaplus_neg = [];
            for i=1:n
                if D(i)<0 then
                    tab_neg = [tab_neg D(i)]
                    tab_thetaplus_neg = [tab_thetaplus_neg thetaplus(i)];
                end   
            end
            prog_app = min(-1*tab_thetaplus_neg./tab_neg); 
            lambda = thetaplus + prog_app*D;
        end
    endfunction
    //*********************************************************************************************************************************
    
    //******************************************* #### ARS on the theta simplex ### ***************************************************
    function [thetaOpt, profitMax]=max_ars(theta_init)

       // initialiazing ARS parameters
       fragsigma=0.5; // Multiplicative value of the standard deviation
       imax=5; // Number of testings of different standard deviations
       napprentissage=50; // # of iterartions of each standard deviation during the learning phase
       nexploit=100; // # of iteration during the exploration phase
       maxnbcalculs=10000; // max # of computation of the objective function
       maxdersigma=5; // max number of uses of the smallest standard deviation
       nbcalculs=0;
       cycle=0;
       compteurdersigma=0;
       dersigmaActuel=%F;
       dersigmaAvant=%F;
       profitMax=-1e30;
       rendapp=-1e30;

       thetaOpt = theta_init;
       disp(thetaOpt);
       cpt_stag = 0; // a stagnation counter to assess when the profit stop being optimized at each step
       while (nbcalculs<=maxnbcalculs & compteurdersigma<=maxdersigma & cpt_stag <= 2) do
            temp_stag = 0;
            cycle=cycle+1; // We end a cycle
          // Learning phase : test of several standard deviations
            sigma=0.5; // Ecart-type
            rendapp=-1e30;
            for i=1:imax do
                thetaOpt2=thetaOpt;
                sigma=fragsigma*sigma;
                disp('Learning phase '+string(i)+' cycle '+string(cycle));
                kmax=round(napprentissage/i);
                for k=1:kmax do       
                // We choose parameters depending on tauopt and sigma :
                lambda = randomSimplexe(thetaOpt2,sigma);        
                // We evaluate the profit for lambda
                rend=Profit(lambda); 
                // If the value is better (+ greater) we keep lambda and sigma
                if rend>rendapp then
                        rendapp=rend;
                        thetaOpt2=lambda;
                        sigmaopt=sigma;
                        stock=lambda;
                        if i==imax then
                            dersigmaActuel=%T;
                        end
                end
            end
                nbcalculs=nbcalculs+kmax;
      end
            // We count the succesive uses of the smallest standard deviation	
            if dersigmaActuel  then
                if dersigmaAvant  then
                    compteurdersigma=compteurdersigma+1;
                else
                    compteurdersigma=1;
                end
            else
                compteurdersigma=0;
            end
            dersigmaAvant=dersigmaActuel;
            dersigmaActuel=%F;
            // If the learning is okay
            if rendapp>profitMax then
                profitMax=rendapp;
                thetaOpt=stock;
                temp_stag = temp_stag+1;
            end
    
            // Exploration
            disp('Exploration, cycle '+string(cycle)+' and sigmaopt =');
            disp(sigmaopt);
            disp(thetaOpt);
            for k=1:nexploit do
                // We choose parameters depending on de tauopt et sigma :
                lambda=randomSimplexe(thetaOpt,sigmaopt);
                // We evaluate the profit for lambda
                rend = Profit(lambda);  
                // If it is better (greater) we keep lambda
                if rend>profitMax then
                    profitMax=rend;
                    thetaOpt=lambda;
                    rendapp=profitMax;
                    temp_stag = temp_stag+1;
                end
            end
            nbcalculs=nbcalculs+nexploit;
            if temp_stag == 0 then
                cpt_stag = cpt_stag+1;
            end    
            //ficinter=ficparam+'dimension'+string(n)+'cycle'+string(cycle)+'.dat';
            //save(ficinter,'tauopt','cycle','profitMax','sigmaopt');
        end
        if cpt_stag >=3 then
        disp("Stagnation")
        end
        disp(string(compteurdersigma)+' successive uses of the smallest standard deviation');
        disp(string(nbcalculs)+' evaluations of the profit');
    endfunction
    [thetaOpt, profitMax]=max_ars([0.25 0.25 0.25 0.25])
    disp("The best proportions of varieties are "+string(thetaOpt)+" and it yields the profit: "+string(profitMax))
endfunction
Optim("massAction",0.0066,[100 300 600 700]) 
