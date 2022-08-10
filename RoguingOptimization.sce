function Optim(varargin) 
    clc;
    // argument in in order : contactType,contactRate,inference,aggregation
    // contactType  = "massAction"; for a beta*U*S density-dependant infection function
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
        end
    end
    rho = 0.2
    m = 20
    mu = 1/20
    omega = 0.1882
    alpha = 0.01
    K = 10000 // Density of plants per hectare
    T = 300
    t0 = 0
    step = 0.01
    Area = 1 //hectare
    selling_price = 380000 // XOF per ton
    maintenance_price = 500000 // XOF per hectare
    // State variables are in plants/hectares and whiteflies/hectare
    
    //******************************************** Initial conditions ****************************************
             H_0 = 9500; L_0=500; I1_0=0; I2_0=0; I3_0=0; I4_0=0; I5_0=0; U0 = 6000; V0=100;


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
        function yprim= season(t,y,gammaq)
        yprim=[-beta_p*y(9)*y(1)
               beta_p*y(9)*y(1) - mu*y(2)
               mu*y(2)-alpha*y(3) - gammaq(1)*y(3)
               alpha*(y(3)-y(4)) - gammaq(2)*y(4)
               alpha*(y(4)-y(5)) - gammaq(3)*y(5)
               alpha*(y(5)-y(6)) - gammaq(4)*y(6)
               alpha*y(6) - gammaq(5)*y(7)
               rho*vertical(y(1)+y(2),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*(y(8)+y(9))*(1-(y(8)+y(9))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)))) - beta_v*y(8)*(y(3)+y(4)+y(5)+y(6)+y(7)) - omega*y(8)
               rho*vertical(y(3)+y(4)+y(5)+y(6)+y(7),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*(y(8)+y(9))*(1-(y(8)+y(9))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)))) + beta_v*y(8)*(y(3)+y(4)+y(5)+y(6)+y(7)) - omega*y(9)
              ]
         endfunction
    end
    //###### Frequency
    if varargin(1) == "frequency" then
        function yprim= season(t,y,gammaq)
        yprim=[-beta_p*y(9)*y(1)/(y(8)+y(9))
               beta_p*y(9)*y(1)/(y(8)+y(9)) - mu*y(2)
               mu*y(2)-alpha*y(3) - gammaq(1)*y(3)
               alpha*(y(3)-y(4)) - gammaq(2)*y(4)
               alpha*(y(4)-y(5)) - gammaq(3)*y(5)
               alpha*(y(5)-y(6)) - gammaq(4)*y(6)
               alpha*y(6) - gammaq(5)*y(7)
               rho*vertical(y(1)+y(2),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*(y(8)+y(9))*(1-(y(8)+y(9))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)))) - beta_v*y(8)*(y(3)+y(4)+y(5)+y(6)+y(7))/(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)) - omega*y(8)
               rho*vertical(y(3)+y(4)+y(5)+y(6)+y(7),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*(y(8)+y(9))*(1-(y(8)+y(9))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)))) + beta_v*y(8)*(y(3)+y(4)+y(5)+y(6)+y(7))/(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))  - omega*y(9)
              ]
         endfunction
    end
    
    //###### Aggregation/inference
    if varargin(1) == "inference" then
        function yprim= season(t,y,gammaq)
        yprim=[-kp*y(9)*aggregation(y(1),1-mp)
               kp*y(9)*aggregation(y(1),1-mp) - mu1*y(2)
               mu1*y(2)-alpha1*y(3) - gammaq(1)*y(3)
               alpha1*y(3)-alpha1*y(4) - gammaq(2)*y(4)
               alpha1*y(4)-alpha1*y(5) - gammaq(3)*y(5)
               alpha1*y(5)-alpha1*y(6) - gammaq(4)*y(6)
               alpha1*y(6) - gammaq(5)*y(7)
               rho*vertical(y(1)+y(2),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*(y(8)+y(9))*(1-(y(8)+y(9))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)))) - kv*(y(3)+y(4)+y(5)+y(6)+y(7))*aggregation(y(8),(1-mv)) - omega*y(8)
               rho*vertical(y(3)+y(4)+y(5)+y(6)+y(7),y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7))*(y(8)+y(9))*(1-(y(8)+y(9))/(m*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)))) + kv*(y(3)+y(4)+y(5)+y(6)+y(7))*aggregation(y(8),(1-mv)) - omega*y(9)
              ]
endfunction
    end
    
    
    //************************************************************************************************************************************
    function YY = Modele(gammaq) //Return all values at T
        t = t0+step:step:T
        Sol = ode([H_0;L_0;I1_0;I2_0;I3_0;I4_0;I5_0;U0;V0],t0,t,list(season,gammaq)); 
        nn = length(Sol(1,:));
        YY = round(Sol(:,nn))
    endfunction
    //************************************************************************************************************************************
    function yy = tuber_yield(gammaq)
        final_values = Modele(gammaq)
        ss = sum(final_values(1:7));
        index = 0;
        for i=3:length(final_values)-2
             index = index+(i-2)*final_values(i)
        end
        mean_index_of_severity = (final_values(1) + final_values(2) + index)/ss
        //disp(ss)
        //disp(mean_index_of_severity)
        yy = (-4.03*mean_index_of_severity + 31.63)*ss/K //To calibrate the yield to the number of plants/ha remaining in the field 
    end
    
    function pp = profit(gammaq)
        pp = Area*(selling_price*tuber_yield(gammaq) - maintenance_price)
    endfunction
    
    //*********************** Adaptive Random Search ***********************************
    function [gammaOpt, maxProfit]=maximizer(gammaMin,gammaMax)
    fragsigma=0.5; // Multiplicative value of the standard deviation
    imax=5; // Number of testings of different standard deviations
    napprentissage=50; // # of iterartions of each standard deviation during the learning phase
    nexploit=75; // # of iteration during the exploration phase
    maxnbcalculs=15000; // max # of computation of the objective function
    maxdersigma=5; // max number of uses of the smallest standard deviation
    nbcalculs=0;
    cycle=0;
    compteurdersigma=0;
    dersigmaActuel=%F;
    dersigmaAvant=%F;
    maxProfit=-1e30;
    coutapp=-1e30;
  
    // Optimization start
    gammaOpt=(gammaMin+gammaMax)/2;
    siz = size(gammaOpt);
  
    while (nbcalculs<=maxnbcalculs & compteurdersigma<=maxdersigma) do
        cycle=cycle+1; // A cycle ends up
        // Learning phase (we try several values of standar deviation)
        sigma=gammaMax-gammaMin; // standard deviation
        coutapp=-1e30;
        for i=1:imax do
          gammaOpt2=gammaOpt;
          sigma=fragsigma*sigma; // standard deviation
         //disp('Learning '+string(i)+' cycle '+string(cycle));
         kmax=round(napprentissage/i);
         for k=1:kmax do       
          // We choose a random point depending on thetaopt et sigma :
         lambda=gammaOpt2+grand(siz(1),siz(2),'nor',0,1).*sigma;
             while (or(lambda<gammaMin) | or(lambda>gammaMax)) do
                 // Try again if out of bounds
               lambda=gammaOpt2+grand(siz(1),siz(2),'nor',0,1).*sigma;
             end
        
                // We evaluate the profit for the randomly chosen lambda
         cout=profit(lambda); 
        // If the profit is better, we keep lambda and sigma
        if cout>coutapp then
          coutapp=cout;
          gammaOpt2=lambda;
          sigmaopt=sigma;
          stock=lambda;
          if i==imax then
            dersigmaActuel=%T;
          end
        end
      end
      nbcalculs=nbcalculs+kmax;
    end
    // We count the number of successive uses of the smallest standard deviation (it's a stopping criteria)
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
    // If the learning phase is okay
    if coutapp>maxProfit then
      maxProfit=coutapp;
      gammaOpt=stock;
    end
    
    // Exploitation
    disp('Exploration, cycle '+string(cycle)+' and sigmaopt =');
    disp(sigmaopt);
    disp(gammaOpt);
    for k=1:nexploit do
      // We choose a random point depending on thetaopt and sigma :
      lambda=gammaOpt+grand(siz(1),siz(2),'nor',0,1).*sigmaopt;
      while (or(lambda<gammaMin) | or(lambda>gammaMax)) do
        // Try again if out of bounds
        lambda=gammaOpt+grand(siz(1),siz(2),'nor',0,1).*sigmaopt;
      end
      
      // We evaluate the profit for the randomly chosen lambda
      cout = profit(lambda); 
      
      // If the profit is better, we keep lambda
      if cout>maxProfit then
        maxProfit=cout;
        gammaOpt=lambda;
        coutapp=maxProfit;
      end
    end
    nbcalculs=nbcalculs+nexploit;
    
    // Intermediate savings (no need if there's not too much computation to be done)
    //ficinter=ficparam+'cycle'+string(cycle)+'.dat';
    //save(ficinter,gammaOpt,cycle,maxProfit,sigmaopt);
  end
  
  disp(string(compteurdersigma)+' successive uses of the smaller standard deviation');
  disp(string(nbcalculs)+' computation of the profit function');

endfunction
//disp(profit([0 0 0 0 0]));abort;
[gammaOpt, maxProfit]=maximizer([0 0 0 0 0],[1 1 1 1 1])
disp("The optimal roguing strategy is:",gammaOpt)
disp("It yields the profit:", maxProfit)
endfunction
Optim("massAction",0.0006) // Implementation with density dependant transmission rate
//Optim("frequency",0.006)
