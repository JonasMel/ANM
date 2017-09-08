%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                             %%%
%%% Different SBP operators used in ANM 2016    %%%
%%% ns is the number of boundary points in      %%%
%%% boundary derivative operator (in D2)        %%%
%%% alpha is how huch we can borrow from D2.    %%%
%%%                                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ordning==2)     % Andra ordningens metod
    SBP2_BV3;
    ns=3;
    alpha=2/5;
    ordningstyp=' Second order';
    
    
elseif (ordning==4)  % F?rstra ordningens upwind 3 Randpunkter
    SBP4_BV3;
    ns=4;
    alpha=.2508560249;;
    ordningstyp=' Fourth order';
    
    
elseif (ordning==6)  % F?rstra ordningens upwind 1 Randpunkt
    SBP6_BV3;
    ns=5;
    alpha=.1878715026;
    ordningstyp=' Sixth order';
    
elseif (ordning==10)  % ANNA 2016 F?rstra ordningens upwind 2 Randpunkt
    SBP10_Block
    ns=9;
    alpha=0.090909090909091;
    ordningstyp=' Tenth order block';

else
    disp('Operators does not exist');
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
