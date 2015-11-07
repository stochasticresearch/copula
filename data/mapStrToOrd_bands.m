function [ ordVal ] = mapStrToOrd_bands( idx, strval )
% MAPSTRTOORD - converts a string to an ordinal number, based on which
% column the data point came from.  This script is specific to the census
% dataset.  

if(idx==2)
    % if we have a question mark 
    if(strcmpi(strval,'Private'))
        ordVal = 1;
    elseif(strcmpi(strval,'Self-emp-not-inc'))
        ordVal = 2;
    elseif(strcmpi(strval,'Self-emp-inc'))
        ordVal = 3;
    elseif(strcmpi(strval,'Federal-gov'))
        ordVal = 4;
    elseif(strcmpi(strval,'Local-gov'))
        ordVal = 5;
    elseif(strcmpi(strval,'State-gov'))
        ordVal = 6;
    elseif(strcmpi(strval,'Without-pay'))
        ordVal = 7;
    elseif(strcmpi(strval,'Never-worked'))
        ordVal = 8;
    else
        ordVal = -999;
    end
elseif(idx==4)
    if(strcmpi(strval,'Bachelors'))
        ordVal = 1;
    elseif(strcmpi(strval,'Some-college'))
        ordVal = 2;
    elseif(strcmpi(strval,'11th'))
        ordVal = 3;
    elseif(strcmpi(strval,'HS-grad'))
        ordVal = 4;
    elseif(strcmpi(strval,'Prof-school'))
        ordVal = 5;
    elseif(strcmpi(strval,'Assoc-acdm'))
        ordVal = 6;
    elseif(strcmpi(strval,'Assoc-voc'))
        ordVal = 7;
    elseif(strcmpi(strval,'9th'))
        ordVal = 8;
    elseif(strcmpi(strval,'7th-8th'))
        ordVal = 9;
    elseif(strcmpi(strval,'12th'))
        ordVal = 10;
    elseif(strcmpi(strval,'Masters'))
        ordVal = 11;
    elseif(strcmpi(strval,'1st-4th'))
        ordVal = 12;
    elseif(strcmpi(strval,'10th'))
        ordVal = 13;
    elseif(strcmpi(strval,'Doctorate'))
        ordVal = 14;
    elseif(strcmpi(strval,'5th-6th'))
        ordVal = 15;
    elseif(strcmpi(strval,'Preschool'))
        ordVal = 16;
    else
        ordVal = -999;
    end
elseif(idx==6)
    if(strcmpi(strval,'Married-civ-spouse'))
        ordVal = 1;
    elseif(strcmpi(strval,'Divorced'))
        ordVal = 2;
    elseif(strcmpi(strval,'Never-married'))
        ordVal = 3;
    elseif(strcmpi(strval,'Separated'))
        ordVal = 4;
    elseif(strcmpi(strval,'Widowed'))
        ordVal = 5;
    elseif(strcmpi(strval,'Married-spouse-absent'))
        ordVal = 6;
    elseif(strcmpi(strval,'Married-AF-spouse'))
        ordVal = 7;
    else
        ordVal = -999;
    end
elseif(idx==7)
    if(strcmpi(strval,'Tech-support'))
        ordVal = 1;
    elseif(strcmpi(strval,'Craft-repair'))
        ordVal = 2;      
    elseif(strcmpi(strval,'Other-service'))
        ordVal = 3;
    elseif(strcmpi(strval,'Sales'))
        ordVal = 4;
    elseif(strcmpi(strval,'Exec-managerial'))
        ordVal = 5;
    elseif(strcmpi(strval,'Prof-specialty'))
        ordVal = 6;
    elseif(strcmpi(strval,'Handlers-cleaners'))
        ordVal = 7;
    elseif(strcmpi(strval,'Machine-op-inspct'))
        ordVal = 8;
    elseif(strcmpi(strval,'Adm-clerical'))
        ordVal = 9;
    elseif(strcmpi(strval,'Farming-fishing'))
        ordVal = 10;
    elseif(strcmpi(strval,'Transport-moving'))
        ordVal = 11;
    elseif(strcmpi(strval,'Priv-house-serv'))
        ordVal = 12;
    elseif(strcmpi(strval,'Protective-serv'))
        ordVal = 13;
    elseif(strcmpi(strval,'Armed-Forces'))
        ordVal = 14;
    else
        ordVal = -999;
    end
elseif(idx==8)
    % relationship: Wife, Own-child, Husband, Not-in-family, Other-relative, Unmarried.
        if(strcmpi(strval,'Wife'))
        ordVal = 1;
    elseif(strcmpi(strval,'Own-child'))
        ordVal = 2;      
    elseif(strcmpi(strval,'Husband'))
        ordVal = 3;
    elseif(strcmpi(strval,'Not-in-family'))
        ordVal = 4;
    elseif(strcmpi(strval,'Other-relative'))
        ordVal = 5;
    elseif(strcmpi(strval,'Unmarried'))
        ordVal = 6;
    else
        ordVal = -999;
    end
elseif(idx==9)
    % race: White, Asian-Pac-Islander, Amer-Indian-Eskimo, Other, Black.
    if(strcmpi(strval,'White'))
        ordVal = 1;
    elseif(strcmpi(strval,'Asian-Pac-Islander'))
        ordVal = 2;      
    elseif(strcmpi(strval,'Amer-Indian-Eskimo'))
        ordVal = 3;
    elseif(strcmpi(strval,'Other'))
        ordVal = 4;
    elseif(strcmpi(strval,'Black'))
        ordVal = 5;
    else
        ordVal = -999;
    end
elseif(idx==10)
    % sex: Female, Male.
    if(strcmpi(strval,'Female'))
        ordVal = 1;
    elseif(strcmpi(strval,'Male'))
        ordVal = 2;      
    else
        ordVal = -999;
    end
elseif(idx==14)
    if(strcmpi(strval,'United-States'))
        ordVal = 1;
    elseif(strcmpi(strval,'Cambodia'))
        ordVal = 2;
    elseif(strcmpi(strval,'England'))
        ordVal = 3;
    elseif(strcmpi(strval,'Puerto-Rico'))
        ordVal = 4;
    elseif(strcmpi(strval,'Canada'))
        ordVal = 5;
    elseif(strcmpi(strval,'Germany'))
        ordVal = 6;
    elseif(strcmpi(strval,'Outlying-US(Guam-USVI-etc)'))
        ordVal = 7;
    elseif(strcmpi(strval,'India'))
        ordVal = 8;
    elseif(strcmpi(strval,'Japan'))
        ordVal = 9;
    elseif(strcmpi(strval,'Greece'))
        ordVal = 10;
    elseif(strcmpi(strval,'South'))
        ordVal = 11;
    elseif(strcmpi(strval,'China'))
        ordVal = 12;
    elseif(strcmpi(strval,'Cuba'))
        ordVal = 13;
    elseif(strcmpi(strval,'Iran'))
        ordVal = 14;
    elseif(strcmpi(strval,'Honduras'))
        ordVal = 15;
    elseif(strcmpi(strval,'Philippines'))
        ordVal = 16;
    elseif(strcmpi(strval,'Italy'))
        ordVal = 17;
    elseif(strcmpi(strval,'Poland'))
        ordVal = 18;
    elseif(strcmpi(strval,'Jamaica'))
        ordVal = 19;
    elseif(strcmpi(strval,'Vietnam'))
        ordVal = 20;
    elseif(strcmpi(strval,'Mexico'))
        ordVal = 21;
    elseif(strcmpi(strval,'Portugal'))
        ordVal = 22;
    elseif(strcmpi(strval,'Ireland'))
        ordVal = 23;
    elseif(strcmpi(strval,'France'))
        ordVal = 24;
    elseif(strcmpi(strval,'Dominican-Republic'))
        ordVal = 25;
    elseif(strcmpi(strval,'Laos'))
        ordVal = 26;
    elseif(strcmpi(strval,'Ecuador'))
        ordVal = 27;
    elseif(strcmpi(strval,'Taiwan'))
        ordVal = 28;
    elseif(strcmpi(strval,'Haiti'))
        ordVal = 29;
    elseif(strcmpi(strval,'Columbia'))
        ordVal = 30;
    elseif(strcmpi(strval,'Hungary'))
        ordVal = 31;
    elseif(strcmpi(strval,'Guatemala'))
        ordVal = 32;
    elseif(strcmpi(strval,'Nicaragua'))
        ordVal = 33;
    elseif(strcmpi(strval,'Scotland'))
        ordVal = 34;
    elseif(strcmpi(strval,'Thailand'))
        ordVal = 35;
    elseif(strcmpi(strval,'Yugoslavia'))
        ordVal = 36;
    elseif(strcmpi(strval,'El-Salvador'))
        ordVal = 37;
    elseif(strcmpi(strval,'Trinadad&Tobago'))
        ordVal = 38;
    elseif(strcmpi(strval,'Peru'))
        ordVal = 39;
    elseif(strcmpi(strval,'Hong'))
        ordVal = 40;
    elseif(strcmpi(strval,'Holand-Netherlands'))
        ordVal = 41;
    else
        ordVal = -999;
    end
elseif(idx==15)
    if(strcmpi(strval,'<=50K'))
        ordVal = 1;
    elseif(strcmpi(strval,'>50k'))
        ordVal = 2;
    else
        ordVal = -999;
    end
else
    warning('Invalid IDX specified!');
    ordVal = -999;
end


end

