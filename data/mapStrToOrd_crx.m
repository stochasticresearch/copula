function [ ordVal ] = mapStrToOrd_crx( idx, strval )
% MAPSTRTOORD - converts a string to an ordinal number, based on which
% column the data point came from.  This script is specific to the census
% dataset.  

if(idx==1)
    % if we have a question mark 
    if(strcmpi(strval,'b'))
        ordVal = 1;
    elseif(strcmpi(strval,'a'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==4)
    if(strcmpi(strval,'u'))
        ordVal = 1;
    elseif(strcmpi(strval,'y'))
        ordVal = 2;
    elseif(strcmpi(strval,'l'))
        ordVal = 3;
    elseif(strcmpi(strval,'t'))
        ordVal = 4;
    else
        ordVal = -999;
    end
elseif(idx==5)
    if(strcmpi(strval,'g'))
        ordVal = 1;
    elseif(strcmpi(strval,'p'))
        ordVal = 2;
    elseif(strcmpi(strval,'gg'))
        ordVal = 3;
    else
        ordVal = -999;
    end
elseif(idx==6)
    if(strcmpi(strval,'c'))
        ordVal = 1;
    elseif(strcmpi(strval,'d'))
        ordVal = 2;      
    elseif(strcmpi(strval,'cc'))
        ordVal = 3;
    elseif(strcmpi(strval,'i'))
        ordVal = 4;
    elseif(strcmpi(strval,'j'))
        ordVal = 5;
    elseif(strcmpi(strval,'k'))
        ordVal = 6;
    elseif(strcmpi(strval,'m'))
        ordVal = 7;
    elseif(strcmpi(strval,'r'))
        ordVal = 8;
    elseif(strcmpi(strval,'q'))
        ordVal = 9;
    elseif(strcmpi(strval,'w'))
        ordVal = 10;
    elseif(strcmpi(strval,'x'))
        ordVal = 11;
    elseif(strcmpi(strval,'e'))
        ordVal = 12;
    elseif(strcmpi(strval,'aa'))
        ordVal = 13;
    elseif(strcmpi(strval,'ff'))
        ordVal = 14;
    else
        ordVal = -999;
    end
elseif(idx==7)
    if(strcmpi(strval,'v'))
        ordVal = 1;
    elseif(strcmpi(strval,'h'))
        ordVal = 2;      
    elseif(strcmpi(strval,'bb'))
        ordVal = 3;
    elseif(strcmpi(strval,'j'))
        ordVal = 4;
    elseif(strcmpi(strval,'n'))
        ordVal = 5;
    elseif(strcmpi(strval,'z'))
        ordVal = 6;
    elseif(strcmpi(strval,'dd'))
        ordVal = 7;
    elseif(strcmpi(strval,'ff'))
        ordVal = 8;
    elseif(strcmpi(strval,'o'))
        ordVal = 9;
    else
        ordVal = -999;
    end
elseif(idx==9)
    if(strcmpi(strval,'t'))
        ordVal = 1;
    elseif(strcmpi(strval,'f'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==10)
    if(strcmpi(strval,'t'))
        ordVal = 1;
    elseif(strcmpi(strval,'f'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==12)
    if(strcmpi(strval,'t'))
        ordVal = 1;
    elseif(strcmpi(strval,'f'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==13)
    if(strcmpi(strval,'g'))
        ordVal = 1;
    elseif(strcmpi(strval,'p'))
        ordVal = 2;
    elseif(strcmpi(strval,'s'))
        ordVal = 3;
    else
        ordVal = -999;
    end
elseif(idx==16)
    if(strcmpi(strval,'+'))
        ordVal = 1;
    elseif(strcmpi(strval,'-'))
        ordVal = 2;
    else
        ordVal = -999;
    end
else
    warning('Invalid IDX specified!');
    ordVal = -999;
end

end

