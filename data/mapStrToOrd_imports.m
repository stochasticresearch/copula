function [ ordVal ] = mapStrToOrd_imports( idx, strval )
% MAPSTRTOORD - converts a string to an ordinal number, based on which
% column the data point came from.  This script is specific to the census
% dataset.  

if(idx==1)
    % if we have a question mark 
    if(strcmpi(strval,'-3'))
        ordVal = 1;
    elseif(strcmpi(strval,'-2'))
        ordVal = 2;
    elseif(strcmpi(strval,'-1'))
        ordVal = 3;
    elseif(strcmpi(strval,'-0'))
        ordVal = 4;
    elseif(strcmpi(strval,'-1'))
        ordVal = 5;
    elseif(strcmpi(strval,'-2'))
        ordVal = 6;
    elseif(strcmpi(strval,'-3'))
        ordVal = 7;
    else
        ordVal = -999;
    end
elseif(idx==3)
    if(strcmpi(strval,'alfa-romero'))
        ordVal = 1;
    elseif(strcmpi(strval,'audi'))
        ordVal = 2;
    elseif(strcmpi(strval,'bmw'))
        ordVal = 3;
    elseif(strcmpi(strval,'chevrolet'))
        ordVal = 4;
    elseif(strcmpi(strval,'dodge'))
        ordVal = 5;
    elseif(strcmpi(strval,'honda'))
        ordVal = 6;
    elseif(strcmpi(strval,'isuzu'))
        ordVal = 7;
    elseif(strcmpi(strval,'jaguar'))
        ordVal = 8;
    elseif(strcmpi(strval,'mazda'))
        ordVal = 9;
    elseif(strcmpi(strval,'mercedes-benz'))
        ordVal = 10;
    elseif(strcmpi(strval,'mercury'))
        ordVal = 11;
    elseif(strcmpi(strval,'mitsubishi'))
        ordVal = 12;
    elseif(strcmpi(strval,'nissan'))
        ordVal = 13;
    elseif(strcmpi(strval,'peugot'))
        ordVal = 14;
    elseif(strcmpi(strval,'plymouth'))
        ordVal = 15;
    elseif(strcmpi(strval,'porsche'))
        ordVal = 16;
    elseif(strcmpi(strval,'renault'))
        ordVal = 17;
    elseif(strcmpi(strval,'saab'))
        ordVal = 18;
    elseif(strcmpi(strval,'subaru'))
        ordVal = 19;
    elseif(strcmpi(strval,'toyato'))
        ordVal = 20;
    elseif(strcmpi(strval,'volkswagen'))
        ordVal = 21;
    elseif(strcmpi(strval,'volvo'))
        ordVal = 22;        
    else
        ordVal = -999;
    end
elseif(idx==4)
    if(strcmpi(strval,'diesel'))
        ordVal = 1;
    elseif(strcmpi(strval,'gas'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==5)
    if(strcmpi(strval,'std'))
        ordVal = 1;
    elseif(strcmpi(strval,'turbo'))
        ordVal = 2;      
    else
        ordVal = -999;
    end
elseif(idx==6)
    if(strcmpi(strval,'four'))
        ordVal = 1;
    elseif(strcmpi(strval,'two'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==7)
    if(strcmpi(strval,'hardtop'))
        ordVal = 1;
    elseif(strcmpi(strval,'wagon'))
        ordVal = 2;
    elseif(strcmpi(strval,'sedan'))
        ordVal = 3;
    elseif(strcmpi(strval,'hatchback'))
        ordVal = 4;
    elseif(strcmpi(strval,'convertible'))
        ordVal = 5;
    else
        ordVal = -999;
    end
elseif(idx==8)
    if(strcmpi(strval,'4wd'))
        ordVal = 1;
    elseif(strcmpi(strval,'fwd'))
        ordVal = 2;
    elseif(strcmpi(strval,'rwd'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==9)
    if(strcmpi(strval,'front'))
        ordVal = 1;
    elseif(strcmpi(strval,'rear'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==15)
    if(strcmpi(strval,'dohc'))
        ordVal = 1;
    elseif(strcmpi(strval,'dohcv'))
        ordVal = 2;
    elseif(strcmpi(strval,'l'))
        ordVal = 3;
    elseif(strcmpi(strval,'ohc'))
        ordVal = 3;
    elseif(strcmpi(strval,'ohcf'))
        ordVal = 3;
    elseif(strcmpi(strval,'ohcv'))
        ordVal = 3;
    elseif(strcmpi(strval,'rotor'))
        ordVal = 3;
    else
        ordVal = -999;
    end
elseif(idx==16)
    if(strcmpi(strval,'eight'))
        ordVal = 1;
    elseif(strcmpi(strval,'five'))
        ordVal = 2;
    elseif(strcmpi(strval,'four'))
        ordVal = 2;
    elseif(strcmpi(strval,'six'))
        ordVal = 2;
    elseif(strcmpi(strval,'three'))
        ordVal = 2;
    elseif(strcmpi(strval,'twelve'))
        ordVal = 2;
    elseif(strcmpi(strval,'two'))
        ordVal = 2;
    else
        ordVal = -999;
    end
elseif(idx==18)
    if(strcmpi(strval,'1bbl'))
        ordVal = 1;
    elseif(strcmpi(strval,'2bbl'))
        ordVal = 2;
    elseif(strcmpi(strval,'4bbl'))
        ordVal = 3;
    elseif(strcmpi(strval,'idi'))
        ordVal = 4;
    elseif(strcmpi(strval,'mfi'))
        ordVal = 5;
    elseif(strcmpi(strval,'mpfi'))
        ordVal = 6;
    elseif(strcmpi(strval,'spdi'))
        ordVal = 7;
    elseif(strcmpi(strval,'spfi'))
        ordVal = 8;
    else
        ordVal = -999;
    end
else
    warning('Invalid IDX specified!');
    ordVal = -999;
end

end

