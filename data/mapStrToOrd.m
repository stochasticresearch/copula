function [ ordVal ] = mapStrToOrd( idx, strval, dataset )

if(strcmpi(dataset,'adult'))
    ordVal = mapStrToOrd_adult(idx,strval);
elseif(strcmpi(dataset,'bands'))
    ordVal = mapStrToOrd_bands(idx,strval);
elseif(strcmpi(dataset,'crx'))
    ordVal = mapStrToOrd_crx(idx,strval);
elseif(strcmpi(dataset,'imports'))
    ordVal = mapStrToOrd_imports(idx,strval);
end

end

