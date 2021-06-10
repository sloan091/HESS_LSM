function name = unpacker(file)
tempname = load(file);
templab = fieldnames(tempname);
if length(templab)>1
    for i = 1:length(templab)
        name.(templab{i,1}) = tempname.(templab{i,1});
    end
else
name  = tempname.(templab{1,1});
end
end