function new_pdb = populateTemp(pdb, col, name)
    col = col(:)';
    pdb = pdbread(pdb);
    [~,~,ix1] = unique([pdb.Model.Atom.resSeq]);
    [~,~,ix2] = unique([pdb.Model.Atom.chainID]);
    C = accumarray(ix2,1);
    g = [];
    for ii = 1:length(C)
        if min(ix1(1:length(C(ii,1)))) ~= 1
            ix1(length(g)+1:C(ii,1)) = ix1(length(g)+1:C(ii,1)) ...
                - ones(C(ii,1),1).*(min(ix1(length(g)+1:C(ii,1)))-1);
        end
        g = length(ix1(length(g)+1:length(C(ii,1))));
    end
    a = 1;
    C1 = [];
    for ii = 1:length(C)
        temp = ix1(a:C(ii,1));
        temp2 = accumarray(temp,1);
        C1((length(C1)+1):(length(C1)+length(temp2))) = temp2;
    end
    vals = repelem(col,C1)';
    for ii = 1:length(ix1)
         pdb.Model.Atom(ii).tempFactor = vals(ii);
    end
    new_pdb = pdbwrite(name,pdb);
end

