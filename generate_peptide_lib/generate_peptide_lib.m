aa = 'ACDEFGHIKLMNPQRSTVWY';
c = 1;
temp = '';
temp2 = 0;
l = 2;
while c < 20^l+1
    temp = '';
    temp2 = 0;
    for ii = 1:l
        randIdx  = randperm(length(aa),1); 
        randChar = aa(randIdx); 
        temp = strcat(temp,randChar);
    end
    if c == 1
        pep_lib(c,1) = {temp};
        c = c + 1;
    end
    if c > 1
        for kk = 1:length(pep_lib)
            f = double((temp == pep_lib{kk,1}));
            temp2 = temp2 + prod(f);
        end
        if temp2 == 0
            pep_lib(c,1) = {temp};
            c = c + 1;
        end
    end
end
writecell(pep_lib,'2mer.txt')

f = fopen('2mer.txt');
pep_lib = textscan(f,'%s');
fileID = fopen('2mer_pdb.txt','w');
for ii = 1:length(pep_lib{1,1})
    temp = char(pep_lib{1,1}(ii));
    fprintf(fileID,'fab %s-OH, chain=X, ss=4 \nattach O,1,1 \nh_add \nsave %s%s \ndelete all \n',temp,temp,'.pdb')
end
fclose(fileID);

