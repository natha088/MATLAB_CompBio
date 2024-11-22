function convert_to_fasta(input_file)
output_file = 'output.fasta';
sequences = readlines(input_file);
fid = fopen(output_file, 'w');
for i = 1:length(sequences)
    if strlength(sequences(i)) > 0
        fprintf(fid, '>%d\n', i); 
        fprintf(fid, '%s\n', sequences(i));
    end
end
fclose(fid);
disp('FASTA file with numerical headers has been created');
end
