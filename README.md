# MATLAB_CompBio
MATLAB Scripts for Various Computational Biology Tasks
## populateTemp
This function can swap out the tempFactor column in a PDB file with any numeric column vector inputted by the user. This function can be used to populate per residue continuous metrics from molecular dynamics such as RMSD, surface area, hydrogen bond lifetimes; or per residue discrete/categorical metrics such as sites of post-translational modification or disease-associated mutations.

populateTemp(pdb, col, name) 
pdb: filename
col: numeric column vector (M x N, where N is number of residues in PDB structure and M is the number of chains).
name: new filename

Usage Example:

populateTemp('6b1g.pdb',[1:80;1:80],'new_6b1g.pdb')


