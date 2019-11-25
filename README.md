# linkage_calc
Source files to create linkage calc, a program that determines if haploid samples in a multi-sample vcf have alleles that are linked, by using a gamitic linkage calculation then permitation tests. Developed initioly for use with the NOTCH2NL genes, but is fully generalized.  

Installation:
make sure you have the required dependincies for htslib (pthread lzip)

make sure you have cmake 3.10.0 or higher installed 

run the command "cmake CMakeLists.txt && make "

then move the binary "linkage" into your path.

