//
// Created by sheepless on 5/25/19.
//


#include "IO.h"

void print_table(bool gene_A[MAXARR], bool gene_B[MAXARR], int samples){
    int n_AB = 0, n_Ab = 0, n_aB = 0, n_ab = 0;
    for (int i = 0; i < samples; ++i) {
        if(gene_A[i]){
            if(gene_B[i]){
                n_AB++;
            } else {
                n_Ab++;
            }
        } else {
            if(gene_B[i]){
                n_aB++;
            } else{
                n_ab++;
            }
        }
    }
    std::cerr <<" \t" << "+"  << " \t " << "-"  << std::endl
              <<"+\t" << n_AB << " \t " << n_Ab << std::endl
              <<"-\t" << n_aB << " \t " << n_ab << std::endl;
}

int write_vals(char* samp[MAXARR], int positions[MAXARR], int data[MAXARR][MAXARR], int fill[2]){
    std::ofstream file;
    file.open(ARGS.outpath + "_raw.tsv");
    int target_index = 0;
    file << "Samples";
    for (int k = 0; k < fill[1]; ++k) {
        file << "\t" << positions[k];
        if(positions[k] == ARGS.target){
            target_index = k;
        }
    }
    file << std::endl;
    for (int i = 0; i < fill[0]; ++i) {
        file << samp[i];
        for (int j = 0; j <fill[1] ; ++j) {
            file << "\t" <<data[j][i];
        }
        file << std::endl;
    }
    file.close();
    return target_index;
}
int write_values(out_data data[MAXARR], int max, std::string & out_path){
    std::ofstream pairwise_data;
    pairwise_data.open( out_path + ".summary.tsv");
    pairwise_data << "Pairing\t Lower_Bound\t Upper_Bound\t "
                  << "Sample_Disequilibrium\t p-value\t "
                  << "alpa-adj\t Decision";
    for (int i = 0; i < max ; ++i) { // loop through all values and print them
        pairwise_data << std::endl
                      << data[i].names     << "\t"
                      << data[i].LB        << "\t"
                      << data[i].UB        << "\t"
                      << data[i].D_stat    << "\t"
                      << data[i].p_value   << "\t"
                      << data[i].adj_alpha << "\t"
                      << data[i].reject;
    }
    pairwise_data.close();
    return 0;
}
