//
// Created by sheepless on 3/26/19.
//
#include "get_relivent_snps.h"


void cleaning_func(int arr1[MAXARR], int arr2[MAXARR], int insamp, bool * clean_arr1, bool * clean_arr2, int * out_n){
    int j = 0;
    for (int i = 0; i < insamp ; ++i) {
        if(arr1[i] == 0 || arr2[i] == 0) {
            *out_n = (*out_n) - 1;
            continue ;
        } else {
            clean_arr1[j] = ( arr1[i] == 1 );
            clean_arr2[j] = ( arr2[i] == 1 );
            j++;
        }
    }
    * out_n = j;
}

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
    std::cout <<" \t" << "+"  << " \t " << "-"  << std::endl
         <<"+\t" << n_AB << " \t " << n_Ab << std::endl
         <<"-\t" << n_aB << " \t " << n_ab << std::endl;

}


float gammitic_disequalibrium(bool gene_A[MAXARR],  bool gene_B[MAXARR], int samples ){
    float   n_AB = 0, n_A = 0, n_B = 0 , n_all = (float) samples; // where A = is ref
    float   p_AB = 0, p_A = 0, p_B = 0;
    // calculate nums of each
    for (int i = 0; i < samples; ++i) {
        if(gene_A[i]){
            n_A++;
            if(gene_B[i]){
                n_B++;
                n_AB++;
            }
        } else if(gene_B[i]){
            n_B++;
        }
    }
    p_A  = (n_A  / n_all);
    p_B  = (n_B  / n_all);
    p_AB = (n_AB / n_all);

    return (p_AB - (p_A * p_B));
}

int permutation_test(bool A[MAXARR], bool B[MAXARR], int max , double  D_stat,  std::string * path_to_log){
    int i ;
    double p_num = 0.0 ,p_value = 0.0;
    srand(time_t(NULL));
    for (i = 0; i < RESAMPLES ; ++i) { // do the resampling for the permutation
        std::random_shuffle(&B[0],&B[max]);
        if (abs(gammitic_disequalibrium(A,B,max)) > abs(D_stat)){
            p_num++;
        }
        p_value = p_num / (double) i;
    }
    std::cout << p_value <<  "\t" << D_stat << "\t" << max << std::endl;
    return 0;
}

int boot_strap(bool A[MAXARR], bool B[MAXARR], int max, float alpha, std::string * path_to_log ){
    float resampled_disequalibriums[BOOT] = {};
    int index ;
    bool A_loop[MAXARR] = {};
    bool B_loop[MAXARR] = {};
    std::cerr << "declared everything\n";
    for (int i = 0; i < BOOT ; ++i) {
        for (int j = 0; j < max; ++j) {     // fill the loop arrays with randome tuples from A and B
            index = ((int) random() % max);
            A_loop[j] = A[index];
            B_loop[j] = B[index];
        }
        resampled_disequalibriums[i] = gammitic_disequalibrium(A_loop,B_loop,max); // fill the alt distrabution
    }
    std::sort(&resampled_disequalibriums[0],&resampled_disequalibriums[BOOT]);
    std::cout << "Lower-Bound= " << resampled_disequalibriums[(int)((alpha/2)*BOOT)] << "\t"
              << "Upper-Bound= " << resampled_disequalibriums[(int)((1-(alpha/2))*BOOT)] << std::endl;
    return 0;
}

int get_stats(int a_raw[MAXARR], int b_raw[MAXARR], int init_num , float * alpha,std::string path_to_log){
    bool A[MAXARR] = {};
    bool B[MAXARR] = {};
    int num = 0;
    bool * A_p = &A[0];
    bool * B_p = &B[0];
    int  * num_p = &num ;
    int i ;
    double test_stat = 0.0 , p_num = 0.0 ,p_value = 0.0;

    cleaning_func(a_raw,b_raw,init_num, A_p, B_p , num_p);
    test_stat = gammitic_disequalibrium(A,B,num);
    print_table(A,B,num); // prints grapphical table of data
    srand(time_t(NULL));
    boot_strap(A, B, num, 0.05 , &path_to_log); // do bootstrap
    permutation_test(A, B, num, test_stat, &path_to_log); // do permutation test
    return 0;
}

int main(){
    // test values
    float alpha = 0.05;
    int a_raw[256] = {1,1,1,1,1,1,1,1,1,1,1,1,0,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0};
    int b_raw[256] = {1,1,1,1,1,2,2,2,2,2,2,0,1,1,1,1,1,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0};
    int init_num = 20;
    if(get_stats(a_raw,b_raw,init_num,&alpha,"") != 0 ) {
        std::cerr << "Permutation test failed\n";
    }
    return 0;
}

