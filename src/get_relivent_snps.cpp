//
// Created by sheepless on 3/19/19.
//
#include "get_relivent_snps.h"
#if DBUG_V
void print_vec(bool v[], int len){
    for(int i = 0 ; i < len ; i++){
        std::cerr << v[i] << " ";
    }
    std::cerr << std::endl;
}
void print_vec(int v[], int len){
    for(int i = 0 ; i < len ; i++){
        std::cerr << v[i] << " ";
    }
    std::cerr << std::endl;
}
#endif
int  get_int_type_fmt(bcf1_t * data, std::string& tag , bcf_idpair_t *pair_of_vals ,int nsamp, int * ints){
    int j = 0, tag_id ;
#if DBUG_V
    if(ARGS.debug_lvl > 1){std::cerr << "nsamp has a value of : \t" << nsamp << std::endl;}
#endif
    if (nsamp > 256) {
        std::cerr << "Can't parse over 256 samples!\n";
        return  -1; }
    // figure out the id of the tag the user input
    while (true){
#if DBUG_V
        if(ARGS.debug_lvl > 1){std::cerr << (std::string)pair_of_vals[j].key << "\t " << tag<< std::endl;}
#endif
        if((std::string) pair_of_vals[data->d.fmt[j].id].key ==  tag){
#if DBUG_V
            if(ARGS.debug_lvl > 1){std::cerr << (std::string)pair_of_vals[data->d.fmt[j].id].key << "\t" << j << std::endl;}
#endif
            tag_id = j;
            break;
        }
        j++;
    }
    // use the id to get the list of int values for each sample (corresponding to the tag) at that pos.
    for (int i = 0; i < nsamp; i++) {
#if DBUG_V
        if(ARGS.debug_lvl > 1){std::cerr << ((unsigned int)(data->d.fmt[tag_id].p[i])) % 128  << "\t"; }
#endif
        ints[i] = ((unsigned int)(data->d.fmt[tag_id].p[i])) % 128 ;
    }
    return 0;
}

int Bonferroni(int m, float alpha){
    float cor_alpha = alpha/(m-1); // alpha under Bonferroni SFER
    for (int i = 0; i < m; ++i) {
        if (DATA[i].names != "self"){
            DATA[i].reject = (DATA[i].p_value < cor_alpha);
        }
    }
    return 0;
}

int Holm_Bonferroni(int m , float alpha){
    std::vector <out_data> temp_hold;
    std::vector < std::pair <float , int > > Pval_Index;
    temp_hold.resize(m);
    Pval_Index.resize(m);
    int k = 0 , true_m, rt_val = 0;
    for(int i = 0; i < m; i++){
        Pval_Index[i] = (std::make_pair(DATA[i].p_value,i));
        temp_hold[i] = DATA[i];
    }
    std::sort(Pval_Index.begin(),Pval_Index.end());
    if((temp_hold[Pval_Index[0].second].names !="self")){
        std::cerr << "Somthing funky happened, the lowest p value is not the comparison ageist itself...\n"
                  << "going to include all \"comparisons\" including potentially self\n";
        rt_val = -1;
    } else {
        k = 1; //skip index 0
    }
    for (int j = 0; j < m ; ++j) {
        DATA[j] = temp_hold[Pval_Index[j].second];
        if (k == 0){
            DATA[j].reject = (DATA[j].p_value < (alpha)/(m + 1 - j));
        }else if(k == 1 && j > 0){
            DATA[j].reject = (DATA[j].p_value < (alpha)/(m - j));
        }
    }
    return rt_val;
}

int do_correction(int correction_type,int num_tests,float alpha){
    if(correction_type == 0){ return Bonferroni(num_tests,alpha);}
    else if (correction_type == 1) { return Holm_Bonferroni(num_tests,alpha);}
    else { return  -9;}
}

int mk_grid(htsFile * bcf, bcf_hdr_t * hdr, int * poss, int * arr, char ** header[MAXARR], int * dem[2]){
    // initialise the data
    int i = 0 , j = 0 , k = 0 , l = 0;
    int cur_conting_pos = 0 , num_not_null = 0;
    int NumSamples = 1 ;
    int var_type[MAXARR] = {} , ref_type[MAXARR] = {}; // make some space for fmt tag values
    int * var_p = &var_type[0]; // pointer to the first pos of the var_type array
    int * ref_p = &ref_type[0]; // pointer to the first pos of the ref_type array
    int end = ARGS.StartptEndpt[1];
    int start = ARGS.StartptEndpt[0];
    bool read_some = false, snp = false;
    unsigned short int store = 0;

    // make some space in the stack
    bcf1_t * line = bcf_init();

    // Get the number of samples in the multisample VCF.
    NumSamples = bcf_hdr_nsamples(hdr);
    if(ARGS.debug_lvl > 0){std::cout << "Number of samples found:\t" << NumSamples << std::endl;}
    while (i < NumSamples && i < MAXARR){   // make sure we dont over flow anything but loop through all
        header[i] = &hdr->samples[i]; // rip the sample names from the header
        if(ARGS.debug_lvl > 0){std::cerr << *(header[i]) << "\t" << i << std::endl;}
        i++;
    }
    // see if we filled the samples or the numver of strings we have
    if ( i == 255){
        std::cerr << "Filled the header line, were there more then 256?... If so split up the samples\n";
        NumSamples = 256;
    }

    //loop through all positions in the VCF
    while(bcf_read(bcf,hdr,line) == 0){
        if(bcf_hdr_id2name(hdr,line->rid) == ARGS.contig) { //see if it is the right contig
            cur_conting_pos = line->pos + 1 ; // make things more readable, and convert to 1 indexed
            if (cur_conting_pos < start) { // if it is before the start... continue
                continue;
            } else if (cur_conting_pos <= end) {
                if (ARGS.debug_lvl > 0) { std::cerr << cur_conting_pos << std::endl; } // testing
                if (bcf_unpack(line,BCF_UN_IND) != 0) { return -2;}// Unpack the data for the individual samples on the lines we care about

                if (get_int_type_fmt(line,ARGS.var_fmt_flag,hdr->id[BCF_DT_ID],NumSamples, var_p ) != 0 ){ return -3;};
                if (get_int_type_fmt(line,ARGS.ref_fmt_flag,hdr->id[BCF_DT_ID],NumSamples, ref_p ) != 0 ){ return -4;};
                num_not_null = 0; // starting a new line, so reset this;
                snp = false;
                for (int m = 0; m < NumSamples; m++) {
                    if(var_type[m] < ref_type[m]){
                        store = 1;
                        num_not_null++;
                    } // set to 1 b/c 0 is null, and this means it is more reference
                    else if (var_type[m] == 0 && ref_type[m] == 0){store = 0;} // they are both 0 so... we really have no data, set it to null
                    else {
                        store = 2;
                        num_not_null++;
                        snp = true;
                    } // we ruled all other cases out, so this is a variant here.
#if DBUG_V
                    if(ARGS.debug_lvl > 1){std::cout<< var_type[m] <<  "," << ref_type[m] << ":";}
#endif
                    if(ARGS.debug_lvl > 0){std::cerr << store << "\t";}
                    arr[m + (MAXARR*j)] = store; // add value to the int array
                }
                // we only want lines with at least one snp
                if((num_not_null >= ARGS.min_not_null && snp) or (cur_conting_pos == ARGS.target )) {
#if DBUG_V
                    if(ARGS.debug_lvl > 1){std::cerr << "adding position " << cur_conting_pos << " and starting next row!\n";}
#endif
                    poss[j] = cur_conting_pos;
                    j++;
                }
                if(DBUG){std::cout << "\n";}

                read_some = true;
            } else { // if it is after the endpoint, break the loop
                break;
            }
            if (j == MAXARR - 1) { // if the array if full break, and tell user
                std::cerr << "256 variants found, breaking!\n";
                break;
            }
        } else if (read_some){
            std::cerr << "Hit the next contig before user endpoint...\n Breaking here to not add incorrect values\n" << std::endl;
            break;
        }
    }
    * dem[0] = NumSamples;
    * dem[1] = j;
    return 0;
}


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

int permutation_test(bool A[MAXARR], bool B[MAXARR], int max , std::string path_to_log , out_data * report){
    int i ;
    double p_num = 0.0;
    float D_loop = 0.0;
    std::ofstream perm_log;
    bool logging = false;
    if(path_to_log != ""){
        perm_log.open( path_to_log + ".permutation.csv");
        logging = true;
    }
    for (i = 0; i < RESAMPLES ; ++i) { // do the resampling for the permutation
        std::random_shuffle(&B[0],&B[max]);
        D_loop = gammitic_disequalibrium(A,B,max);
        if (abs(D_loop) >= abs(report->D_stat)){
            p_num++;
        }
        if(logging){
            perm_log << D_loop << std::endl;
        }
    }
    report->p_value = (float) (p_num / (double) RESAMPLES);
    if((ARGS.debug_lvl > 0) || !(logging)) {
        std::cerr << report->p_value << "\t" << report->D_stat << "\t" << max << std::endl;
    }
    if(logging){
        perm_log.close();
    }
    return 0;
}

int boot_strap(bool A[MAXARR], bool B[MAXARR], int max, float alpha, std::string path_to_log ,out_data * report){
    auto resampled_disequalibriums = (float *) malloc(RESAMPLES * sizeof(float));
    int index ;
    bool A_loop[MAXARR] = {};
    bool B_loop[MAXARR] = {};
    bool logging = false;
    float lower = 0.0 , upper = 0.0;
    std::ofstream boot_log ;
    if (path_to_log != ""){ // see if we want to log the alt distrabution
        boot_log.open(path_to_log + ".bootstrap.csv");
        logging = true;
    }
    for (int i = 0; i < RESAMPLES ; ++i) {
        for (int j = 0; j < max; ++j) {     // fill the loop arrays with randome tuples from A and B
            index = ((int) random() % max);
            A_loop[j] = A[index];
            B_loop[j] = B[index];
        }
        resampled_disequalibriums[i] = gammitic_disequalibrium(A_loop,B_loop,max); // fill the alt distrabution
        if(logging){
            boot_log << resampled_disequalibriums[i] << std::endl; // add the values
        }
    }
    std::sort(&resampled_disequalibriums[0],&resampled_disequalibriums[RESAMPLES]);
    lower = resampled_disequalibriums[(int) ((alpha / 2) * RESAMPLES)];
    upper = resampled_disequalibriums[(int) ((1 - (alpha / 2)) * RESAMPLES)];
    if( (ARGS.debug_lvl > 0)  || !(logging)) {     // see if we should print to screen if so, do
        std::cerr << "Lower-Bound= " << lower << "\t"
                  << "Upper-Bound= " << upper << std::endl;
    }
    if(logging){
        boot_log.close();
    }
    report->LB = lower;
    report->UB = upper;
    return 0;
}

int  get_stats(int a_raw[MAXARR], int b_raw[MAXARR], int init_num, float * alpha, std::string path_to_log, int index_to_fill){
    bool A[MAXARR] = {};
    bool B[MAXARR] = {};
    int num = 0;
    bool * A_p = &A[0];
    bool * B_p = &B[0];
    int  * num_p = &num ;
#if DBUG_V
    if(ARGS.debug_lvl > 1){
       print_vec(a_raw, init_num);
       print_vec(b_raw, init_num);
    }
#endif
    cleaning_func(a_raw,b_raw,init_num, A_p, B_p , num_p);
#if DBUG_V
    if(ARGS.debug_lvl > 1){
       print_vec(A_p, num);
       print_vec(B_p, num);
    }
#endif
    DATA[index_to_fill].D_stat = gammitic_disequalibrium(A,B,num);
    print_table(A,B,num); // prints grapphical table of data
    srand(time_t(NULL));
    boot_strap(A, B, num, * alpha , path_to_log  + DATA[index_to_fill].names, &DATA[index_to_fill]); // do bootstrap
    permutation_test(A, B, num , path_to_log  + DATA[index_to_fill].names, &DATA[index_to_fill]); // do permutation test
    DATA[index_to_fill].reject = (DATA[index_to_fill].p_value < *alpha);
    DATA[index_to_fill].state = "done";
    return 0;
}

void *child_ps(void* in_vals) {     // (int Alpha[MAXARR], int Beta[MAXARR], int num , user_arguments & args, out_data * return_vals){
    auto vals_cast = (PS_Data * ) in_vals;
    get_stats(vals_cast->A,vals_cast->B,vals_cast->num,&(vals_cast->arg.alpha),vals_cast->arg.log_path,vals_cast->data_index);
    pthread_exit(NULL);
}
int main(int argc, char ** argv){
    // declare variables
    htsFile *       input_vcf = NULL ; // can also be an indexed bcf
    bcf_hdr_t *     input_vcf_hdr = NULL ; // just the header
    int             tsv_array[MAXARR][MAXARR] = {};
    int             chr_positions[MAXARR] = {};
    int             dementions[2] = {};
    int  *          chr_positions_p = &chr_positions[0];
    int  *          tsv_arr_p = &tsv_array[0][0]; // pointer to a MAXARR by MAXARR array of ints
    int  *          dementions_p[2] = {&dementions[0],&dementions[1]};
    char * *        head_line[MAXARR] = {};
    int             grid_flag = -1, data_used = 0 , ps_running = 0 ;
    int             ps_done = 0, ps_last_started = 0 , target_index = 0 , test_errors = 0;
    unsigned int    ps_wait = WAIT_PS;
    bool            first_run = true;

    //start program
    parse(argc , argv); // files global ARGS
    if (ARGS.debug_lvl > 0){
        std::cout << "Contig:\t" << ARGS.contig << std::endl
                  << "Threads:\t" << ARGS.threads << std::endl
                  << "Start-End:\t" << ARGS.StartptEndpt[0] << " - " << ARGS.StartptEndpt[1] << std::endl
                  << "Input:\t" << argv[ARGS.bcf_idpt] << std::endl;
    }
    input_vcf = bcf_open(argv[ARGS.bcf_idpt], "r"); // initialise a new bcf object as a file with the path from user

    if(input_vcf == NULL){std::cerr << "Failed to open the file!\n";}
    input_vcf_hdr = bcf_hdr_read(input_vcf); // get the header data
    if(input_vcf_hdr == NULL){std::cerr << "Failed to get header!\n";}

    grid_flag = mk_grid(input_vcf,input_vcf_hdr, chr_positions_p , tsv_arr_p, head_line, dementions_p);
    if ( grid_flag == 0){
        std::cout << "Finished analyzing VCF starting to write\n";
    } else if (grid_flag == -1){
        std::cerr << "Failed to use VCF\n";
    } else if (grid_flag == -2){
        std::cerr << "Couldn't unpack the sample lines\n";
    } else if (grid_flag == -3){
        std::cerr << "Failed to get variant depth for one of the samples\n";
    } else if (grid_flag == -4){
        std::cerr << "Failed to get reference depth for one of the samples\n";
    }
    if(ARGS.debug_lvl > 0){std::cerr << "Obtaind an array of "<< dementions[0] <<" by " << dementions[1]<<std::endl;}
    target_index = write_vals(*head_line,chr_positions,tsv_array,dementions); // writhe the data to the file!!!
    bcf_hdr_destroy(input_vcf_hdr);
    bcf_close(input_vcf);

    data_used = dementions[1];
    do{
        if(!(first_run)) {
            for (int i = ((ps_last_started - ps_running)); i < ps_last_started; ++i) {
                if (DATA[i].state == "done") {
                    if(!(DATA[i].names == "self")){ // if it is not the one that didn't really run, we need to free some space
                        DATA[i].names = std::to_string(chr_positions[target_index])+"-vs-"+ std::to_string(chr_positions[i]);
                    }
                    if (ARGS.debug_lvl > 0) { std::cerr << DATA[i].names << " is DONE!!!" << std::endl; }
                    ps_wait = WAIT_PS;
                    ps_done++;
                    ps_running--;

                } else if (DATA[i].state == "running" ) {
#if DBUG_VV
                    if(ARGS.debug_lvl > 2){
                        std::cerr << "Prossess causing us to brake is:\n"
                                  << DATA[i].names << std::endl
                                  << "With a lower bound:\t" << DATA[i].LB << std::endl
                                  << "A D-stat of:\t" << DATA[i].D_stat << std::endl
                                  << "An upper bound:\t" << DATA[i].UB << std::endl
                                  << "and p-value of:\t" << DATA[i].p_value << std::endl;
                    }else
#endif
                    if (ARGS.debug_lvl > 0) { std::cerr << DATA[i].names << " is still running!!!" << std::endl;}
                    if (ARGS.threads <= ps_running || ps_running + ps_done >= data_used){
                        pthread_join(DATA[i].ps, nullptr); // we don't have any threads, so no need to continue till we have one
                        if( DATA[i].state == "done" ){
                            ps_done++;
                            ps_running--;
                        }
                    }
                    break; // the first running process was found and we can add some, break here
                } else{
                    // you should not get here... so if you do, exit hard
                    return -1;
                }
            }
        } else {
            first_run = false;
        }
#if DBUG_VV
        if(ARGS.debug_lvl > 2){std::cerr << "Data_used:\t"         << data_used
                           << "\nps_running:\t"      << ps_running
                           << "\nps_last_started:\t" << ps_last_started
                           << "\nps_done:\t"         << ps_done << std::endl;}
#endif
        if(ps_done == data_used && ps_running == 0){
            break;
        } else if (ps_done == data_used && ps_running > 0){
            std::cerr <<  "all processes done, but still have threads running....\n";
        }
        if(ps_running < ARGS.threads && ps_last_started < data_used && ps_last_started != target_index){
            if(ARGS.debug_lvl > 0){std::cerr << "Kicking off \t" << std::to_string(chr_positions[target_index]) + "-vs-" + std::to_string(chr_positions[ps_last_started])<< std::endl;}
            DATA[ps_last_started].state = "running";
            DATA[ps_last_started].ptid = ps_last_started;
            DATA[ps_last_started].names = std::to_string(chr_positions[target_index]) + "-vs-" + std::to_string(chr_positions[ps_last_started]);
            // start the child process and add one to the number of prosesses running.
            PS_Data in_vals = ps_fill(tsv_array[target_index], tsv_array[ps_last_started], dementions[0],ps_last_started);
            if (pthread_create(&(DATA[ps_last_started].ps), nullptr, child_ps,(void *) &in_vals) == 0){
                std::cout << "launched thread with tid:\t" << DATA[ps_last_started].ptid << std::endl;
            };
            ps_wait = WAIT_PS;
            ps_last_started++;
            ps_running++;
        } else if (ps_last_started == target_index){
            if(ARGS.debug_lvl > 0){std::cerr << "Don't run vs yourself....\n";}
            DATA[ps_last_started].names = "self";
            DATA[ps_last_started].state = "done";
            DATA[ps_last_started].p_value = 0.0;
            DATA[ps_last_started].D_stat = 0;
            DATA[ps_last_started].UB = 0.0;
            DATA[ps_last_started].LB = 0.0;
            ps_last_started++;
            ps_running++;
        }
    } while (sleep(ps_wait) == 0);
    ///todo report errors
    test_errors = do_correction(ARGS.FWEC,data_used,ARGS.alpha);
    write_values(&DATA[0],data_used,ARGS.outpath);

    return 0;
}

