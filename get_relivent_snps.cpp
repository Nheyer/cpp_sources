//
// Created by sheepless on 3/19/19.
//


#include "get_relivent_snps.h"


user_arguments parse(int Arglength, char ** ArgComands ){
    int i = 1;
    int middle = 0, range = 0;
    user_arguments save;

    while (i < Arglength) {
        std::string arg_i = (std::string) ArgComands[i];
        if (arg_i == "-i") {
            save.bcf_idpt = i + 1;
            i += 2; // skip one b/c we used it here.
        } else if (arg_i == "-t") {
            save.target = std::stoi(ArgComands[i + 1]);
            middle = std::stoi(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        } else if (arg_i == "-r") {
            range = std::stoi(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        } else if (arg_i == "-@") {
            save.threads = std::stoi(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        } else if (arg_i == "-c"){
            save.contig = (std::string)ArgComands[i + 1];
            i += 2; // skip one b/c we used it here.
        }else if (arg_i == "-o"){
            save.outpath = (std::string)(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        }else if (arg_i == "-R") {
            save.ref_fmt_flag = ArgComands[i + 1];
            i += 2; // skip one b/c we used it here.
        } else if (arg_i == "-V"){
            save.var_fmt_flag = ArgComands[i + 1];
            i += 2; // skip one b/c we used it here.
        } else if (arg_i == "-n"){
            save.min_not_null = std::stoi(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        }else {
                std::cerr << "Unexpected Input:\t" << arg_i << std::endl;
                i++;
            }
        }

    if(middle < range){
        std::cerr << "Range is greater then the mid point! . . . . .\n Setting start to zero instead\n";
        save.StartptEndpt[0] = 0;
    } else {
        save.StartptEndpt[0] = middle - (range / 2);
    }
    save.StartptEndpt[1] = middle + (range / 2);


    return save;
};

int  get_int_type_fmt(bcf1_t * data, std::string tag , bcf_idpair_t *pair_of_vals ,int nsamp, int * ints){
    int j = 0, tag_id ;
    if(DBUG_V){std::cerr << "nsamp has a value of : \t" << nsamp << std::endl;}
    if (nsamp > 256) {
        std::cerr << "Can't parse over 256 samples!\n";
        return  -1; }
    // figure out the id of the tag the user input
    while (true){
        if(DBUG_V){std::cerr << (std::string)pair_of_vals[j].key << "\t " << tag<< std::endl;}
        if((std::string) pair_of_vals[data->d.fmt[j].id].key ==  tag){
            if(DBUG_V){std::cerr << (std::string)pair_of_vals[data->d.fmt[j].id].key << "\t" << j << std::endl;}
            tag_id = j;
            break;
        }
        j++;
    }
    //if(DBUG_V){std::cerr << data->d.fmt.id<< std::endl;}
    // use the id to get the list of int values for each sample (corresponding to the tag) at that pos.
    for (int i = 0; i < nsamp; i++) {
        if(DBUG_V){std::cerr << ((unsigned int)(data->d.fmt[tag_id].p[i])) % 128  << "\t"; }
        ints[i] = ((unsigned int)(data->d.fmt[tag_id].p[i])) % 128 ;
    }
    return 0;
}


int mk_grid(htsFile * bcf, bcf_hdr_t * hdr,user_arguments * args, int * poss, int * arr, char ** header[MAXARR], int * dem[2]){
    // initialise the data
    int i = 0 , j = 0 , k = 0 , l = 0;
    int cur_conting_pos = 0 , num_not_null = 0;
    int NumSamples = 1 ;
    int var_type[MAXARR] = {} , ref_type[MAXARR] = {}; // make some space for fmt tag values
    int * var_p = &var_type[0]; // pointer to the first pos of the var_type array
    int * ref_p = &ref_type[0]; // pointer to the first pos of the ref_type array
    int end = args->StartptEndpt[1];
    int start = args->StartptEndpt[0];
    bool read_some = false, snp = false;
    unsigned short int store = 0;

    // make some space in the stack
    bcf1_t * line = bcf_init();

    // Get the number of samples in the multisample VCF.
    NumSamples = bcf_hdr_nsamples(hdr);
    if(DBUG){std::cout << "Number of samples found:\t" << NumSamples << std::endl;}
    while (i < NumSamples && i < MAXARR){   // make sure we dont over flow anything but loop through all
        header[i] = &hdr->samples[i]; // rip the sample names from the header
        if(DBUG){std::cerr << *(header[i]) << "\t" << i << std::endl;}
        i++;
    }
    // see if we filled the samples or the numver of strings we have
    if ( i == 255){
        std::cerr << "Filled the header line, were there more then 256?... If so split up the samples\n";
        NumSamples = 256;
    }

    //loop through all positions in the VCF
    while(bcf_read(bcf,hdr,line) == 0){
        if(bcf_hdr_id2name(hdr,line->rid) == args->contig) { //see if it is the right contig
            cur_conting_pos = line->pos + 1 ; // make things more readable, and convert to 1 indexed
            if (cur_conting_pos < start) { // if it is before the start... continue
                continue;
            } else if (cur_conting_pos <= end) {
                if (DBUG) { std::cerr << cur_conting_pos << std::endl; } // testing
                if(bcf_unpack(line,BCF_UN_IND) != 0) { return -2;}// Unpack the data for the individual samples on the lines we care about

                if (get_int_type_fmt(line,args->var_fmt_flag,hdr->id[BCF_DT_ID],NumSamples, var_p ) != 0 ){ return -3;};
                if (get_int_type_fmt(line,args->ref_fmt_flag,hdr->id[BCF_DT_ID],NumSamples, ref_p ) != 0 ){ return -4;};
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
                    if(DBUG_V){std::cout<< var_type[m] <<  "," << ref_type[m] << ":";}
                    if(DBUG){std::cerr << store << "\t";}
                    arr[m + (MAXARR*j)] = store; // add value to the int array
                }
                // we only want lines with at least one snp
                if((num_not_null >= args->min_not_null && snp) or (cur_conting_pos == args->target )) {
                    if(DBUG_V){std::cerr << "adding position " << cur_conting_pos << " and starting next row!\n";}
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

int write_vals(user_arguments * args, char* samp[MAXARR], int positions[MAXARR], int data[MAXARR][MAXARR], int fill[2]){
    std::ofstream file;
    file.open("./outs" + args->outpath + "_raw.tsv");
    int target_index = 0;
    file << "Samples";
    for (int k = 0; k < fill[1]; ++k) {
        file << "\t" << positions[k];
        if(positions[k] == args->target){
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
    std::cerr <<" \t" << "+"  << " \t " << "-"  << std::endl
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

int permutation_test(bool A[MAXARR], bool B[MAXARR], int max , std::string path_to_log , out_data * report){
    int i ;
    double p_num = 0.0 ;
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
        if (abs(D_loop) > abs(report->D_stat)){
            p_num++;
        }
        if(logging){
            perm_log << D_loop << std::endl;
        }
    }
    report->p_value = (float) (p_num / (double) i);
    if(DBUG || !(logging)) {
        std::cerr << report->p_value << "\t" << report->D_stat << "\t" << max << std::endl;
    }
    if(logging){
        perm_log.close();
    }
    return 0;
}

int boot_strap(bool A[MAXARR], bool B[MAXARR], int max, float alpha, std::string path_to_log ,out_data * report){
    static float resampled_disequalibriums[RESAMPLES] = {};
    int index ;
    bool A_loop[MAXARR] = {};
    bool B_loop[MAXARR] = {};
    bool logging = false;
    float lower = 0.0 , upper = 0.0;
    std::ofstream boot_log ;
    if(path_to_log != ""){ // see if we want to log the alt distrabution
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
    if(DBUG || !(logging)) {     // see if we should print to screen if so, do
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

out_data get_stats(int a_raw[MAXARR], int b_raw[MAXARR], int init_num , float * alpha , std::string path_to_log , std::string name){
    out_data pair_vals ;
    bool A[MAXARR] = {};
    bool B[MAXARR] = {};
    int num = 0;
    bool * A_p = &A[0];
    bool * B_p = &B[0];
    int  * num_p = &num ;
    int i ;

    cleaning_func(a_raw,b_raw,init_num, A_p, B_p , num_p);
    pair_vals.D_stat = gammitic_disequalibrium(A,B,num);
    print_table(A,B,num); // prints grapphical table of data
    srand(time_t(NULL));
    boot_strap(A, B, num, * alpha , path_to_log  + "/logs/" + name, &pair_vals); // do bootstrap
    permutation_test(A, B, num , path_to_log + "/logs/" + name, &pair_vals); // do permutation test
    pair_vals.reject = (pair_vals.p_value < *alpha);
    pair_vals.state = "done";
    pair_vals.names = name;
    return pair_vals;
}
int write_values(out_data data[MAXARR], int max, std::string out_path){
    std::ofstream pairwise_data;
    pairwise_data.open("outs/" + out_path + ".summary.tsv");
    pairwise_data << "Pairing \t Lower_Bound \t Upper_Bound \t Sample_Disequilibrium \t p-value \t Decision";
    for (int i = 0; i < max ; ++i) { // loop through all values and print them
        pairwise_data << std::endl
                      << data[i].names   << "\t"
                      << data[i].LB      << "\t"
                      << data[i].UB      << "\t"
                      << data[i].D_stat  << "\t"
                      << data[i].p_value << "\t"
                      << data[i].reject;
    }
    pairwise_data.close();
    return 0;
}


int main(int argc, char ** argv){
    // declare variables
    out_data data[MAXARR];
    user_arguments arguments;
    htsFile *   input_vcf = NULL ; // can also be an indexed bcf
    bcf_hdr_t * input_vcf_hdr = NULL ; // just the header
    int         tsv_array[MAXARR][MAXARR] = {};
    int         chr_positions[MAXARR] = {};
    int         dementions[2] = {};
    int  *      chr_positions_p = &chr_positions[0];
    int  *      tsv_arr_p = &tsv_array[0][0]; // pointer to a MAXARR by MAXARR array of ints
    int  *      dementions_p[2] = {&dementions[0],&dementions[1]};
    char * *    head_line[MAXARR] = {};
    int         grid_flag = -1, data_used = 0 , ps_running = 0 ;
    int         ps_done = 0, ps_last_started = 0 , target_index = 0;
    bool        first_run = true;

    //start program
    arguments = parse(argc , argv);
    if (DBUG){
        std::cout << "Contig:\t" << arguments.contig << std::endl
                  << "Threads:\t" << arguments.threads << std::endl
                  << "Start-End:\t" << arguments.StartptEndpt[0] << " - " << arguments.StartptEndpt[1] << std::endl
                  << "Input:\t" << argv[arguments.bcf_idpt] << std::endl;
    }
    input_vcf = bcf_open(argv[arguments.bcf_idpt], "r"); // initialise a new bcf object as a file with the path from user

    if(input_vcf == NULL){std::cerr << "Failed to open the file!\n";}
    input_vcf_hdr = bcf_hdr_read(input_vcf); // get the header data
    if(input_vcf_hdr == NULL){std::cerr << "Failed to get header!\n";}

    grid_flag = mk_grid(input_vcf,input_vcf_hdr, & arguments, chr_positions_p , tsv_arr_p, head_line, dementions_p);
    // cheack for errors in the above functinon
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
    if(DBUG){std::cerr << "Obtaind an array of "<< dementions[0] <<" by " << dementions[1]<<std::endl;}
    target_index = write_vals( &arguments, *head_line,chr_positions,tsv_array,dementions); // writhe the data to the file!!!
    bcf_hdr_destroy(input_vcf_hdr);
    bcf_close(input_vcf);
    /// todo  turn this into its own function [DRY]
    data_used = dementions[1];
    do{
        if(!(first_run)) {
            for (int i = ((ps_last_started - ps_running) - 1); i < ps_running; ++i) {
                if (data[i].state == "done") {
                    if (DBUG) { std::cerr << data[i].names << " is DONE!!!" << std::endl; }
                    data[i].names = std::to_string(chr_positions[target_index])
                                    + "-vs-"
                                    + std::to_string(chr_positions[i]);
                    ps_done++;
                    ps_running--;

                } else if (data[i].state == "running") {
                    if (DBUG) { std::cerr << data[i].names << "is still running!!!" << std::endl; }
                    break; // the first running process was found , break here
                }
            }
        } else {
            first_run = false;
        }
        if(DBUG){std::cerr << "Data_used:\t"         << data_used
                           << "\nps_running:\t"      << ps_running
                           << "\nps_last_started:\t" << ps_last_started
                           << "\nps_done:\t"         << ps_done << std::endl;}
        if(ps_done == data_used && ps_running == 0){
            break;
        } else if (ps_done == data_used && ps_running > 0){
            std::cerr <<  "all processes done, but still have threads running....\n";
        }
        if(ps_running < arguments.threads && ps_last_started < data_used && ps_last_started != target_index){
            data[ps_last_started].state = "running";
            data[ps_last_started] = get_stats(tsv_array[target_index],
                                              tsv_array[ps_last_started],
                                              dementions[0],
                                              &arguments.alpha,
                                              arguments.log_path,
                                              std::to_string(chr_positions[target_index])
                                              + "-vs-"
                                              + std::to_string(chr_positions[ps_last_started]));
            ps_last_started++;
            ps_running++;
            if(DBUG){std::cerr << data[ps_last_started].state << std::endl;}
        } else if (ps_last_started == target_index){
            if(DBUG){std::cerr << "Don't run vs yourself....\n";}
            data[ps_last_started].names = "self";
            data[ps_last_started].state = "done";
            data[ps_last_started].p_value = 0.0;
            data[ps_last_started].D_stat = 0;
            data[ps_last_started].UB = 0.0;
            data[ps_last_started].LB = 0.0;
            ps_last_started++;
            ps_done++;
        }
    } while (sleep(20) == 0);
    write_values(&data[0],data_used,arguments.outpath);

    return 0;
}

