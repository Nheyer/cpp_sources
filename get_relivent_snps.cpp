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

int write_vals(std::string path, char* samp[MAXARR], int positions[MAXARR], int data[MAXARR][MAXARR], int fill[2]){
    std::ofstream file;
    file.open(path);
    file << "Samples";
    for (int k = 0; k < fill[1]; ++k) {
        file << "\t" << positions[k];
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
    return 0;
}



// function to calculate the disequalibriam coefficient
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





int main(int argc, char ** argv){
    // declare variables
    user_arguments arguments;
    htsFile *   input_vcf = NULL ; // can also be an indexed bcf
    bcf_hdr_t * input_vcf_hdr = NULL ; // just the header
    int         tsv_array[MAXARR][MAXARR] = {{}};
    int         chr_positions[MAXARR] = {};
    int         dementions[2] = {};
    int  *      chr_positions_p = &chr_positions[0];
    int  *      tsv_arr_p = &tsv_array[0][0]; // pointer to a MAXARR by MAXARR array of ints
    int  *      dementions_p[2] = {&dementions[0],&dementions[1]};
    char * *    head_line[MAXARR] = {};
    int         grid_flag = -1 , error_code = 0;

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
    error_code = grid_flag;
    if(DBUG){std::cerr << "Obtaind an array of "<< dementions[0] <<" by " << dementions[1]<<std::endl;}

    write_vals( arguments.outpath, *head_line,chr_positions,tsv_array,dementions); // writhe the data to the file!!!

    bcf_hdr_destroy(input_vcf_hdr);
    bcf_close(input_vcf);
    return error_code;
}

