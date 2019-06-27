//
// Created by sheepless on 5/25/19.
//

#include "Arg_parser.h"

void print_help(){
    std::cerr << "~~ Usage ~~ \n\n"
              << "-h\t          \t Print this help message\n"
              << "----------Mandatory Flags----------\n"
              << "-i\t<PATH>    \t Path to ether a vcf file or a indexed bcf file\n"
              << "-t\t<INT>     \t The position that is to be tested against for linkage\n"
              << "-c\t<STR>     \t The contig of the positions\n"
              << "-r\t<INT>     \t The range around the center to test\n"
              << "-A\t<INT><INT>\t instead of testing ageist a target compare all\n"
              << "  \t          \t between a range pairwise (dont use -t or -r)\n"
              << "----------Optional Flags----------\n\n"
              << "-n\t<INT>     \t Minimum number of none null values in non-target needed to test it \t default == 10\n"
              << "-V\t<STR>     \t Format flag for total varent depth\t default == AO\n"
              << "-R\t<STR>     \t Format flag for total refference depth\t default == RO\n"
              << "-o\t<PATH>    \t Base file name, and path to output data table \t default == ./\n"
              << "-L\t<PATH>    \t Directory to place log files\t default == ./\n"
              << "-C\t<INT>     \t type of strong error wise correction to use HOLM_BONFERRONI = 1, default == BONFERRONI = 0 \n"
              << "-@\t<INT>     \t Maximum number of cores to use for the permutation and bootstrap tests\n"
              << "-v\t          \t print first level of debug info, remake with the DEBUG_V and Debug_VV \n"
              << "  \t          \t in the header for more info\n";
    exit(-1);
}

void parse(int Arglength, char ** ArgComands ){
    int i = 1;
    int middle = 0, range = 0;
    bool need_to_set_start_end = true;

    while (i < Arglength) {
        std::string arg_i = (std::string) ArgComands[i];
        if (arg_i == "-i") {ARGS.bcf_idpt = i + 1; i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-t") {
            ARGS.target = std::stoi(ArgComands[i + 1]);
            middle = std::stoi(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        }
        else if (arg_i == "-A") {
            need_to_set_start_end = false;
            ARGS.pairwise_targets = true;
            ARGS.StartptEndpt[0] = std::abs(std::stoi(ArgComands[i + 1]));
            ARGS.StartptEndpt[1] = std::abs(std::stoi(ArgComands[i + 2]));
            i += 3;
        }
        else if (arg_i == "-r") {range             = std::stoi(ArgComands[i + 1])    ; i += 2;}
        else if (arg_i == "-@") {ARGS.threads      = std::stoi(ArgComands[i + 1])    ; i += 2;}
        else if (arg_i == "-c") {ARGS.contig       = (std::string)ArgComands[i + 1]  ; i += 2;}
        else if (arg_i == "-o") {ARGS.outpath      = (std::string)(ArgComands[i + 1]); i += 2;}
        else if (arg_i == "-R") {ARGS.ref_fmt_flag = ArgComands[i + 1]               ; i += 2;}
        else if (arg_i == "-V") {ARGS.var_fmt_flag = ArgComands[i + 1]               ; i += 2;}
        else if (arg_i == "-n") {ARGS.min_not_null = std::stoi(ArgComands[i + 1])    ; i += 2;}
        else if (arg_i == "-a") {ARGS.alpha        = std::stoi(ArgComands[i + 1])    ; i += 2;}
        else if (arg_i == "-L") {ARGS.log_path     = ArgComands[i + 1]               ; i += 2;}
        else if (arg_i == "-C") {ARGS.FWEC         = std::stoi(ArgComands[i + 1])    ; i += 2;}
        else if (arg_i[1] == 'v'){ARGS.debug_lvl   = arg_i.length() - 1              ; i += 2;}
        else if (arg_i == "-h") {print_help();}
        else {std::cerr << "Unexpected Input:\t" << arg_i << std::endl; i++;}
    }
    if(!need_to_set_start_end){
        return;
    }else if(ARGS.target < 0 || range < 0 || ARGS.contig == ""){
        std::cerr << "Missing at least one of the required arguments:\n";
        print_help();
    } else if(middle < range){
        std::cerr << "Range is greater then the mid point! . . . . .\n Setting start to zero instead\n";
        ARGS.StartptEndpt[0] = 0;
    } else {
        ARGS.StartptEndpt[0] = middle - (range / 2);
    }
    ARGS.StartptEndpt[1] = middle + (range / 2);
}

