//
// Created by sheepless on 3/19/19.
//

#ifndef RECEARCH_GET_RELIVENT_SNPS_H
#define RECEARCH_GET_RELIVENT_SNPS_H
#define DBUG      true
#define DBUG_V    true
#define DBUG_VV   false
#define MAXARR    256
#define RESAMPLES 10000000
#define WAIT_PS   5

//---packages---//
#include <iostream>
#include <string>
#include "htslib/vcf.h"
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <stdint.h>
#include <pthread.h>
#include <zconf.h>

void print_help(){
    std::cerr << "~~ Usage ~~ \n\n"
              << "-h\t\tPrint this help message\n"
              << "----------Mandatory Flags----------\n"
              << "-i\t<PATH>\t Path to ether a vcf file or a indexed bcf file\n"
              << "-t\t<INT>\t The position that is to be tested against for linkage\n"
              << "-c\t<STR>\t The contig of the positions\n"
              << "-r\t<INT>\t The range around the center to test\n"
              << "\n----------Optional Flags----------\n\n"
              << "-n\t<INT>\t Minimum number of none null values in non-target needed to test it \t default == 10\n"
              << "-V\t<STR>\t Format flag for total varent depth\t default == AO\n"
              << "-R\t<STR>\t Format flag for total refference depth\t default == RO\n"
              << "-L\t<PATH>\t Directory to place log files\t default == ./\n"
              << "-@\t<INT>\t Maximum number of cores to use for the permutation and bootstrap tests\n";
    exit(-1);
}
//---structures---//
struct user_arguments{
    int         bcf_idpt        = 0;
    std::string contig          = "";
    int         StartptEndpt[2] = {0,0};
    int         target          = 0;
    int         threads         = 5;
    int         min_not_null    = 10;
    std::string outpath         = "";
    std::string ref_fmt_flag    = "RO";
    std::string var_fmt_flag    = "AO";
    std::string log_path        = "./";
    float       alpha           = 0.05;

};

user_arguments parse(int Arglength, char ** ArgComands ){
    int i = 1;
    int middle = 0, range = 0;
    user_arguments save;

    while (i < Arglength) {
        std::string arg_i = (std::string) ArgComands[i];
        if (arg_i == "-i") { save.bcf_idpt = i + 1; i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-t") {
            save.target = std::stoi(ArgComands[i + 1]);
            middle = std::stoi(ArgComands[i + 1]);
            i += 2; // skip one b/c we used it here.
        }
        else if (arg_i == "-r") {range             = std::stoi(ArgComands[i + 1]);     i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-@") {save.threads      = std::stoi(ArgComands[i + 1]);     i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-c") {save.contig       = (std::string)ArgComands[i + 1];   i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-o") {save.outpath      = (std::string)(ArgComands[i + 1]); i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-R") {save.ref_fmt_flag = ArgComands[i + 1];                i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-V") {save.var_fmt_flag = ArgComands[i + 1];                i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-n") {save.min_not_null = std::stoi(ArgComands[i + 1]);     i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-a") {save.alpha        = std::stoi(ArgComands[i + 1]);     i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-L") {save.log_path     = ArgComands[i + 1];                i += 2;} // skip one b/c we used it here.
        else if (arg_i == "-h") {print_help(); ;}
        else {
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

struct out_data{
    std::string names   = "";
    float       LB      = 0.0;
    double      D_stat  = 0.0;
    float       UB      = 0.0;
    float       p_value = 0.0;
    bool        reject  = false;
    std::string state   = "";
    pthread_t   ps;
    int         ptid    = -1;
};
struct PS_Data{
    int A[MAXARR] = {};
    int B[MAXARR] = {};
    int num = 0;
    user_arguments arg;
    int data_index = -1;
};
PS_Data ps_fill(int *A, int *B, int num, user_arguments &arg, int data_p){
    PS_Data rt;
    rt.num = num;
    rt.arg = arg;
    rt.data_index = data_p;
    for (int i = 0; i < num ; ++i) {
        rt.A[i] = A[i];
        rt.B[i] = B[i];
    }
    return rt;
}
//-- GLOBAL-VARS --//
static out_data DATA[MAXARR];
#endif //RECEARCH_GET_RELIVENT_SNPS_H
