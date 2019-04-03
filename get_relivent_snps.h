//
// Created by sheepless on 3/19/19.
//

#ifndef RECEARCH_GET_RELIVENT_SNPS_H
#define RECEARCH_GET_RELIVENT_SNPS_H
#define DBUG false
#define DBUG_V false
#define MAXARR 256
#define RESAMPLES 10000000
#define BOOT 10000

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

//---structures---//
struct user_arguments{
    int bcf_idpt = 0;
    std::string contig = "";
    int StartptEndpt[2] = {0,0};
    int target  = 0;
    int threads = 5;
    int min_not_null = 10;
    std::string outpath;
    std::string ref_fmt_flag = "RO";
    std::string var_fmt_flag = "AO";

};

std::string get_contig(user_arguments * st){
    return st->contig;
}

#endif //RECEARCH_GET_RELIVENT_SNPS_H
