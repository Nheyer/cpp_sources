//
// Created by sheepless on 3/19/19.
//

#ifndef RECEARCH_GET_RELIVENT_SNPS_H
#define RECEARCH_GET_RELIVENT_SNPS_H
#define DBUG true
#define DBUG_V false
#define MAXARR 256
#define RESAMPLES 100000


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
struct out_data{
    pthread_t   process;
    std::string names   = "";
    float       LB      = 0.0;
    double      D_stat  = 0.0;
    float       UB      = 0.0;
    float       p_value = 0.0;
    bool        reject  = false;
    std::string state   = "";
};
#endif //RECEARCH_GET_RELIVENT_SNPS_H
