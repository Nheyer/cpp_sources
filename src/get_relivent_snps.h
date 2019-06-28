//
// Created by nheyer on 3/19/19.
//

#ifndef RECEARCH_GET_RELIVENT_SNPS_H
#define RECEARCH_GET_RELIVENT_SNPS_H
#define DBUG            false
#define DBUG_V          false
#define DBUG_VV         false
#define MAXARR          256
#define RESAMPLES       10000000
#define WAIT_PS         1
// define BONFERRONI      0
// define HOLM_BONFERRONI 1


//---packages---//
#include "Arg_parser.cpp"
#include "IO.cpp"
#include "htslib/vcf.h"
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <stdint.h>
#include <pthread.h>
#include <zconf.h>
#include <cmath>


//---structures---//


struct PS_Data{
    int A[MAXARR] = {};
    int B[MAXARR] = {};
    int num = 0;
    user_arguments arg;
    int data_index = -1;
};
PS_Data ps_fill(int *A, int *B, int num,  int data_p){
    PS_Data rt;
    rt.num = num;
    rt.arg = ARGS;
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
