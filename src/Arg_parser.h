//
// Created by sheepless on 5/25/19.
//

#ifndef MAIN_ARG_PARSER_H
#define MAIN_ARG_PARSER_H
#include <iostream>
#include <string>
struct user_arguments{
    int          bcf_idpt        = 0       ;
    std::string  contig          = ""      ;
    int          StartptEndpt[2] = {-1,-1} ;
    int          target          = -1      ;
    int          threads         = 5       ;
    int          min_not_null    = 10      ;
    std::string  outpath         = "./outs/"      ;
    std::string  ref_fmt_flag    = "RO"    ;
    std::string  var_fmt_flag    = "AO"    ;
    std::string  log_path        = "./logs/"    ;
    float        alpha           = 0.05    ;
    int          FWEC            = 0       ;
    int          debug_lvl       = 0       ;
};

user_arguments ARGS;
#endif //MAIN_ARG_PARSER_H
