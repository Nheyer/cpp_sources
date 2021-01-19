

#ifndef MAIN_ARG_PARSER_H
#define MAIN_ARG_PARSER_H
#include <iostream>
#include <string>
struct user_arguments{
    int          bcf_idpt        = 0       ;
    std::string  contig          = ""      ;
    int          StartptEndpt[2] = {-1,-1} ;
    int          target          = -2      ; // -1 reserved for samp name targeting
    int          threads         = 5       ;
    int          min_not_null    = 10      ;
    std::string  outpath         = "./"    ;
    const char *  ref_fmt_flag    = "RO"   ;
    const char *  var_fmt_flag    = "AO"   ;
    std::string  log_path        = "./"    ;
    float        alpha           = 0.05    ;
    int          FWEC            = 0       ;
    int          debug_lvl       = 0       ;
    bool         pairwise_targets= false   ;
    std::string  regx_match                ;

};

user_arguments ARGS;
#endif //MAIN_ARG_PARSER_H
