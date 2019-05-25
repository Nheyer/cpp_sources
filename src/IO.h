//
// Created by sheepless on 5/25/19.
//

#ifndef MAIN_IO_H
#define MAIN_IO_H

#include <fstream>
#include <iostream>
#include <string>
#include <pthread.h>

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



#endif //MAIN_IO_H
