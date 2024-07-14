/*****************************************************************//**
 * \file   test_common.h
 * \brief  including all pakages of the project
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"../GF/GFq.h"					// including <cmath>, <vector>
#include"../GF/Statistics.h"
#include"../my_lib/my_double.h"
#include"../my_lib/Vector_static.h"
#include"../channel/channel.h"
#include"../channel/simulation.h"		// containing all .h file in "code"
#include"../reprocess/reprocess_include_all.h"
#include"../Viterbi/Viterbi_include_all.h"
#include<ctime>
#include<stdio.h>
#include<unordered_map>
#include<unordered_set>
#include<fstream>
using namespace std;

int frame_max;
int frame_min;
int frame_error_max;
my_double SNR_dB_start;	
my_double SNR_dB_step;
my_double SNR_dB_end;