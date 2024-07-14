#pragma once
/*****************************************************************//**
 * \file   Hybrid_Chase2_OSD.h
 * \brief  simple hybrid of chase and osd
 * 
 * \author 26259
 * \date   March 2023
 *********************************************************************/
#include"LRP.h"
#include"MRIP.h"

class Hybrid_Chase2_OSD {
public:

	static int num_Chase2_stop;
	static int num_early_stop;							// only means for early stop at Chase2 decoding
	static unsigned long long total_used_list_num;		// counting used list num
	static bool is_ML;

	template<class code_type>
	static Matrix<GF2> decode_v(code_type& code_instance, const Matrix<my_double>& r, int eta = -1, int order = -1) {
		is_ML = false;
		unsigned long long total_used_list_num_before = Chase::total_used_list_num;
		Matrix<GF2> can_Chase2_v = Chase::decode2_v(code_instance, r, eta);
		total_used_list_num += Chase::total_used_list_num - total_used_list_num_before;
		if (Chase::is_ML == true) {
			is_ML = true; 
			num_Chase2_stop++;
			num_early_stop++;
			return can_Chase2_v;
		}

		OSD_r osd(code_instance.get_parity_matrix(), code_instance.get_d());
		Matrix<GF2> can_OSD_v = osd.decode_v(r, order);
		total_used_list_num += osd.total_used_list_num;

		if (osd.is_ML == true) {
			is_ML = true;
			num_early_stop++;
			return can_OSD_v;
		}

		// hybrid them
		return Chase::lambda_min < osd.lambda_min ? can_Chase2_v : can_OSD_v;		// problematic
	}
};

int Hybrid_Chase2_OSD::num_Chase2_stop = 0;
int Hybrid_Chase2_OSD::num_early_stop = 0;
unsigned long long Hybrid_Chase2_OSD::total_used_list_num = 0;
bool Hybrid_Chase2_OSD::is_ML = false;

