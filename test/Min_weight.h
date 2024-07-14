/*****************************************************************//**
 * \file   Min_weight.h
 * \brief  test for the minimum weight of a code. Count the minnimum weight codewords of a code.
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"test_common.h"


class test_Min_weight{
public:

	static void BCH_dual_test() {
		const int m = 6;
		GF2e<m>::init();
		BCH<m, 7> bch;

		bch.print_info();
		Matrix<GF2> G = bch.get_generator_matrix();
		cout << "G" << G;

		Matrix<GF2> H = bch.get_parity_matrix();
		cout << "H" << H;

		cout << "G.multiply_transpose_of(H).isZero() = " << G.multiply_transpose_of(H).isZero() << endl;
	}
};
