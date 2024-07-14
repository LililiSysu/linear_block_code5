/*****************************************************************//**
 * \file   LL_OSD_for_Jiabao.h
 * \brief  additional test function of LL_OSD for Jiabao's work 
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"test_common.h"

class test_LL_OSD_for_Jiabao {
public:

	static void LL_OSD_genius_aid_simulation(bool is_first_time, my_double SNR_dB, int first_order, int second_order) {
		int simulation_times = 10000000;		// it is not until 10000 times did the operation conunt become stable, using pass_standard

		// encoding parameter
		const int m = 7;
		const int t = 14;
		GF2e<m>::init();
		BCH<m, t> bch;
		int k = bch.get_k();
		const int n = (1 << m) - 1;
		const int d_min = 2 * t + 1;
		const int k_prime = n - d_min + 1;
		if (is_first_time) {
			bch.print_info();			// BCH(31,16,7)
			cout << "---------------" << endl;
			cout << "SNR(dB)\tFER" << endl;
		}

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = bch.encode(u);
		//cout << "v" << v;

		Matrix<my_double> c = BPSK::modulation(v);

		// channel
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * k / (double)n * SNR));	// BPSK modulation {-1,1} has energy 1
		AWGN::sigma = sigma;

		// decoding parameter
		LLOSD_genius_aid lga(k, k_prime, first_order, second_order);

		// performance recording
		int error_frame = 0;
		const int max_error_frame = 200;
		//Matrix<int> error_frame_ind(1, 200, 'v');

		clock_t start, end;

		start = clock();
		for (int i = 0; i < simulation_times; ++i) {
			Matrix<my_double> recv = AWGN::pass_standard(c);

			if (lga.is_decode_successful(v, recv));
			else {
				// record the error frame

				//error_frame_ind.push_back(i);
				/*cout << "--------- error frame index i = " << i << "--------" << endl;
				if (early_stop_before != early_stop_after) {
					cout << "early stop using" << endl;
					cout << "recv" << recv;
					cout << "used_list = " << total_used_list_num_after - total_used_list_num_before << endl;
				}*/

				error_frame++;
				if (error_frame == max_error_frame) {
					simulation_times = i + 1;
					break;
				}
			}
		}
		end = clock();

		double time_consume = ((double)end - start) / CLOCKS_PER_SEC / (double)simulation_times;
		double pass_AWGN_time = 2.85e-6 / 31 * n;
		cout << SNR_dB << "\t";

		cout << scientific << setprecision(2);
		cout << error_frame / (my_double)simulation_times << endl;
		cout.unsetf(ios::floatfield);
		cout << setprecision(6);
	}
	static void LL_OSD_genius_aid_simulation_multi_SNR(int first_order, int second_order) {
		bool if_output_to_file = true;
		my_double start_SNR = 2.5, SNR_gap = 1, end_SNR = 4.5;

		// set cout into file
		if (if_output_to_file) {
			char file_name[55] = { 0 };
			sprintf_s(file_name, 55, "LL_OSD_genius_aid_Seg(%d,%d).txt", first_order, second_order);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		Matrix<my_double> test_SNR(start_SNR, SNR_gap, end_SNR + 0.01, 'd');
		int len = test_SNR.size();
		bool is_first_time = true;
		for (int i = 0; i < len; ++i) {
			LL_OSD_genius_aid_simulation(is_first_time, test_SNR(i), first_order, second_order);
			is_first_time = false;
		}
	}

	static void TEP_jiabao_test() {

		int list_first;
		int K_first = 10;
		int order_first = 1;
		int** TEP_first = jiabao_TEP::generate(list_first, K_first, order_first);
		for (int i = 0; i < list_first; i++) {
			for (int j = 0; j < order_first; j++) {
				cout << TEP_first[i][j] << ' ';
			}
			cout << endl;
		}

		int list_second;
		int K_second = 20;
		int order_second = 2;
		int** TEP_second = jiabao_TEP::generate(list_second, K_second, order_second);
		for (int i = 0; i < list_second; i++) {
			for (int j = 0; j < order_second; j++) {
				cout << TEP_second[i][j] << ' ';
			}
			cout << endl;
		}
		cout << "---------" << endl;

		int** TEP_combine = jiabao_TEP::combine(TEP_first, list_first, order_first, TEP_second, list_second, order_second);
		int list_total = list_first * list_second + list_first + list_second;
		int order_total = order_first + order_second;

		for (int i = 0; i < list_total; i++) {
			for (int j = 0; j < order_total; j++) {
				cout << TEP_combine[i][j] << ' ';
			}
			cout << endl;
		}

		jiabao_TEP::delete_matrix(TEP_first, list_first);
		jiabao_TEP::delete_matrix(TEP_second, list_second);
		jiabao_TEP::delete_matrix(TEP_combine, list_total);
	}
	static void TEP_jiabao_simple_look_test() {
		int list_first;
		int K_first = 6;
		int order_first = 2;
		int** TEP_first = jiabao_TEP::generate(list_first, K_first, order_first);
		for (int i = 0; i < list_first; i++) {
			for (int j = 0; j < order_first; j++) {
				cout << TEP_first[i][j] << '\t';
			}
			cout << endl;
		}

		jiabao_TEP::delete_matrix(TEP_first, list_first);
	}
	static void TEP_generate_segmentation_test() {
		int list_total;
		int k = 4;
		int order_first = 1;
		int K = 6;
		int order_second = 1;

		int** TEP = jiabao_TEP::generate_segmentation(k, order_first, K, order_second, list_total);
		cout << "-----------" << endl;
		cout << "list_total (outside) = " << list_total << endl;
		int order_total = order_first + order_second;
		for (int i = 0; i < list_total; i++) {
			for (int j = 0; j < order_total; j++) {
				cout << TEP[i][j] << '\t';
			}
			cout << endl;
		}

		jiabao_TEP::delete_matrix(TEP, list_total);
	}
	static void generate_order_pattern_test() {
		int order_first = 1;
		int order_second = 2;

		int sum_order_pattern = (order_first + 1) * (order_second + 1);
		Matrix<int> pattern(sum_order_pattern, 2, '0');			// the length of pattern is only 2

		int pattern_col = 0;
		for (int i = 0; i <= order_first; i++) {
			for (int j = 0; j <= order_second; j++) {
				pattern(pattern_col, 0) = i;
				pattern(pattern_col, 1) = j;
				pattern_col++;
			}
		}
		cout << "pattern" << pattern;
	}
};

