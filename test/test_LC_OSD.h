/*****************************************************************//**
 * \file   LC_OSD_test.h
 * \brief  pick up LC OSD
 * 
 * \author 26259
 * \date   March 2024
 *********************************************************************/

#pragma once

#include"test_common.h"

extern int frame_max;
extern int frame_min;
extern int frame_error_max;
extern my_double SNR_dB_start;
extern my_double SNR_dB_step;
extern my_double SNR_dB_end;

class test_LC_OSD {
public:
	static bool if_output_to_file;

	static void init_frame_SNR() {

		bool run_performance = true;
		if_output_to_file = false;

		if (run_performance) {
			frame_max = 50000000;			// 50000000
			frame_min = 10000;			// 10000
			frame_error_max = 200;		// 200
			SNR_dB_start = 1;		// 0
			SNR_dB_step = 0.5;		// 1,  0.5 
			SNR_dB_end = 3.5;		// 8, 7, 4.5, 4
		}
		else {
			frame_max = 10000;
			frame_min = 10000;
			frame_error_max = 1000000;
			SNR_dB_start = 5;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
	}

	static int LC_OSD_FER() {
		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 7;
		const int t = 10;
		// LC_OSD parameter
		const int delta = 8;
		const int Lmax = 16384;


		//// BCH code parameter
		//const int m = 6;
		//const int t = 3;
		//// LC_OSD parameter
		//const int delta = 6;
		//const int Lmax = 1024;

		GF2e<m>::init();
		eBCH<m, t> ebch;
		Matrix<GF2> G = ebch.get_generator_matrix();
		Matrix<GF2> H = ebch.get_parity_matrix();
		int n = ebch.get_n();
		int k = ebch.get_k();
		int d = ebch.get_d();

		/*cout << "G = " << endl;
		G.print();
		cout << "H = " << endl;
		H.print();*/

		cout << "-------------------------" << endl << endl;

		Matrix<GF2> u(1, k);											// message
		for (int i = 0; i < k; ++i) {
			u[i] = my::rand_01();
		}
		Matrix<GF2> c = u * G;									// transmitted codeword
		cout << "c = " << endl;
		c.print();
		cout << "-------------------------" << endl << endl;

		cout << "eBCH(" << n << ", " << k << ", " << d << ")" << endl;


		Matrix<my_double> SNR_dB_vec(SNR_dB_start, SNR_dB_step, SNR_dB_end, 'd');
		//Vector_ext::equally_distant(SNR_dB_start, SNR_dB_end, SNR_dB_step, SNR_dB_vec);// 	SNR_dB_vec[0] = 6;
		//cout << "SNR_dB_vec" << SNR_dB_vec;


		LC_OSD_r lc_osd(H, d, delta);
		cout << "LC_OSD_r(" << delta << ", " << Lmax << ")" << endl;


		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			sprintf_s(file_name, 155, "(%d,%d,%d)eBCH_LC_OSD(%d,%d).txt", n, k, d, delta, Lmax);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}
		else {
			printf("\n");
		}
		printf("SNR(dB)   FER          Errors    Frames        GF2_ope      float_ope    ave_list     ML_list      max_list_rate    time(s)\n");
		clock_t start, end;

		my::set_seed_adv(1);

		for (int i = 0; i < SNR_dB_vec.size(); ++i) {

			double SNR = (double)pow(10, SNR_dB_vec[i] / 10);
			double sigma = sqrt(n / (2 * k * SNR));
			int frame_error = 0;
			int frame_finished = 0;
			GF2_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0; 

			lc_osd.set_sigma(sigma);
			lc_osd.total_used_list_num = 0;
			lc_osd.ML_used_list_num = 0;
			lc_osd.times_reach_max_list_num = 0;

			start = clock();
			while (frame_finished < frame_max) {

				/*Matrix<my_double, 1, n> r({
					-0.929966,-0.777951,-0.852718,1.35274,0.0624572,-0.267671,-0.690881,1.07161,-0.37672,-1.59395,-0.84877,-1.48822,-1.39876,-0.137612,0.991339,
				});*/										// recieved word


				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = (c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}


				Matrix<GF2> c_hat = lc_osd.decode_v(r, Lmax);

				if (c_hat != c) {
					frame_error++;
					//cout << "decode error r";
					//r.print();
				}

				frame_finished++;
				if (frame_error >= frame_error_max && frame_finished >= frame_min) {
					break;
				}
			}

			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC;
			//cout << "time_consume = " << time_consume << " s" << endl;

			printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E        %0.3E\n", \
				(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
				(double)GF2_auxiliary_storage::operation_number / frame_finished, \
				(double)my_double_auxiliary_storage::operation_number / frame_finished, \
				(double)lc_osd.total_used_list_num / frame_finished, \
				(double)lc_osd.ML_used_list_num / frame_finished, \
				(double)lc_osd.times_reach_max_list_num / frame_finished, \
				time_consume / frame_finished);
		}

		return 0;
	}
};

bool test_LC_OSD::if_output_to_file = false;
