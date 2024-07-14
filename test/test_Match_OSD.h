/*****************************************************************//**
 * \file   test_Match_OSD.h
 * \brief  new idea of matching TEP for OSD
 * 
 * \author 26259
 * \date   March 2024
 *********************************************************************/

#pragma once

#include"test_common.h"


class test_Match_OSD {
public:
	static int Match_OSD_initialize() {

		// BCH code parameter
		const int m = 4;
		const int t = 2;
		// Match OSD parameter
		int order = 1;
		int delta = 2;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int n = G.col();
		int k = G.row();
		int d = bch.get_d();

		//Match_OSD osd(H, d, order, delta, false);
		Match_OSD_standard osd(G, H, d, order, delta);

		Matrix<my_double> r(1, n);
		my_double sigma = 10;								// nearly without noise
		for (int p = 0; p < n; ++p) {
			// transmit all 0 vector
			r[p] = 1 + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
		}

		osd.solve(r);

		return 0;
	}

	static bool if_output_to_file;
	static void init_frame_SNR() {

		bool run_performance = true;
		if_output_to_file = false;

		if (run_performance) {
			frame_max = 50000000;			// 50000000
			frame_min = 10000;			// 10000
			frame_error_max = 200;		// 200
			SNR_dB_start = 0;		// 0
			SNR_dB_step = 0.5;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
		else {
			frame_max = 10000;
			frame_min = 10000;
			frame_error_max = 1000000;
			SNR_dB_start = 6;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
	}

	// OSD test FER
	static int Match_OSD_FER()
	{
		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 4;
		const int t = 2;
		// OSD parameter
		int order = 1;
		int delta = 2;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int n = G.col();
		int k = G.row();
		int d = bch.get_d();

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

		cout << "BCH(" << n << ", " << k << ", " << d << ")" << endl;

		Matrix<my_double> SNR_dB_vec(SNR_dB_start, SNR_dB_step, SNR_dB_end, 'd');
		//Vector_ext::equally_distant(SNR_dB_start, SNR_dB_end, SNR_dB_step, SNR_dB_vec);// 	SNR_dB_vec[0] = 6;
		//cout << "SNR_dB_vec" << SNR_dB_vec;



		Match_OSD_standard osd(G, H, d, order, delta);
		//Match_OSD osd(G, H, d, order, delta);

		//osd.ET_type = 2;		// ET_type to be updated in IBU and IBUc of OSD

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			if (osd.type == 21)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_Match_OSD(%d, %d).txt", n, k, d, order, delta);		// it semms harder and no state decrease.
			else if (osd.type == 24)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_Match_OSD_dual(%d, %d).txt", n, k, d, order, delta);
			else
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_unknown(%d).txt", n, k, d, order);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}
		else {
			printf("\n");
		}
		//printf("SNR(dB)   FER          Errors    Frames        GF2_ope      float_ope    time(s)\n");
		printf("SNR(dB)   FER          Errors    Frames        GF2_ope      float_ope    ave_list     ET_rate      time(s)\n");
		//printf("________________________________________________________________________________\n\n");
		clock_t start, end;

		my::set_seed_adv(1);

		for (int i = 0; i < SNR_dB_vec.size(); ++i) {

			double SNR = (double)pow(10, SNR_dB_vec[i] / 10);
			double sigma = sqrt(n / (2 * k * SNR));
			int frame_error = 0;
			int frame_finished = 0;
			GF2_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0;
			unsigned long long TEP_total = 0;
			int ET_frames_total = 0;

			//osd.sigma = sigma;						// ET, to be extended in 'OSD_v4_dual'

			start = clock();
			while (frame_finished < frame_max) {

				/*Matrix<my_double, 1, n> r({
					-0.929966,-0.777951,-0.852718,1.35274,0.0624572,-0.267671,-0.690881,1.07161,-0.37672,-1.59395,-0.84877,-1.48822,-1.39876,-0.137612,0.991339,
				});*/										// recieved word


				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = (c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}

				//cout << "frame_finished = " << frame_finished << endl;


				//Vector<int, n> nat;
				//Vector_ext::natual(nat);
				//cout << "nat" << nat;
				//cout << "r";
				//r.print();


				osd.solve(r);				// estimated codeword
				if (osd.c_hat != c) {
					frame_error++;
					//cout << "decode error r";
					//r.print();
				}
				TEP_total += osd.TEP_num;
				//ET_frames_total += (osd.is_early_termination == true);

				//osd_dual.solve(r);				// estimated codeword
				//if (osd_dual.c_hat != c) {
				//	frame_error++;
				//	cout << "decode error r";
				//	r.print();
				//	osd.solve(r);				// estimated codeword
				//	if (osd.c_hat == c) {					
				//		cout << "decode should be correct" << endl;
				//	}
				//}

				frame_finished++;
				if (frame_error >= frame_error_max && frame_finished >= frame_min) {
					break;
				}
			}

			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC;
			//cout << "time_consume = " << time_consume << " s" << endl;

			//printf("SNR_dB= %2.1f ;FER= %0.3E ;Errors= %5d ;Frames= %10d ;GF2_ope= %0.3E ;float_ope= %0.3E ;time(s)= %0.3E\n", \
			//	(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
			//	(double)GF2_auxiliary_storage::operation_number / frame_finished, \
			//	(double)my_double_auxiliary_storage::operation_number / frame_finished, time_consume / frame_finished);
			printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E\n", \
				(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
				(double)GF2_auxiliary_storage::operation_number / frame_finished, \
				(double)my_double_auxiliary_storage::operation_number / frame_finished, \
				(double)TEP_total / frame_finished, \
				(double)ET_frames_total / frame_finished, \
				time_consume / frame_finished);
		}

		return 0;
	}

};

bool test_Match_OSD::if_output_to_file = false;