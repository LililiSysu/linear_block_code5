/*****************************************************************//**
 * \file   test_GE_skip_OSD.h
 * \brief  test class for 
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
extern my_double SNR_dB_end;		// can be omitted

class test_GE_skip_OSD {
public:
	// The following statement is allowed
	//int a = 5;
	//int b[3][2] = { {3,1},{4,8},{10,0} };
	//static const bool if_output_to_file = false;

	static bool if_output_to_file;

	static void init_frame_SNR() {

		bool run_performance = true;
		if_output_to_file = true;

		if (run_performance) {
			frame_max = 50000000;			// 50000000
			frame_min = 10000;			// 10000
			frame_error_max = 200;		// 200
			SNR_dB_start = 1;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
		else {
			frame_max = 100000;
			frame_min = 100000;
			frame_error_max = 1000000;
			SNR_dB_start = 4;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
	}

	// OSD test FER
	static int FER()
	{
		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 7;
		const int t = 4;
		// OSD parameter
		int order = 1;

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


		//Pc_SGM_OSD osd(G, H, d, order);
		Adaptive_GE_OSD osd(G, H, d, order);

		osd.ET_type = 1;		// ET_type to be updated in IBU and IBUc of OSD

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			if (osd.type == 31)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_Pc_SGM_OSD(%d).txt", n, k, d, order);
			else if (osd.type == 34)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_Pc_SGM_OSD_dual(%d).txt", n, k, d, order);
			else if (osd.type == 41)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_Adaptive_GE_OSD(%d).txt", n, k, d, order);
			else if (osd.type == 44)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_Adaptive_GE_OSD_dual(%d).txt", n, k, d, order);
			else
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_unknown(%d).txt", n, k, d, order);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}
		else {
			printf("\n");
		}
		//printf("SNR(dB)   FER          Errors    Frames        GF2_ope      float_ope    time(s)\n");
		printf("SNR(dB)   FER          Errors    Frames        GF2_ope      float_ope    ave_list     ET_rate      GE_skip_rate     time(s)\n");
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
			int GE_reduction_frame_total = 0;

			osd.sigma = sigma;						// ET, to be extended in 'OSD_v4_dual'



			start = clock();
			while (frame_finished < frame_max) {

				/*Matrix<my_double, 1, n> r({
					-0.929966,-0.777951,-0.852718,1.35274,0.0624572,-0.267671,-0.690881,1.07161,-0.37672,-1.59395,-0.84877,-1.48822,-1.39876,-0.137612,0.991339,
				});*/										// recieved word

				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = ((int)c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel, supress operation counting
				}


				//Vector<int, n> nat;
				//Vector_ext::natual(nat);
				//cout << "nat" << nat;
				//cout << "r";
				//r.print();

				//osd.call_GE_only(r);
				osd.solve(r);				// estimated codeword
				if (osd.c_hat != c) {
					frame_error++;
					//cout << "decode error r";
					//r.print();
				}
				TEP_total += osd.TEP_num;
				ET_frames_total += (osd.is_early_termination == true);
				GE_reduction_frame_total += (osd.is_GE_skip == true);

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

			//printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E\n", \
			//	(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
			//	(double)GF2_auxiliary_storage::operation_number / frame_finished, \
			//	(double)my_double_auxiliary_storage::operation_number / frame_finished, \
			//	(double)TEP_total / frame_finished, \
			//	(double)ET_frames_total / frame_finished, \
			//	time_consume / frame_finished);

			printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E        %0.3E\n", \
				(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
				(double)GF2_auxiliary_storage::operation_number / frame_finished, \
				(double)my_double_auxiliary_storage::operation_number / frame_finished, \
				(double)TEP_total / frame_finished, \
				(double)ET_frames_total / frame_finished, \
				(double)GE_reduction_frame_total / frame_finished, \
				time_consume / frame_finished);
		}

		return 0;
	}

	static int test_qfunc() {
		double x;
		for (x = -6; x < 6; x += 0.1) {
			cout << "my::qfunc(" << x << ") = " << my::qfunc(x) << endl;
		}
		return 0;
	}

	static int test_GE_v2() {
		// BCH code parameter
		const int m = 6;
		const int t = 3;
		// OSD parameter
		int order = 1;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int n = G.col();
		int k = G.row();
		int d = bch.get_d();

		//cout << "G" << G << endl;

		Matrix<GF2> G_sys_1 = G;
		Matrix<int> permute_record_1(1, n, 'N');
		G_sys_1.GJE_4_GF2_left(permute_record_1);

		//cout << "G_sys_1" << G_sys_1;
		//cout << "permute_record_1" << permute_record_1;

		Matrix<GF2> G_sys_2 = G;
		Matrix<int> permute_record_2(1, n, 'N');
		G_sys_2.GE_left_identity_match_jiabao(permute_record_2);

		//cout << "G_sys_2" << G_sys_2;
		//cout << "permute_record_2" << permute_record_2;

		cout << "(G_sys_2 == G_sys_1) : " << (G_sys_2 == G_sys_1) << endl;
		cout << "(permute_record_2 == permute_record_1) : " << (permute_record_2 == permute_record_1) << endl;

		for (int i = 0; i < 1000; ++i) {
			Matrix<int> rand_permutation(1, n, 'N');
			rand_permutation.permute_rand();


			Matrix<GF2> G_sys_1 = G;
			G_sys_1.permute_col(rand_permutation);
			Matrix<int> permute_record_1(1, n, 'N');
			G_sys_1.GJE_4_GF2_left(permute_record_1);

			Matrix<GF2> G_sys_2 = G;
			G_sys_2.permute_col(rand_permutation);
			Matrix<int> permute_record_2(1, n, 'N');
			G_sys_2.GE_left_identity_match_jiabao(permute_record_2);

			if (G_sys_2 != G_sys_1) {
				cout << "(permute_record_2 == permute_record_1) : " << (permute_record_2 == permute_record_1) << endl;
				cout << "rand_permutation" << rand_permutation;
			}
		}

		{
			//cout << "H" << H << endl;

			Matrix<GF2> H_sys_1 = H;
			Matrix<int> permute_record_1(1, n, 'N');
			H_sys_1.GJE_4_GF2_right(permute_record_1);

			//cout << "H_sys_1" << H_sys_1;
			//cout << "permute_record_1" << permute_record_1;

			Matrix<GF2> H_sys_2 = H;
			Matrix<int> permute_record_2(1, n, 'N');
			H_sys_2.GE_right_identity_match_jiabao(permute_record_2);

			//cout << "H_sys_2" << H_sys_2;
			//cout << "permute_record_2" << permute_record_2;

			cout << "(H_sys_2 == H_sys_1) : " << (H_sys_2 == H_sys_1) << endl;
			cout << "(permute_record_2 == permute_record_1) : " << (permute_record_2 == permute_record_1) << endl;
		}

		for (int i = 0; i < 1000; ++i) {
			Matrix<int> rand_permutation(1, n, 'N');
			rand_permutation.permute_rand();


			Matrix<GF2> H_sys_1 = H;
			H_sys_1.permute_col(rand_permutation);
			Matrix<int> permute_record_1(1, n, 'N');
			H_sys_1.GJE_4_GF2_right(permute_record_1);

			Matrix<GF2> H_sys_2 = H;
			H_sys_2.permute_col(rand_permutation);
			Matrix<int> permute_record_2(1, n, 'N');
			H_sys_2.GE_right_identity_match_jiabao(permute_record_2);

			if (H_sys_2 != H_sys_1) {
				cout << "(permute_record_2 == permute_record_1) : " << (permute_record_2 == permute_record_1) << endl;
				cout << "rand_permutation" << rand_permutation;
			}
		}

		return 0;
	}
};

bool test_GE_skip_OSD::if_output_to_file = false;
