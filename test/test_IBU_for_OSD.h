/*****************************************************************//**
 * \file   IBU_for_OSD.h
 * \brief  From linear_block_code_4
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

class test_IBU_for_OSD_v4 {
public:
	static bool if_output_to_file;
	static void init_frame_SNR() {

		bool run_performance = true;
		if_output_to_file = true;

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
			SNR_dB_start = 6;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
	}

	// OSD test FER
	static int OSD_v4_FER()
	{
		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 7;
		const int t = 10;
		// OSD parameter
		int order = 4;

		GF2e<m>::init();
		eBCH<m, t> bch;
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

		Matrix<GF2> u(1,k);											// message
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
	

		Matrix<int> nat(1, n, 'N');
		Matrix<GF2> G_sys = G;
		G_sys.GJE_4_GF2(nat);

		Matrix<GF2> H_sys = H;
		H_sys.GJE_4_GF2(nat);
		int s = m;
		//cout << "H_sys" << H_sys;
	
		OSD_v4 osd(G, d, order, true);
		//OSD_v4 osd(H, d, order, false);
		//OSD_IBUc osd(G_sys, d, s, order, false);		//osd.ibuc.reset_s(1);
		//OSD_IBUc osd(G_sys, d, s, order, true);		//osd.ibuc.reset_s(s);
		//OSD_v4_dual osd(H, d, order);
		//OSD_IBUc_dual osd(H_sys, d, s, order, false);		//osd.ibuc.reset_s(1);
		//OSD_IBUc_dual osd(H_sys, d, s, order, true);		//osd.ibuc.reset_s(s);
		

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			if (osd.type == 1)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD(%d).txt", n, k, d, order);		// it semms harder and no state decrease.
			else if(osd.type == 2)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_IBU(%d).txt", n, k, d, order);
			else if (osd.type == 3)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_IBUc(%d).txt", n, k, d, order);
			else if (osd.type == 4)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_dual(%d).txt", n, k, d, order);
			else if (osd.type == 5)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_IBU_dual(%d).txt", n, k, d, order);
			else if (osd.type == 6)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_IBUc_dual(%d).txt", n, k, d, order);
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

			start = clock();
			while (frame_finished < frame_max) {

				/*Matrix<my_double, 1, n> r({
					-0.929966,-0.777951,-0.852718,1.35274,0.0624572,-0.267671,-0.690881,1.07161,-0.37672,-1.59395,-0.84877,-1.48822,-1.39876,-0.137612,0.991339,
				});*/										// recieved word


				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = (c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}


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
				ET_frames_total += (osd.is_early_termination == true);

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

	// IBUc test
	static int IBUc_test()
	{
		const int n = 15;
		const int k = 7;
		const int d = 5;
		const int s = 4;
		Matrix<GF2> Gs(k, n,{
			1,0,0,0,0,0,0,1,0,0,0,1,0,1,1,
			0,1,0,0,0,0,0,1,1,0,0,1,1,1,0,
			0,0,1,0,0,0,0,0,1,1,0,0,1,1,1,
			0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,
			0,0,0,0,1,0,0,0,1,0,1,1,1,0,0,
			0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,
			0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,
			});
		Matrix<GF2> H(n - k, n, {
			1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,
			0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,
			0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,
			0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,
			0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,
			0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,
			0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,
			0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,
			});
		cout << "Gs" << Gs;
		cout << "H" << H;
		Matrix<GF2> H_T = H.Transpose();
		cout << "Gs * H_T" << Gs * H_T;

		cout << "-------------------------" << endl << endl;

		Matrix<GF2> c(1, n, { 1,   1,   0,   0,   0,   1,   1,   0,   1,   1,   1,   1,   1,   0,   0, });

		Matrix<my_double> r(1, n, {
			-1.71704,	-0.0355986,	0.953004,	1.03202,	1.27793,	-0.130623,	-0.81576,
			1.57263,	-1.14507,	-0.151133,	-0.600191,	-0.54,	0.550899,	-0.165649,	1.15493,
			});										// recieved word

		cout << "r" << r;
		Matrix<int> nat(1, n, 'N');
		//Vector_ext::natual(nat);
		cout << "nat" << nat;

		// IBU method
		IBUc ibuc(n, k, s);
		ibuc.solve(Gs, r);

		cout << "-------------------------" << endl << endl;
		cout << "-------------------------" << endl << endl;

		cout << "ibu.Gs_tilde" << ibuc.Gs;
		cout << "ibu.permute_record_BU" << ibuc.r_abs_Gs_record;

		Matrix<my_double> ibu_r_tilde = r;
		ibu_r_tilde.permute(ibuc.r_abs_Gs_record);
		cout << "ibuc.r_tilde" << ibu_r_tilde;
		cout << "ibuc.r_abs_tilde_record" << ibuc.r_abs_Gs_record;
		Matrix<GF2> ibu_Gs_back = ibuc.Gs;
		ibu_Gs_back.permute_col_back(ibuc.r_abs_Gs_record);
		cout << "ibu_Gs_back * H_T" << ibu_Gs_back * H_T;

		cout << "-------------------------" << endl << endl;

		// GE method
		Matrix<my_double> r_abs(1, n);
		for (int i = 0; i < n; ++i) {
			r_abs[i] = my::abs(r[i]);
		}
		Matrix<int> permute_record(1, n, 'N');
		Matrix<my_double> r_abs_bar = r_abs;
		r_abs_bar.quick_sort_recur_gt_with_ind(0,n - 1, permute_record);
		Matrix<GF2> Gs_tilde = Gs;
		Gs_tilde.permute_col(permute_record);
		Gs_tilde.GJE_4_GF2(permute_record);

		cout << "Gs_tilde" << Gs_tilde;
		cout << "permute_record" << permute_record;
		Matrix<my_double> r_tilde = r;
		r_tilde.permute(permute_record);
		cout << "r_tilde" << r_tilde;
		Matrix<GF2> Gs_back = Gs_tilde;
		Gs_back.permute_col_back(permute_record);
		cout << "Gs_back * H_T" << Gs_back * H_T;
		return 0;
	}

	// GE, IBU, IBUc test operation number
	static int GE_IBU_IBUc()
	{
		// BCH code parameter
		const int m = 5;
		const int t = 3;
		// OSD parameter
		int order = 2;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int n = G.col();
		int k = G.row();
		int d = 2 * t + 1;

		cout << "G = " << endl;
		G.print();
		cout << "H = " << endl;
		H.print();

		cout << "-------------------------" << endl << endl;

		Matrix<GF2> u(1, k);											// message
		for (int i = 0; i < k; ++i) {
			u[i] = my::rand_01();
		}
		Matrix<GF2> c = u * G;									// transmitted codeword
		cout << "c = " << endl;
		c.print();
		cout << "-------------------------" << endl << endl;


		Matrix<my_double> SNR_dB_vec(SNR_dB_start, SNR_dB_step, SNR_dB_end, 'd');
		//Vector_ext::equally_distant(SNR_dB_start, SNR_dB_end, SNR_dB_step, SNR_dB_vec);// 	SNR_dB_vec[0] = 6;
		//cout << "SNR_dB_vec" << SNR_dB_vec;

		cout << "BCH(" << n << ", " << k << ", " << d << ")" << endl;

		Matrix<int> nat(1, n, 'N');
		Matrix<GF2> G_sys = G;
		G_sys.GJE_4_GF2(nat);

		Matrix<GF2> H_sys = H;
		H_sys.GJE_4_GF2(nat);
		int s = m;
		//cout << "H_sys" << H_sys;

		OSD_v4 osd(G, d, order, true);
		OSD_IBUc osd_IBU(G_sys, d, s, order, false);
		OSD_IBUc osd_IBUc(G_sys, d, s, order, true);		//osd.ibuc.reset_s(s);

		OSD_v4_dual osd_dual(H, d, order);
		OSD_IBUc_dual osd_IBU_dual(H_sys, d, s, order, false);
		OSD_IBUc_dual osd_IBUc_dual(H_sys, d, s, order, true);		//osd.ibuc.reset_s(s);

		bool if_output_to_file = false;

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };

			sprintf_s(file_name, 155, "(%d,%d,%d)BCH_GE_IBU_IBUc_complexity.txt", n, k, d);		// it semms harder and no state decrease.
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		clock_t start, end;

		my::set_seed_adv(1);

		for (int i = 0; i < SNR_dB_vec.size(); ++i) {

			double SNR = (double)pow(10, SNR_dB_vec[i] / 10);
			double sigma = sqrt(n / (2 * k * SNR));
			int frame_finished = 0;

			unsigned long long GE_operations = 0;
			unsigned long long IBU_operations = 0;
			unsigned long long IBUc_operations = 0;

			unsigned long long GE_dual_operations = 0;
			unsigned long long IBU_dual_operations = 0;
			unsigned long long IBUc_dual_operations = 0;

			//my_double_auxiliary_storage::operation_number = 0;

			start = clock();
			while (frame_finished < frame_min) {

				/*Matrix<my_double, 1, n> r({
					-0.929966,-0.777951,-0.852718,1.35274,0.0624572,-0.267671,-0.690881,1.07161,-0.37672,-1.59395,-0.84877,-1.48822,-1.39876,-0.137612,0.991339,
				});*/										// recieved word


				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = (c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}

				GF2_auxiliary_storage::operation_number = 0;
				osd.call_GE_only(r);
				GE_operations += GF2_auxiliary_storage::operation_number;

				GF2_auxiliary_storage::operation_number = 0;
				osd_IBU.call_IBUc_only(r);
				IBU_operations += GF2_auxiliary_storage::operation_number;

				GF2_auxiliary_storage::operation_number = 0;
				osd_IBUc.call_IBUc_only(r);
				IBUc_operations += GF2_auxiliary_storage::operation_number;

				GF2_auxiliary_storage::operation_number = 0;
				osd_dual.call_GE_only(r);
				GE_dual_operations += GF2_auxiliary_storage::operation_number;

				GF2_auxiliary_storage::operation_number = 0;
				osd_IBU_dual.call_IBUc_only(r);
				IBU_dual_operations += GF2_auxiliary_storage::operation_number;

				GF2_auxiliary_storage::operation_number = 0;
				osd_IBUc_dual.call_IBUc_only(r);
				IBUc_dual_operations += GF2_auxiliary_storage::operation_number;

				frame_finished++;
			}

			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC;
			//cout << "time_consume = " << time_consume << " s" << endl;

			printf("GE_operations= %0.3E \nIBU_operations= %0.3E \nIBUc_operations= %0.3E \nGE_dual_operations= %0.3E \nIBU_dual_operations= %0.3E \nIBUc_dual_operations= %0.3E \ntime(s)= %0.3E\n", \
				(double)GE_operations / frame_finished, (double)IBU_operations / frame_finished, \
				(double)IBUc_operations / frame_finished, (double)GE_dual_operations / frame_finished, \
				(double)IBU_dual_operations / frame_finished, (double)IBUc_dual_operations / frame_finished, time_consume);
		}
		return 0;
	}

	static int GJE_right() {
		const int m = 6;
		const int t = 4;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int n = G.col();
		int k = G.row();

		cout << "G" << G;
		cout << "H" << H;

		Matrix<GF2> Gs;
		Matrix<GF2> Hs;

		//Matrix<int> perumtation_first(1, n, { 6, 5,10,13,14, 7, 9, 3,12, 2, 1, 4, 0,11, 8});
		Matrix<int> perumtation_first(1, n, 'N');
		Matrix<int> permutation_second(1, n, 'N'); // Vector_ext::natual<n>(permutation_second);

		my::set_seed_adv(1);
		for (int p = 0; p < 10000; ++p) {					

			perumtation_first.permute_rand();
			//cout << "perumtation_first" << endl;
			//perumtation_first.print();

			Gs = G;
			Gs.permute_col(perumtation_first);
			Gs.GJE_4_GF2(permutation_second);
			//cout << "permutation_second" << permutation_second;
			//cout << "Gs" << Gs;

			Hs = H;
			Hs.permute_col(perumtation_first);
			Hs.GJE_4_GF2_right(permutation_second);
			Matrix<GF2> Gs_new(k, n, 'i');

			// set Gs = [I_k | P^T]
			for (int i = 0; i < k; ++i) {
				for (int j = k; j < n; ++j) {
					Gs_new(i, j) = Hs(j - k, i);
				}
			}
			//cout << "Gs_new" << Gs_new;

			if ((Gs - Gs_new).isZero() == false) {
				cout << "perumtation_first" << perumtation_first;
				cout << "Gs - Gs_new" << Gs - Gs_new;		// no print, all success
			}
		}

		return 0;
	}
};

bool test_IBU_for_OSD_v4::if_output_to_file = false;
