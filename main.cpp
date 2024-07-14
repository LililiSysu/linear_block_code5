/*****************************************************************//**
 * \file   main.cpp
 * \brief  implementation of linear block code, starting from (7,4) code
 * 
 * \author lilili
 * \date   September 2022
 *********************************************************************/

//#define RUN_MSG
#define DEBUG
#define use_my_double
#define count_operation_number
//#define count_add_mul_number
//#define count_compare_operation_number
//#define Complex_output_supress_approx_0
//#define need_trace_table

#include"test/test_include_all.h"

int main_orig()
{
	cout << boolalpha;

	//test_OSD_Chase_orig::matrix_identity_test();
	//test_OSD_Chase_orig::matrix_mul();
	//test_OSD_Chase_orig::matrix_inv();
	//test_OSD_Chase_orig::matrix_end_max();
	//test_OSD_Chase_orig::matrix_low_triangle();
	//test_OSD_Chase_orig::matrix_up_triangle();
	//test_OSD_Chase_orig::matrix_permute();
	//test_OSD_Chase_orig::bch_parity_matrix();
	//test_OSD_Chase_orig::GF_pack();
	//test_OSD_Chase_orig::bch_info();
	//test_OSD_Chase_orig::ebch_info();
	//test_OSD_Chase_orig::ebch_parity_matrix();

	//test_OSD_Chase_orig::run_operation();
	//test_OSD_Chase_orig::run_FER_auto();
	//test_OSD_Chase_orig::run_FER_specified_SNR();
	//test_OSD_Chase_orig::run_operation_all_in_one();
	//test_OSD_Chase_orig::compute_OSD_dependent_col();	// need to change OSD function

	//test_GS_Viterbi_newOSD::plain_old_data();	
	//test_GS_Viterbi_newOSD::matrix_push_back_test();
	//test_GS_Viterbi_newOSD::rank_test();
	//test_GS_Viterbi_newOSD::polynomial_test();
	//test_GS_Viterbi_newOSD::GF_test();
//	test_GS_Viterbi_newOSD::BCH_test();
	//test_GS_Viterbi_newOSD::RS_test();
	//test_GS_Viterbi_newOSD::RS_test2();
	//test_GS_Viterbi_newOSD::RS_test3();
	//test_GS_Viterbi_newOSD::RS_test4();
	//test_GS_Viterbi_newOSD::RS_test5();
	//test_GS_Viterbi_newOSD::RS_test6();
	//test_GS_Viterbi_newOSD::RS_test7();
	//test_GS_Viterbi_newOSD::polynomial_v2_test();
	//test_GS_Viterbi_newOSD::koetter_test();
	//test_GS_Viterbi_newOSD::Q_next_test();
	//test_GS_Viterbi_newOSD::RR_root_finding_test();
	//test_GS_Viterbi_newOSD::n_choose_k_test();
	//test_GS_Viterbi_newOSD::num_test();
	//test_GS_Viterbi_newOSD::mem_test();
	//test_GS_Viterbi_newOSD::GF2e_to_GF_test();
	//test_GS_Viterbi_newOSD::RS_simulate_test();
	//test_GS_Viterbi_newOSD::mod_2_special_test();
	
	//test_GS_Viterbi_newOSD::polynomial_map_test();
	//test_GS_Viterbi_newOSD::polynomial_umap_test();
	//test_GS_Viterbi_newOSD::Lagrange_test();
	//test_GS_Viterbi_newOSD::BCH_test2();
	//test_GS_Viterbi_newOSD::BCH_test3();
	//test_GS_Viterbi_newOSD::BCH_test4();
	//test_GS_Viterbi_newOSD::BCH_test5();
	//test_GS_Viterbi_newOSD::unordered_map_testing();
	//test_GS_Viterbi_newOSD::BCH_non_Gaussian_elimination();
	//test_GS_Viterbi_newOSD::repeatedly_test_BCH_non_Gaussian_elimination();
	//test_GS_Viterbi_newOSD::BCH_non_Gaussian_elimination_simulation();
	//test_GS_Viterbi_newOSD::BCH_non_Gaussian_elimination_simulation2();
	//test_GS_Viterbi_newOSD::BCH_non_Gaussian_elimination_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::viterbi3_test_any_PM();
	//test_GS_Viterbi_newOSD::erase_rows_test();
	//test_GS_Viterbi_newOSD::BCH_test_viterbi3_repeatedly_time();
	//test_GS_Viterbi_newOSD::Sorted_vector_find_test();
	//test_GS_Viterbi_newOSD::OSD_new_test();
	//test_GS_Viterbi_newOSD::OSD_encode_test();
	//test_GS_Viterbi_newOSD::OSD_new_simulation();
	//test_GS_Viterbi_newOSD::viterbi_preprocessing_test();
	//test_GS_Viterbi_newOSD::viterbi_preprocessing_test2(INT_MAX);
	//test_GS_Viterbi_newOSD::viterbi_preprocessing_test_repeatedly();
	//test_GS_Viterbi_newOSD::viterbi_preprocessing_test_orig();
	//test_GS_Viterbi_newOSD::BCH_generate_test();
	//test_GS_Viterbi_newOSD::viterbi_preprocessing_row_transform_test();
	//test_GS_Viterbi_newOSD::viterbi_preprocessing_row_transform_test2();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test2();
	//test_GS_Viterbi_newOSD::Matrix_flex_col_test();
	//test_GS_Viterbi_newOSD::find_opt_PM_test();
	//test_GS_Viterbi_newOSD::find_opt_PM_test2();
	//test_GS_Viterbi_newOSD::find_opt_PM_test3();
	//test_GS_Viterbi_newOSD::find_opt_PM_test4();
	//test_GS_Viterbi_newOSD::find_opt_PM_test5();
	//test_GS_Viterbi_newOSD::find_opt_PM_test6();
	//test_GS_Viterbi_newOSD::find_opt_PM_test7();
	//test_GS_Viterbi_newOSD::find_opt_PM_test8();
	//test_GS_Viterbi_newOSD::find_opt_PM_test9();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test3();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test4();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test5();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test6();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_test7();
	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_simulation();
	//test_GS_Viterbi_newOSD::OSD_r_simulation_old();
	//test_GS_Viterbi_newOSD::OSD_r_simulation_multi_SNR_old();
	//test_GS_Viterbi_newOSD::LC_OSD_r_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::Chase2_simulation_old();
	//test_GS_Viterbi_newOSD::my_float_test();
	//test_GS_Viterbi_newOSD::rand_test();
	//test_GS_Viterbi_newOSD::sort_test();
	//test_GS_Viterbi_newOSD::generate_systematic_generator_any_pos_test();
	//test_GS_Viterbi_newOSD::viterbi_demonstrate();
	//test_GS_Viterbi_newOSD::non_0_pos_test();
	//test_GS_Viterbi_newOSD::LL_OSD_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::LL_OSD_hybrid_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::OSD_r_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::Chase2_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::GMD_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::Hybrid_Chase2_OSD_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::Hamming_code_example();
	//test_GS_Viterbi_newOSD::list_viterbi_zero_columns();
	//test_GS_Viterbi_newOSD::BCH_BM_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::RS_BM_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::RS_BM_test();

	//test_GS_Viterbi_newOSD::RS_7_3_code_test();
	//test_GS_Viterbi_newOSD::RS_7_5_code_test();
	//test_GS_Viterbi_newOSD::RS_7_5_dual_code_test();
	//test_GS_Viterbi_newOSD::linear_block_code_6_3_3_test();
	//test_GS_Viterbi_newOSD::Hamming_code_7_4_3_test();
	//test_GS_Viterbi_newOSD::GF2t_test();
	//test_GS_Viterbi_newOSD::GF2e_to_bits_row_extention();
	//test_GS_Viterbi_newOSD::GF2e_trace_test();
	//test_GS_Viterbi_newOSD::BCH_dual_test();
	//test_GS_Viterbi_newOSD::GF2e_trace_test_only();

	//test_Data_sturcture::heap_test();
	//test_Data_sturcture::vector_test();
	//test_Data_sturcture::flip_TEP_diff_test();
	//test_Data_sturcture::data_ope_test();
	//test_Data_sturcture::choose_order_test();
	//test_Data_sturcture::Euclidean_distance_test();
	//test_Data_sturcture::P2_matrix_test();
	//test_Data_sturcture::RS_Generator_Lagrange_interpolation_test();
	//test_Data_sturcture::matrix_mul_transpose_test();
	//test_Data_sturcture::my_compute_test();
	//test_Data_sturcture::BCH_generate_test();
	//test_Data_sturcture::nnsBCH_generate_test();
	//test_Data_sturcture::npBCH_generate_test();

	//test_Data_sturcture::TEP_jiabao_test();
	//test_Data_sturcture::my_vector_test();
	//test_Data_sturcture::List_remove_only_test();
	//test_Data_sturcture::for_loop_test();

	//test_Data_sturcture::triple_operator_test();
	//test_Data_sturcture::Vec_s_test();
	//test_Data_sturcture::Heap_max_s_test();
	//test_Data_sturcture::sList_s_test();
	//test_Data_sturcture::dList_s_test();		
	//test_Data_sturcture::dList_s_test2();
	
#ifdef need_trace_table
	//test_bin::generate_Generator_of_dual_RS_test();
	//test_bin::count_probability_of_generating_zero_vec();
	//test_bin::trace_generate_linear_combination_test();
#endif

	//test_NewGE_in_OSD::GF2e_inv_test();
	//test_NewGE_in_OSD::generating_all_n_choose_k_pattern_test();
	//test_NewGE_in_OSD::new_generation_of_G_residual_test();
	//test_NewGE_in_OSD::binary_tree_testing();

	//test_NewGE_in_OSD::new_gramma_test();
	//test_NewGE_in_OSD::compare_generator_method();
	//test_NewGE_in_OSD::PreStored_Matrix_test();
	//test_NewGE_in_OSD::PreStored_Matrix_test_multi_trial();	
	//test_NewGE_in_OSD::PreStored_Matrix_test_time();
	//test_NewGE_in_OSD::PreStored_Matrix_red_test_time();
	//test_NewGE_in_OSD::GE_test_time();
	//test_NewGE_in_OSD::GE_left_identity_4_GF2_test();
	//test_NewGE_in_OSD::cycle_code_generator_test();
	//test_NewGE_in_OSD::density_of_BCH_generator_matrix();
	//test_NewGE_in_OSD::density_of_BCH_generator_matrix_2();
	//test_NewGE_in_OSD::check_out_eBCH_test();
	
	//test_NewGE_in_OSD::num_test();
	//test_NewGE_in_OSD::number_sequence_test();
	//test_NewGE_in_OSD::GE_worst_case_ope_estimation_test();
	//test_NewGE_in_OSD::BU_worst_case_ope_estimation_test();
	//test_NewGE_in_OSD::PreStored_Matrix_red_recursive_partition_test_ini();
	//test_NewGE_in_OSD::PreStored_Matrix_red_recursive_partition_test_multi_trial();
	//test_NewGE_in_OSD::n_ary_incremental_scan_test();
	//test_NewGE_in_OSD::n_randomly_divide_into_k_and_res_test();	
	
	//test_NewGE_in_OSD::PreStored_Matrix_red_recursive_partition_test_time();
	//test_NewGE_in_OSD::PreStored_Matrix_red_loose_bin_test_ini();
	//test_NewGE_in_OSD::PreStored_Matrix_red_loose_bin_test_multi_trial();	
	
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_test_ini();
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_test_multi_trial();		
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_fill_test_ini();
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_fill_test_multi_trial();	
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_cycle_test_ini();
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_cycle_test_multi_trial();	
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_extend_cycle_test_ini();
	//test_NewGE_in_OSD::PreStored_Matrix_red_equal_dist_extend_cycle_test_multi_trial();
	//test_NewGE_in_OSD::PreStored_Matrix_red_cycle_ini_test();
	//test_NewGE_in_OSD::PreStored_Matrix_red_cycle_test_multi_trial();
	//test_NewGE_in_OSD::GE_IBU_IBUc_test_all();
	//test_NewGE_in_OSD::GJE_test();
	//test_NewGE_in_OSD::GE_IBU_IBUc_test_all_2();

//	test_NewGE_in_OSD::OSD_v2_simulation_multi_SNR();
	//test_NewGE_in_OSD::OSD_v2_fixed_received_vector();
	//test_NewGE_in_OSD::OSD_PSM_aided_GE_simulation_multi_SNR();
	//test_NewGE_in_OSD::OSD_PSM_aided_GE_4_cyc_simulation_multi_SNR();
	//test_NewGE_in_OSD::OSD_PSM_aided_GE_4_ext_cyc_simulation_multi_SNR();
	//test_NewGE_in_OSD::OSD_v2_and_PSM_aided_double_ope_test();
	//test_NewGE_in_OSD::OSD_v2_simulation_complexity_distribution(5.8);

	//test_IBU::IBU_init_test();
	//test_IBU::IBU_test();
	//test_IBU::IBUc_init_test();
	//test_IBU::IBUc_test();

	//test_Min_weight::BCH_dual_test();
	
	return 0;
}

int main_MA() {

	// it seems cannot work for matrix polynomial now, problem to be solved

	//Matrix_analysis::polybomial_determinant_test();
	//Matrix_analysis::polynomial_matrix_test();
	//Matrix_analysis::polynomial_matrix_test_simple();
	//Matrix_analysis::add_orthogonal_complementary_test();
	//Matrix_analysis::ritz_test();
	//Matrix_analysis::ritz_inner_test();
	Matrix_analysis::edobtm_test();

	//Matrix_analysis::Hermitian_test();
	//Matrix_analysis::givens_Hessenberg_test();
	//Matrix_analysis::givens_Tridiagonal_test();
	//Matrix_analysis::ritz_inner_test2();
	//Matrix_analysis::set_cols_test();
	//Matrix_analysis::ritz_outer_test();
	//Matrix_analysis::ritz_test2();
	//Matrix_analysis::svd_val_test();
	//Matrix_analysis::svd_test();

	return 0;
}

int main_arg(int argc, char* argv[])
{
	cout << boolalpha;

	test_LL_OSD_for_Jiabao::LL_OSD_genius_aid_simulation_multi_SNR(atoi(argv[1]), atoi(argv[2]));
	//test_LL_OSD_for_Jiabao::LL_OSD_genius_aid_simulation_multi_SNR(3,8);
	//test_LL_OSD_for_Jiabao::TEP_jiabao_simple_look_test();
	//test_LL_OSD_for_Jiabao::TEP_generate_segmentation_test();
	//test_LL_OSD_for_Jiabao::generate_order_pattern_test();

	return 0;
}

int main/*_new*/(){

	//test_Match_OSD::Match_OSD_initialize();
	//test_Match_OSD::Match_OSD_FER();

	//test_GS_Viterbi_newOSD::viterbi_optimized_PM_simulation();
	//test_GS_Viterbi_newOSD::LC_OSD_r_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::LL_OSD_simulation_multi_SNR();
	//test_GS_Viterbi_newOSD::generate_systematic_generator_any_pos_test();
	//test_GS_Viterbi_newOSD::BCH_non_Gaussian_elimination_simulation_multi_SNR();		// LL_OSD_Viterbi
	//test_GS_Viterbi_newOSD::LL_OSD_hybrid_simulation_multi_SNR();						// LL_OSD_Hybrid_flip_Viterbi
	//test_GS_Viterbi_newOSD::BCH_BM_simulation_multi_SNR();							// BM
	//test_GS_Viterbi_newOSD::BCH_info();
	//test_GS_Viterbi_newOSD::BCH_info2();
	//test_GS_Viterbi_newOSD::BCH_info3();

	//test_BCH_permutation::eBCH_affine_permutation();
	//test_BCH_permutation::eBCH_affine_permutation2();
	//test_BCH_permutation::BCH_weight();
	//test_BCH_permutation::BCH_dual_weight();
	//test_BCH_permutation::BCH_punctured_wieght();
	//test_BCH_permutation::BCH_punctured_wieght_multi();
	//test_BCH_permutation::BCH_punctured_wieght_flip_multi();
	//test_BCH_permutation::nnsBCH_weight();
	
	//test_IBU_OSD::OSD_v4_FER();
	//test_IBU_OSD::IBUc_test();
	//test_IBU_OSD::GE_IBU_IBUc();
	//test_IBU_OSD::GJE_right();
	//test_IBU_OSD::GE_IBU_IBUc_v2();
	test_IBU_OSD::OSD_v4_FER_new_test();
	
	//test_LC_OSD::LC_OSD_FER();
	//test_LL_OSD::LL_OSD_FER_old();
	//test_LL_OSD::LL_OSD_FER();
	//test_LL_OSD::LL_OSD_Viterbi_FER();

	//test_GE_skip_OSD::FER();
	//test_GE_skip_OSD::test_qfunc();
	//test_GE_skip_OSD::test_GE_v2();

	//test_GRAND::FER();

	return 0;
}

//
//int main()
//{
//	Matrix<int> A(3, 2, {
//		0,1,
//		0,0,
//		1,1,
//	});
//
//	cout << "A = " << A;
//	A = A.erase_all_0_rows();
//	cout << "A = " << A;
//	return 0;
//}
