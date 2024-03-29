#ifndef FUNCTION_POINTERS_C
#define FUNCTION_POINTERS_C

#include"../Macro/macro.h"

#include"function_pointers.h"
#include"../G2/g2.h"
#include"../G2/g2_check.h"
#include"../G2/g2_un.h"
#include"../G2/g2_upd.h"
#include"../Su2/su2.h"
#include"../Su2/su2_upd.h"
#include"../SuN/sun.h"
#include"../SuN/sun_upd.h"



void init_function_pointers(void)
  {
  #ifdef SuN_GROUP 
  one  = &one_SuN;
  zero = &zero_SuN;
  equal     = &equal_SuN;
  equal_dag = &equal_dag_SuN;
  plus_equal     = &plus_equal_SuN;
  plus_equal_dag = &plus_equal_dag_SuN;
  minus_equal     = &minus_equal_SuN;
  minus_equal_dag = &minus_equal_dag_SuN;
  lin_comb       = &lin_comb_SuN;
  lin_comb_dag1  = &lin_comb_dag1_SuN;
  lin_comb_dag2  = &lin_comb_dag2_SuN;
  lin_comb_dag12 = &lin_comb_dag12_SuN;
  times_equal_real = &times_equal_real_SuN;
  times_equal     = &times_equal_SuN;
  times_equal_dag = &times_equal_dag_SuN;
  times       = &times_SuN;
  times_dag1  = &times_dag1_SuN;
  times_dag2  = &times_dag2_SuN;
  times_dag12 = &times_dag12_SuN;
  rand_matrix = &rand_matrix_SuN;
  retr = &retr_SuN;
  imtr = &imtr_SuN;
  unitarize = &unitarize_SuN;
  print_on_screen = &print_on_screen_SuN;
  print_on_file   = &print_on_file_SuN;
  print_on_binary_file   = &print_on_binary_file_SuN;
  read_from_file   = &read_from_file_SuN;
  read_from_binary_file   = &read_from_binary_file_SuN;
  read_from_binary_file_swap   = &read_from_binary_file_swap_SuN;
  single_heatbath = &single_heatbath_SuN;
  single_overrelaxation = &single_overrelaxation_SuN;
  cool = &cool_SuN;
  #endif


  #ifdef Su2_GROUP 
  one  = &one_Su2;
  zero = &zero_Su2;
  equal     = &equal_Su2;
  equal_dag = &equal_dag_Su2;
  plus_equal     = &plus_equal_Su2;
  plus_equal_dag = &plus_equal_dag_Su2;
  minus_equal     = &minus_equal_Su2;
  minus_equal_dag = &minus_equal_dag_Su2;
  lin_comb       = &lin_comb_Su2;
  lin_comb_dag1  = &lin_comb_dag1_Su2;
  lin_comb_dag2  = &lin_comb_dag2_Su2;
  lin_comb_dag12 = &lin_comb_dag12_Su2;
  times_equal_real = &times_equal_real_Su2;
  times_equal     = &times_equal_Su2;
  times_equal_dag = &times_equal_dag_Su2;
  times       = &times_Su2;
  times_dag1  = &times_dag1_Su2;
  times_dag2  = &times_dag2_Su2;
  times_dag12 = &times_dag12_Su2;
  rand_matrix = &rand_matrix_Su2;
  retr = &retr_Su2;
  imtr = &imtr_Su2;
  unitarize = &unitarize_Su2;
  print_on_screen = &print_on_screen_Su2;
  print_on_file   = &print_on_file_Su2;
  print_on_binary_file   = &print_on_binary_file_Su2;
  read_from_file   = &read_from_file_Su2;
  read_from_binary_file   = &read_from_binary_file_Su2;
  read_from_binary_file_swap   = &read_from_binary_file_swap_Su2;
  single_heatbath = &single_heatbath_Su2;
  single_overrelaxation = &single_overrelaxation_Su2;
  cool = &cool_Su2;
  #endif
 
  #ifdef G2_GROUP
  one  = &one_G2;
  zero = &zero_G2;
  equal     = &equal_G2;
  equal_dag = &equal_dag_G2;
  plus_equal     = &plus_equal_G2;
  plus_equal_dag = &plus_equal_dag_G2;
  minus_equal     = &minus_equal_G2;
  minus_equal_dag = &minus_equal_dag_G2;
  lin_comb       = &lin_comb_G2;
  lin_comb_dag1  = &lin_comb_dag1_G2;
  lin_comb_dag2  = &lin_comb_dag2_G2;
  lin_comb_dag12 = &lin_comb_dag12_G2;
  times_equal_real = &times_equal_real_G2;
  times_equal     = &times_equal_G2;
  times_equal_dag = &times_equal_dag_G2;
  times       = &times_G2;
  times_dag1  = &times_dag1_G2;
  times_dag2  = &times_dag2_G2;
  times_dag12 = &times_dag12_G2;
  rand_matrix = &rand_matrix_G2;
  retr = &retr_G2;
  imtr = &imtr_G2;
  unitarize = &unitarize_G2;
  print_on_screen = &print_on_screen_G2;
  print_on_file   = &print_on_file_G2;
  print_on_binary_file   = &print_on_binary_file_G2;
  read_from_file   = &read_from_file_G2;
  read_from_binary_file   = &read_from_binary_file_G2;
  read_from_binary_file_swap   = &read_from_binary_file_swap_G2;
  single_heatbath = &single_heatbath_G2;
  single_overrelaxation = &single_overrelaxation_G2;
  cool = &cool_G2;
  #endif
  }
#endif
