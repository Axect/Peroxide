Lib used for benchmarking: Criterion
Matrix size: 1000x1000

Running benches/parallel_rayon/matrix_benchmark.rs

ser_matrix_bench        time:   [535.12 µs 544.51 µs 556.68 µs]
Found 11 outliers among 100 measurements (11.00%)
  4 (4.00%) high mild
  7 (7.00%) high severe

par_matrix_bench        time:   [5.0912 ms 5.1431 ms 5.1995 ms]
Found 7 outliers among 100 measurements (7.00%)
  1 (1.00%) low mild
  5 (5.00%) high mild
  1 (1.00%) high severe

ser_py_matrix_bench     time:   [4.3100 ms 4.3309 ms 4.3544 ms]
Found 7 outliers among 100 measurements (7.00%)
  2 (2.00%) high mild
  5 (5.00%) high severe

par_py_matrix_bench     time:   [11.667 ms 11.789 ms 11.920 ms]
Found 10 outliers among 100 measurements (10.00%)
  6 (6.00%) high mild
  4 (4.00%) high severe

ser_matrix_change_shape_bench
                        time:   [7.3630 ms 7.4075 ms 7.4608 ms]
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) high mild
  4 (4.00%) high severe

par_matrix_change_shape_bench
                        time:   [10.276 ms 10.385 ms 10.499 ms]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe

ser_matrix_extract_row_bench
                        time:   [613.39 µs 622.44 µs 633.72 µs]
Found 7 outliers among 100 measurements (7.00%)
  7 (7.00%) high severe

par_matrix_extract_row_bench
                        time:   [5.4321 ms 5.4851 ms 5.5399 ms]
Found 4 outliers among 100 measurements (4.00%)
  4 (4.00%) high mild

ser_matrix_from_index_bench
                        time:   [2.4174 ms 2.4490 ms 2.4851 ms]
Found 14 outliers among 100 measurements (14.00%)
  1 (1.00%) high mild
  13 (13.00%) high severe

par_matrix_from_index_bench
                        time:   [2.3912 ms 2.4090 ms 2.4304 ms]
Found 9 outliers among 100 measurements (9.00%)
  2 (2.00%) high mild
  7 (7.00%) high severe

ser_matrix_to_vec_bench time:   [2.4800 ms 2.5082 ms 2.5423 ms]
Found 10 outliers among 100 measurements (10.00%)
  4 (4.00%) high mild
  6 (6.00%) high severe

par_matrix_to_vec_bench time:   [6.4041 ms 6.4618 ms 6.5250 ms]
Found 6 outliers among 100 measurements (6.00%)
  5 (5.00%) high mild
  1 (1.00%) high severe

ser_matrix_to_diag_bench
                        time:   [2.4335 ms 2.4526 ms 2.4750 ms]
Found 14 outliers among 100 measurements (14.00%)
  6 (6.00%) high mild
  8 (8.00%) high severe

par_matrix_to_diag_bench
                        time:   [13.514 ms 13.684 ms 13.868 ms]
Found 10 outliers among 100 measurements (10.00%)
  7 (7.00%) high mild
  3 (3.00%) high severe

Benchmarking ser_matrix_submat_bench: Warming up for 3.0000 s
Warning: Unable to complete 100 samples in 5.0s. You may wish to increase target time to 8.3s, enable flat sampling, or reduce sample count to 50.
ser_matrix_submat_bench time:   [1.6077 ms 1.6243 ms 1.6451 ms]
Found 16 outliers among 100 measurements (16.00%)
  3 (3.00%) high mild
  13 (13.00%) high severe

par_matrix_submat_bench time:   [10.611 ms 10.761 ms 10.942 ms]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe

ser_matrix_add_vec_bench
                        time:   [7.3077 ms 7.3485 ms 7.3946 ms]
Found 12 outliers among 100 measurements (12.00%)
  2 (2.00%) high mild
  10 (10.00%) high severe

par_matrix_add_vec_bench
                        time:   [11.331 ms 11.480 ms 11.636 ms]
Found 2 outliers among 100 measurements (2.00%)
  2 (2.00%) high mild

ser_matrix_norm_bench   time:   [5.1600 ms 5.1864 ms 5.2165 ms]
Found 7 outliers among 100 measurements (7.00%)
  1 (1.00%) high mild
  6 (6.00%) high severe

par_matrix_norm_bench   time:   [2.6565 ms 2.6810 ms 2.7091 ms]
Found 5 outliers among 100 measurements (5.00%)
  2 (2.00%) high mild
  3 (3.00%) high severe

Benchmarking ser_matrix_norm_bench #2: Warming up for 3.0000 s
Warning: Unable to complete 100 samples in 5.0s. You may wish to increase target time to 8.9s, enable flat sampling, or reduce sample count to 50.
ser_matrix_norm_bench #2
                        time:   [1.7262 ms 1.7391 ms 1.7541 ms]
Found 15 outliers among 100 measurements (15.00%)
  10 (10.00%) high mild
  5 (5.00%) high severe

par_matrix_norm_bench #2
                        time:   [6.7071 ms 6.7883 ms 6.8703 ms]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high mild

ser_matrix_norm_bench #3
                        time:   [9.7582 ms 9.9006 ms 10.057 ms]
Found 12 outliers among 100 measurements (12.00%)
  5 (5.00%) high mild
  7 (7.00%) high severe

par_matrix_norm_bench #3
                        time:   [9.3004 ms 9.4088 ms 9.5239 ms]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high mild

ser_matrix_inner_prod_bench
                        time:   [5.2730 ms 5.3590 ms 5.4583 ms]
Found 14 outliers among 100 measurements (14.00%)
  3 (3.00%) high mild
  11 (11.00%) high severe

par_matrix_inner_prod_bench
                        time:   [5.0987 ms 5.1644 ms 5.2402 ms]
Found 7 outliers among 100 measurements (7.00%)
  3 (3.00%) high mild
  4 (4.00%) high severe

ser_matrix_hadamard_bench
                        time:   [5.6521 ms 5.6870 ms 5.7262 ms]
Found 12 outliers among 100 measurements (12.00%)
  3 (3.00%) high mild
  9 (9.00%) high severe

par_matrix_hadamard_bench
                        time:   [14.155 ms 14.335 ms 14.527 ms]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe

ser_matrix_take_row_bench
                        time:   [3.7894 ms 3.8234 ms 3.8613 ms]
Found 15 outliers among 100 measurements (15.00%)
  7 (7.00%) high mild
  8 (8.00%) high severe

par_matrix_take_row_bench
                        time:   [8.4008 ms 8.5171 ms 8.6523 ms]
Found 9 outliers among 100 measurements (9.00%)
  6 (6.00%) high mild
  3 (3.00%) high severe

ser_matrix_fpmap_bench  time:   [3.2526 ms 3.2739 ms 3.2977 ms]
Found 12 outliers among 100 measurements (12.00%)
  2 (2.00%) high mild
  10 (10.00%) high severe

par_matrix_fpmap_bench  time:   [10.604 ms 10.765 ms 10.937 ms]
Found 11 outliers among 100 measurements (11.00%)
  8 (8.00%) high mild
  3 (3.00%) high severe

ser_matrix_reduce_bench time:   [2.6748 ms 2.6964 ms 2.7201 ms]
Found 9 outliers among 100 measurements (9.00%)
  6 (6.00%) high mild
  3 (3.00%) high severe

par_matrix_reduce_bench time:   [6.2453 ms 6.3198 ms 6.4034 ms]
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
