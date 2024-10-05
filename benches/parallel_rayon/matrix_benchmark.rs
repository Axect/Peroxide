use criterion::{black_box, criterion_group, criterion_main, Criterion};
use peroxide::{
    fuga::*,
    traits::{
        fp::ParallelFPMatrix,
        math::{ParallelInnerProduct, ParallelMatrixProduct, ParallelNormed, ParallelVector},
    },
};

pub fn par_matrix_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 630.92 µs
    cr.bench_function("ser_matrix_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row)))
    });

    // Result: 1000x1000 matrix: 9.6995 ms
    cr.bench_function("par_matrix_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row)))
    });
}

pub fn par_py_matrix_benchmark(cr: &mut Criterion) {
    let v: Vec<Vec<f64>> = (0..1000)
        .into_iter()
        .map(|i: i32| {
            (0..1000)
                .into_iter()
                .map(|j| 2.0 * (i as f64) * (j as f64))
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>();

    // Result: 1000x1000 matrix: 3.9858 ms
    cr.bench_function("ser_py_matrix_bench", |b| {
        b.iter(|| black_box(py_matrix(v.clone())))
    });

    // Result: 1000x1000 matrix: 24.654 ms
    cr.bench_function("par_py_matrix_bench", |b| {
        b.iter(|| black_box(par_py_matrix(v.clone())))
    });
}

pub fn par_matrix_change_shape_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 7.6453 ms
    cr.bench_function("ser_matrix_change_shape_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).change_shape()))
    });

    // Result: 1000x1000 matrix: 9.8586 ms
    cr.bench_function("par_matrix_change_shape_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_change_shape()))
    });
}

pub fn par_matrix_extract_row_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 571.90 µs
    cr.bench_function("ser_matrix_extract_row_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).row(100)))
    });

    // Result: 1000x1000 matrix: 10.440 ms
    cr.bench_function("par_matrix_extract_row_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_row(100)))
    });
}

pub fn par_matrix_from_index_benchmark(cr: &mut Criterion) {
    let f = |x: usize, y: usize| 2.0 * (x as f64) * (y as f64);
    let size: (usize, usize) = (1000, 1000);

    // Result: 1000x1000 matrix: 2.3662 ms
    cr.bench_function("ser_matrix_from_index_bench", |b| {
        b.iter(|| black_box(Matrix::from_index(f, size)))
    });

    // Result: 1000x1000 matrix:  2.3355 ms
    cr.bench_function("par_matrix_from_index_bench", |b| {
        b.iter(|| black_box(Matrix::from_index(f, size)))
    });
}

pub fn par_matrix_to_vec_and_diag_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 5.0333 ms
    cr.bench_function("ser_matrix_to_vec_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).to_vec()))
    });

    // Result: 1000x1000 matrix: 10.625 ms
    cr.bench_function("par_matrix_to_vec_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_to_vec()))
    });

    // Result: 1000x1000 matrix: 3.0501 ms
    cr.bench_function("ser_matrix_to_diag_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).to_diag()))
    });

    // Result: 1000x1000 matrix: 16.362 ms
    cr.bench_function("par_matrix_to_diag_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_to_diag()))
    });
}

pub fn par_matrix_submat_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 1.5658 ms
    cr.bench_function("ser_matrix_submat_bench", |b| {
        b.iter(|| {
            black_box(matrix(v.clone(), 1000, 1000, Shape::Row).submat((100, 100), (900, 900)))
        })
    });

    // Result: 1000x1000 matrix: 13.000 ms
    cr.bench_function("par_matrix_submat_bench", |b| {
        b.iter(|| {
            black_box(
                par_matrix(v.clone(), 1000, 1000, Shape::Row).par_submat((100, 100), (900, 900)),
            )
        })
    });
}

pub fn par_matrix_add_vec_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    let w: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 3.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 7.1774 ms
    cr.bench_function("ser_matrix_add_vec_bench", |b| {
        b.iter(|| {
            black_box(matrix(v.clone(), 1000, 1000, Shape::Row).add_vec(&matrix(
                w.clone(),
                1000,
                1000,
                Shape::Row,
            )))
        })
    });

    // Result: 1000x1000 matrix: 11.565 ms
    cr.bench_function("par_matrix_add_vec_bench", |b| {
        b.iter(|| {
            black_box(
                matrix(v.clone(), 1000, 1000, Shape::Row).par_add_vec(&matrix(
                    w.clone(),
                    1000,
                    1000,
                    Shape::Row,
                )),
            )
        })
    });
}

// Check: better parallel results (ran test 6 times)
pub fn par_matrix_norm_lpq_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: [5.5969 ms 5.7555 ms 5.9515 ms 6.0843 ms 6.3072 ms 6.5636 ms]
    cr.bench_function("ser_matrix_norm_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).norm(Norm::Lpq(4.0, 2.0))))
    });

    // Result: 1000x1000 matrix: [3.1796 ms 3.2714 ms 3.3714 ms 3.6123 ms 3.7398 ms 3.8761 ms]
    cr.bench_function("par_matrix_norm_bench", |b| {
        b.iter(|| {
            black_box(matrix(v.clone(), 1000, 1000, Shape::Row).par_norm(Norm::Lpq(4.0, 2.0)))
        })
    });
}

pub fn par_matrix_norm_f_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 1.7414 ms
    cr.bench_function("ser_matrix_norm_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).norm(Norm::F)))
    });

    // Result: 1000x1000 matrix: 9.0499 ms
    cr.bench_function("par_matrix_norm_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).par_norm(Norm::F)))
    });
}

pub fn par_matrix_norm_l1_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 9.0287 ms
    cr.bench_function("ser_matrix_norm_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).norm(Norm::L1)))
    });

    // Result: 1000x1000 matrix: 10.393 ms
    cr.bench_function("par_matrix_norm_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).par_norm(Norm::L1)))
    });
}

// Check: better parallel results (ran test 6 times)
pub fn par_matrix_inner_prod_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    let w: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 3.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix:  [5.1075 ms 5.1505 ms 5.2013 ms 5.7617 ms 6.0196 ms 6.3009 ms]
    cr.bench_function("ser_matrix_inner_prod_bench", |b| {
        b.iter(|| {
            black_box(matrix(v.clone(), 1000, 1000, Shape::Row).dot(&matrix(
                w.clone(),
                1000,
                1000,
                Shape::Row,
            )))
        })
    });

    // Result: 1000x1000 matrix: [4.9931 ms 5.0244 ms 5.0642 ms 5.0322 ms 5.0819 ms 5.1404 ms]
    cr.bench_function("par_matrix_inner_prod_bench", |b| {
        b.iter(|| {
            black_box(matrix(v.clone(), 1000, 1000, Shape::Row).par_dot(&matrix(
                w.clone(),
                1000,
                1000,
                Shape::Row,
            )))
        })
    });
}

pub fn par_matrix_hadamard_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    let w: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 3.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix:  5.8664 ms
    cr.bench_function("ser_matrix_hadamard_bench", |b| {
        b.iter(|| {
            black_box(matrix(v.clone(), 1000, 1000, Shape::Row).hadamard(&matrix(
                w.clone(),
                1000,
                1000,
                Shape::Row,
            )))
        })
    });

    // Result: 1000x1000 matrix: 13.982 ms
    cr.bench_function("par_matrix_hadamard_bench", |b| {
        b.iter(|| {
            black_box(
                matrix(v.clone(), 1000, 1000, Shape::Row).par_hadamard(&matrix(
                    w.clone(),
                    1000,
                    1000,
                    Shape::Row,
                )),
            )
        })
    });
}

pub fn par_matrix_take_row_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();

    // Result: 1000x1000 matrix: 3.7958 ms
    cr.bench_function("ser_matrix_take_row_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).take_row(500)))
    });

    // Result: 1000x1000 matrix: 9.7897 ms
    cr.bench_function("par_matrix_take_row_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_take_row(500)))
    });
}

pub fn par_matrix_fpmap_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();
    let f = |x: f64| 2.0 * x;

    // Result: 1000x1000 matrix: 3.3385 ms
    cr.bench_function("ser_matrix_fpmap_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).fmap(f)))
    });

    // Result: 1000x1000 matrix: 9.6310 ms
    cr.bench_function("par_matrix_fpmap_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_fmap(f)))
    });
}

pub fn par_matrix_reduce_benchmark(cr: &mut Criterion) {
    let v: Vec<f64> = (0..1000000)
        .into_iter()
        .map(|i: i32| 2.0 * (i as f64))
        .collect::<Vec<f64>>();
    let f = |x: f64, y: f64| 2.0 * x * y;

    // Result: 1000x1000 matrix:  2.8077 ms
    cr.bench_function("ser_matrix_reduce_bench", |b| {
        b.iter(|| black_box(matrix(v.clone(), 1000, 1000, Shape::Row).reduce(0.0, f)))
    });

    // Result: 1000x1000 matrix: 9.4877 ms
    cr.bench_function("par_matrix_reduce_bench", |b| {
        b.iter(|| black_box(par_matrix(v.clone(), 1000, 1000, Shape::Row).par_reduce(0.0, f)))
    });
}

criterion_group!(
    benches,
    par_matrix_benchmark,
    par_py_matrix_benchmark,
    par_matrix_change_shape_benchmark,
    par_matrix_extract_row_benchmark,
    par_matrix_from_index_benchmark,
    par_matrix_to_vec_and_diag_benchmark,
    par_matrix_submat_benchmark,
    par_matrix_add_vec_benchmark,
    par_matrix_norm_lpq_benchmark,
    par_matrix_norm_f_benchmark,
    par_matrix_norm_l1_benchmark,
    par_matrix_inner_prod_benchmark,
    par_matrix_hadamard_benchmark,
    par_matrix_take_row_benchmark,
    par_matrix_fpmap_benchmark,
    par_matrix_reduce_benchmark,
);
criterion_main!(benches);
