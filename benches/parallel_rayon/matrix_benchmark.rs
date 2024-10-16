use criterion::{black_box, criterion_group, criterion_main, Criterion};
use peroxide::{
    fuga::*,
    traits::math::{ParallelInnerProduct, ParallelNormed},
};

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

criterion_group!(
    benches,
    par_matrix_from_index_benchmark,
    par_matrix_norm_lpq_benchmark,
    par_matrix_norm_l1_benchmark,
    par_matrix_inner_prod_benchmark,
);
criterion_main!(benches);
