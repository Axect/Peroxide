use criterion::{black_box, criterion_group, criterion_main, Criterion};
use peroxide::fuga::*;

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

criterion_group!(
    benches,
    /*par_matrix_benchmark,*/ par_py_matrix_benchmark
);
criterion_main!(benches);
