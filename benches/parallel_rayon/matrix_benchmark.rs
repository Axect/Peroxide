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

criterion_group!(
    benches,
    /*par_matrix_benchmark,*/
    /*par_py_matrix_benchmark,*/
    /*par_matrix_change_shape_benchmark,*/
    /*par_matrix_extract_row_benchmark,*/
    /*par_matrix_from_index_benchmark,*/
    /*par_matrix_to_vec_and_diag_benchmark,*/
    par_matrix_submat_benchmark,
);
criterion_main!(benches);
