extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_type_cast() {
    let mut a = DataFrame::new(vec![]);
    a.push("x", Series::new(vec![1, 2, 3, 4]));
    a.push("y", Series::new(vec![true, false, false, true]));

    let mut b = DataFrame::new(vec![]);
    b.push("x", Series::new(vec![1usize, 2, 3, 4]));
    b.push("y", Series::new(vec![true, false, false, true]));

    a.as_types(vec![USIZE, U8]);
    a.as_types(vec![USIZE, Bool]);

    assert_eq!(a, b);
}

// =========================================================================
// Phase 1: Series helpers
// =========================================================================

#[test]
fn test_series_select_indices() {
    let a = Series::new(vec![10, 20, 30, 40, 50]);
    let b = a.select_indices(&[0, 2, 4]);
    assert_eq!(b, Series::new(vec![10, 30, 50]));
}

#[test]
fn test_series_select_indices_string() {
    let a = Series::new(vec!["a".to_string(), "b".to_string(), "c".to_string()]);
    let b = a.select_indices(&[1, 2]);
    assert_eq!(b, Series::new(vec!["b".to_string(), "c".to_string()]));
}

#[test]
fn test_series_select_indices_empty() {
    let a = Series::new(vec![1, 2, 3]);
    let b = a.select_indices(&[]);
    assert_eq!(b.len(), 0);
}

#[test]
fn test_series_to_f64_vec() {
    let a = Series::new(vec![1i32, 2, 3]);
    let v = a.to_f64_vec().unwrap();
    assert_eq!(v, vec![1.0, 2.0, 3.0]);
}

#[test]
fn test_series_to_f64_vec_float() {
    let a = Series::new(vec![1.5f64, 2.5, 3.5]);
    let v = a.to_f64_vec().unwrap();
    assert_eq!(v, vec![1.5, 2.5, 3.5]);
}

#[test]
fn test_series_to_f64_vec_non_numeric() {
    let a = Series::new(vec!["hello".to_string()]);
    assert!(a.to_f64_vec().is_err());
}

// =========================================================================
// Phase 2: Shape / Info
// =========================================================================

fn sample_df() -> DataFrame {
    let mut df = DataFrame::new(vec![]);
    df.push("x", Series::new(vec![1, 2, 3, 4]));
    df.push("y", Series::new(vec![0.1, 0.2, 0.3, 0.4]));
    df.push("z", Series::new(vec!['a', 'b', 'c', 'd']));
    df
}

#[test]
fn test_nrow() {
    let df = sample_df();
    assert_eq!(df.nrow(), 4);
}

#[test]
fn test_ncol() {
    let df = sample_df();
    assert_eq!(df.ncol(), 3);
}

#[test]
fn test_shape() {
    let df = sample_df();
    assert_eq!(df.shape(), (4, 3));
}

#[test]
fn test_dtypes() {
    let df = sample_df();
    assert_eq!(df.dtypes(), vec![I32, F64, Char]);
}

#[test]
fn test_is_empty() {
    let df = sample_df();
    assert!(!df.is_empty());

    let empty_df = DataFrame::new(vec![]);
    assert!(empty_df.is_empty());
}

#[test]
fn test_contains() {
    let df = sample_df();
    assert!(df.contains("x"));
    assert!(df.contains("y"));
    assert!(!df.contains("w"));
}

// =========================================================================
// Phase 3: Row Operations
// =========================================================================

#[test]
fn test_head() {
    let df = sample_df();
    let h = df.head(2);
    assert_eq!(h.nrow(), 2);
    assert_eq!(h.ncol(), 3);
    assert_eq!(h["x"], Series::new(vec![1, 2]));
    assert_eq!(h["y"], Series::new(vec![0.1, 0.2]));
}

#[test]
fn test_head_larger_than_nrow() {
    let df = sample_df();
    let h = df.head(100);
    assert_eq!(h.nrow(), 4);
}

#[test]
fn test_tail() {
    let df = sample_df();
    let t = df.tail(2);
    assert_eq!(t.nrow(), 2);
    assert_eq!(t["x"], Series::new(vec![3, 4]));
    assert_eq!(t["y"], Series::new(vec![0.3, 0.4]));
}

#[test]
fn test_tail_larger_than_nrow() {
    let df = sample_df();
    let t = df.tail(100);
    assert_eq!(t.nrow(), 4);
}

#[test]
fn test_slice() {
    let df = sample_df();
    let s = df.slice(1, 2);
    assert_eq!(s.nrow(), 2);
    assert_eq!(s["x"], Series::new(vec![2, 3]));
}

#[test]
fn test_slice_beyond_end() {
    let df = sample_df();
    let s = df.slice(2, 100);
    assert_eq!(s.nrow(), 2);
    assert_eq!(s["x"], Series::new(vec![3, 4]));
}

// =========================================================================
// Phase 4: Column Operations
// =========================================================================

#[test]
fn test_select() {
    let df = sample_df();
    let selected = df.select(&["x", "z"]);
    assert_eq!(selected.ncol(), 2);
    assert_eq!(selected.column_names(), vec!["x", "z"]);
    assert_eq!(selected["x"], Series::new(vec![1, 2, 3, 4]));
}

#[test]
#[should_panic(expected = "Column 'w' not found")]
fn test_select_missing_column() {
    let df = sample_df();
    df.select(&["w"]);
}

#[test]
fn test_rename() {
    let mut df = sample_df();
    df.rename("x", "a");
    assert!(df.contains("a"));
    assert!(!df.contains("x"));
    assert_eq!(df["a"], Series::new(vec![1, 2, 3, 4]));
}

#[test]
#[should_panic(expected = "Column 'w' not found")]
fn test_rename_missing_column() {
    let mut df = sample_df();
    df.rename("w", "a");
}

#[test]
fn test_column_names() {
    let df = sample_df();
    assert_eq!(df.column_names(), vec!["x", "y", "z"]);
}

#[test]
fn test_select_dtypes() {
    let df = sample_df();
    let numeric = df.select_dtypes(&[I32, F64]);
    assert_eq!(numeric.ncol(), 2);
    assert_eq!(numeric.column_names(), vec!["x", "y"]);
}

#[test]
fn test_select_dtypes_no_match() {
    let df = sample_df();
    let result = df.select_dtypes(&[Bool]);
    assert_eq!(result.ncol(), 0);
}

// =========================================================================
// Phase 5: Series Statistics
// =========================================================================

#[test]
fn test_series_sum() {
    let a = Series::new(vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(a.sum().unwrap(), 10.0);
}

#[test]
fn test_series_sum_int() {
    let a = Series::new(vec![1i32, 2, 3, 4]);
    assert_eq!(a.sum().unwrap(), 10.0);
}

#[test]
fn test_series_mean() {
    let a = Series::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    assert_eq!(a.mean().unwrap(), 3.0);
}

#[test]
fn test_series_var() {
    let a = Series::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    assert_eq!(a.var().unwrap(), 2.5);
}

#[test]
fn test_series_sd() {
    let a = Series::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let sd = a.sd().unwrap();
    assert!((sd - 2.5f64.sqrt()).abs() < 1e-10);
}

#[test]
fn test_series_min() {
    let a = Series::new(vec![3, 1, 4, 1, 5]);
    let m = a.min().unwrap();
    assert_eq!(m, Scalar::new(1i32));
}

#[test]
fn test_series_min_f64() {
    let a = Series::new(vec![3.0, 1.5, 4.0]);
    let m = a.min().unwrap();
    assert_eq!(m, Scalar::new(1.5f64));
}

#[test]
fn test_series_max() {
    let a = Series::new(vec![3, 1, 4, 1, 5]);
    let m = a.max().unwrap();
    assert_eq!(m, Scalar::new(5i32));
}

#[test]
fn test_series_max_string() {
    let a = Series::new(vec!["apple".to_string(), "banana".to_string(), "cherry".to_string()]);
    let m = a.max().unwrap();
    assert_eq!(m, Scalar::new("cherry".to_string()));
}

#[test]
fn test_series_stats_non_numeric() {
    let a = Series::new(vec!["hello".to_string()]);
    assert!(a.sum().is_err());
    assert!(a.mean().is_err());
}

// =========================================================================
// Phase 6: DataFrame Statistics
// =========================================================================

#[test]
fn test_df_describe() {
    let mut df = DataFrame::new(vec![]);
    df.push("a", Series::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]));
    df.push("b", Series::new(vec![10.0, 20.0, 30.0, 40.0, 50.0]));
    df.push("c", Series::new(vec!['x', 'y', 'z', 'w', 'v']));

    let desc = df.describe();
    // Should have stat column + 2 numeric columns
    assert_eq!(desc.ncol(), 3);
    assert!(desc.contains("stat"));
    assert!(desc.contains("a"));
    assert!(desc.contains("b"));
    assert!(!desc.contains("c")); // char column excluded

    // Check count row (first data row)
    let a_col: Vec<f64> = desc["a"].to_vec();
    assert_eq!(a_col[0], 5.0); // count
    assert_eq!(a_col[1], 3.0); // mean
}

#[test]
fn test_df_sum() {
    let mut df = DataFrame::new(vec![]);
    df.push("a", Series::new(vec![1.0, 2.0, 3.0]));
    df.push("b", Series::new(vec![10, 20, 30]));
    df.push("c", Series::new(vec!['x', 'y', 'z']));

    let s = df.sum();
    assert_eq!(s.ncol(), 2); // only numeric columns
    let a_sum: Vec<f64> = s["a"].to_vec();
    assert_eq!(a_sum[0], 6.0);
    let b_sum: Vec<f64> = s["b"].to_vec();
    assert_eq!(b_sum[0], 60.0);
}

#[test]
fn test_df_mean() {
    let mut df = DataFrame::new(vec![]);
    df.push("a", Series::new(vec![1.0, 2.0, 3.0]));
    df.push("b", Series::new(vec![10, 20, 30]));

    let m = df.mean();
    assert_eq!(m.ncol(), 2);
    let a_mean: Vec<f64> = m["a"].to_vec();
    assert_eq!(a_mean[0], 2.0);
    let b_mean: Vec<f64> = m["b"].to_vec();
    assert_eq!(b_mean[0], 20.0);
}
