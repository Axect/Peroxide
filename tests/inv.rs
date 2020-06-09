extern crate peroxide;
#[allow(unused_imports)]
use peroxide::fuga::*;

#[cfg(feature = "dataframe")]
#[test]
fn test_inverse() {
    for i in 6 .. 11 {
        let df = DataFrame::read_nc_by_header(&format!("test_data/inv/inverse_{}.nc", i), vec!["a", "b"]).unwrap();
        let a: Matrix = matrix(df["a"].clone(), i, i, Col);
        let b: Matrix = matrix(df["b"].clone(), i, i, Col);
        let c = a.inv().unwrap();
        assert_eq!(b, c);
    }
}
