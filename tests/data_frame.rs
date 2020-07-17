extern crate peroxide;
#[allow(unused_imports)]
use peroxide::fuga::*;
// use std::fs::remove_file;

#[test]
#[cfg(feature = "dataframe")]
fn conversion_with_matrix() {
    let mat1 = matrix(seq(1, 6, 1), 3, 2, Col);
    let df = DataFrame::from_matrix(mat1.clone());
    let mat2 = df.to_matrix();
    assert_eq!(mat1, mat2);
}

#[test]
#[cfg(feature = "dataframe")]
fn set_header_test() {
    let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
    df.set_header(vec!["a", "b", "c"]);
    let head = df.headers().map(|x| x.as_str()).collect::<Vec<&str>>();
    assert_eq!(head, vec!["a", "b", "c"]);
}

//#[test]
//#[cfg(feature = "dataframe")]
//fn read_write_nc_test() {
//    let m = matrix(seq(1,6,1), 3, 2, Col);
//
//    // Prepare dataframe
//    let mut df = DataFrame::from_matrix(m);
//    df.set_header(vec!["x", "y"]);
//    df.print();
//
//    // Write NetCDF
//    match df.write_nc("example_data/ex.nc") {
//        Ok(_) => (),
//        Err(_) => panic!("Can't write dataframe"),
//    }
//
//    // Read NetCDF
//    let dg = match DataFrame::read_nc("example_data/ex.nc") {
//        Ok(m) => m,
//        Err(_) => panic!("Can't read dataframe"),
//    };
//    dg.print();
//
//    assert_eq!(df, dg);
//
//    remove_file("example_data/ex.nc");
//}
