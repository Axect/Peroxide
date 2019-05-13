extern crate peroxide;
use peroxide::*;

fn main() {
    let a = seq(1,100,1);
    //a.write_single_pickle("example_data/pickle_example_vec.pickle").expect("Can't write pickle");

    let mut b = Matrix::from_index(|i, j| (i + j) as f64, (100, 2));
    match b.shape {
        Row => {
            b = b.change_shape(); // To column shape
        },
        Col => ()
    }

    let c = MATLAB::new("1 2; 3 4");

    //b.write_single_pickle("example_data/pickle_example_mat.pickle").expect("Can't write pickle");
    let mut w = SimpleWriter::new();
    w.insert_header(vec!["x", "y"])
        .insert_matrix(b)
        .insert_matrix(c)
        .insert_vector(a)
        .set_path("example_data/pickle_example.pickle")
        .write_pickle();
}