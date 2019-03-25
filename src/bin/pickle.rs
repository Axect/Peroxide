extern crate peroxide;
use peroxide::*;

fn main() {
    let a = seq(1,100,1);
    a.write_pickle("example_data/pickle_example_vec").expect("Can't write pickle");

    let mut b = Matrix::from_index(|i, j| i + j, (100, 2));
    match b.shape {
        Row => {
            b = b.change_shape(); // To column shape
        },
        Col => ()
    }

    b.write_pickle("example_data/pickle_example_mat").expect("Can't write pickle");
}