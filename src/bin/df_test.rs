extern crate peroxide;
use peroxide::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    #[cfg(feature = "dataframe")]
    {
        let mut a = DataFrame::new();

        a.insert("x", vec![1f64, 2f64, 3f64]);
        a.insert("y", vec![4f64, 5f64]);
        a.insert("z", vec![6f64]);

        a.print();
    }
    Ok(())
}
