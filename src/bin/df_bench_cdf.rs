extern crate peroxide;
use peroxide::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
    df["x"] = vec![0f64; 100_0000];
    df["y"] = vec![0f64; 100_0000];
    df["z"] = vec![0f64; 100_0000];

    df.write_nc("example_data/df_bench.nc")?;


    let dg = DataFrame::read_nc("example_data/df_bench.nc", vec!["x", "y", "z"])?;
    dg.print();

    Ok(())
}
