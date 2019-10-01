extern crate peroxide;
use peroxide::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
    df["x"] = vec![0f64; 100_0000];
    df["y"] = vec![0f64; 100_0000];
    df["z"] = vec![0f64; 100_0000];

    df.write_csv("example_data/df_bench.csv")?;
    //df.write_cdf("example_data/df_bench.cdf")?;


    let dg = DataFrame::read_csv("example_data/df_bench.csv", ',')?;
    //let dg = DataFrame::read_cdf("example_data/df_bench.cdf", vec!["x", "y", "z"])?;
    dg.print();

    Ok(())
}
