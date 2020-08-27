extern crate peroxide;
use peroxide::fuga::*;
use std::error::Error;
use std::fs;

fn main() -> Result<(), Box<dyn Error>> {
    #[cfg(feature = "dataframe")]
    {
        let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
        df["x"] = vec![0f64; 1000_000];
        df["y"] = vec![0f64; 1000_000];
        df["z"] = vec![0f64; 1000_000];

        df.write_nc("example_data/df_bench.nc")?;

        // read_nc(path: &str, header: Vec<&str>)
        let dg = DataFrame::read_nc("example_data/df_bench.nc")?;
        // or read_nc with specific header
        let dg = DataFrame::read_nc_by_header("example_data/df_bench.nc", vec!["x", "z"])?;
        dg.print();

        // For test, we remove nc file
        fs::remove_file("example_data/df_bench.nc")?;
    }

    Ok(())
}
