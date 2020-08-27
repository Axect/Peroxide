extern crate peroxide;
use peroxide::fuga::*;

fn main() -> Result<(), Box<dyn Error>> {
    #[cfg(feature="dataframe")]
    {
        let x = seq(1, 10, 1);
        let y = x.fmap(|t| t.powi(2));
        let mut df = DataFrame::with_header(vec!["x", "y"]);
        df["x"] = x;
        df["y"] = y;
        df.write_nc("example_data/nc_test.nc")?;

        let dg = DataFrame::read_nc("example_data/nc_test.nc")?;
        dg.head_print(4);
        dg.tail_print(4);
    }
    Ok(())
}
