extern crate peroxide;
use peroxide::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    #[cfg(feature = "dataframe")]
    {
        let mut a: DataFrame = DataFrame::new();
        a.insert("x", c!(1,2,3,4,5));
        a.insert("y", c!(4,5,6));
        a.insert("z", c!(7,8,9));
        a.print();
        a.write_csv("example_data/df_test.csv")?;
        a.write_cdf("example_data/df_test.cdf")?;

        let b = DataFrame::read_csv("example_data/df_test.csv", ',')?;
        b.print();
    }
    Ok(())
}
