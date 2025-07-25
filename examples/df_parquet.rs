use peroxide::fuga::*;

fn main() -> Result<(), Box<dyn Error>> {
    #[cfg(feature = "parquet")]
    {
        let a = Series::new(vec![1, 2, 3, 4]);
        let b = Series::new(vec![0.1, 0.2, 0.3, 0.4]);
        let c = Series::new(vec![true, false, false, true]);
        let d = Series::new(
            vec!["a", "b", "c", "d"]
                .into_iter()
                .map(|x| x.to_string())
                .collect(),
        );

        let mut df = DataFrame::new(vec![a, b, c, d]);
        df.set_header(vec!["a", "b", "c", "d"]);
        println!("Write:");
        df.print();

        df.write_parquet("example_data/df_parquet.parquet", UNCOMPRESSED)?;

        println!("\nRead:");
        let mut dg = DataFrame::read_parquet("example_data/df_parquet.parquet")?;
        dg.print();

        println!("\nConvert:");
        dg[2].as_type(Bool);
        dg.print();
    }

    Ok(())
}
