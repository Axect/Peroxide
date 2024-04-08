extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_stats() -> Result<(), Box<dyn std::error::Error>> {
    let intervals = vec![0f64, 3f64, 6f64];
    let weights = vec![2f64, 1f64];
    let unif = WeightedUniform::new(weights, intervals)?;
    assert_eq!(unif.mean(), 2.5);
    assert_eq!(unif.var(), 2.75);

    Ok(())
}

#[test]
fn test_max_pool_1d() -> Result<(), Box<dyn std::error::Error>> {
    let w = WeightedUniform::from_max_pool_1d(f, (-2f64, 3f64), 10, 1e-5)?;
    println!("{:?}", w.intervals());
    let x = linspace(w.intervals()[0].0, w.intervals()[8].1, 11);
    let y1 = x.fmap(f);
    let y2 = x.fmap(|t| w.weight_at(t));
    for (a, b) in y1.iter().zip(y2.iter()) {
        println!("a: {}, b: {}", a, b);
        assert!(a <= b);
    }

    Ok(())
}

fn f(x: f64) -> f64 {
    if x.abs() < 1f64 {
        1f64 - x.abs()
    } else {
       0f64
    }
}
