use peroxide::fuga::*;
use peroxide::numerical::root::{Pt, Intv};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let problem = Quadratic;
    let broyden = BroydenMethod { max_iter: 100, tol: 1e-6 };
    let root = broyden.find(&problem)?;
    let result = problem.function(root)?;

    println!("root: {:?}", root);
    println!("result: {:?}", result);

    Ok(())
}

struct Quadratic;

impl RootFindingProblem<2, 2, Intv<2>> for Quadratic {
    fn function(&self, x: Pt<2>) -> anyhow::Result<Pt<2>> {
        Ok([
            x[0] * x[0] + x[1] * x[1] - 1.0,
            x[0] + x[1] - 2f64.sqrt()
        ])
    }

    fn initial_guess(&self) -> Intv<2> {
        ([0.0, 0.1], [-0.1, 0.2])
    }
}
