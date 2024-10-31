use anyhow::Result;
use peroxide::fuga::*;

#[test]
fn test_cubic_root() -> Result<()> {
    let problem = Cubic;
    let bisect = BisectionMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let newton = NewtonMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let false_pos = FalsePositionMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let root_bisect = bisect.find(&problem)?;
    let root_newton = newton.find(&problem)?;
    let root_false_pos = false_pos.find(&problem)?;

    let result_bisect = problem.eval(root_bisect)?[0];
    let result_newton = problem.eval(root_newton)?[0];
    let result_false_pos = problem.eval(root_false_pos)?[0];

    assert!(result_bisect.abs() < 1e-6);
    assert!(result_newton.abs() < 1e-6);
    assert!(result_false_pos.abs() < 1e-6);

    Ok(())
}

struct Cubic;

impl Cubic {
    fn eval(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        Ok([(x[0] - 1f64).powi(3)])
    }
}

impl RootFindingProblem<1, 1, (f64, f64)> for Cubic {
    fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        self.eval(x)
    }
    fn initial_guess(&self) -> (f64, f64) {
        (0.0, 2.0)
    }
}

impl RootFindingProblem<1, 1, f64> for Cubic {
    fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        self.eval(x)
    }
    fn initial_guess(&self) -> f64 {
        0.0
    }
    fn derivative(&self, x: [f64; 1]) -> Result<Jaco<1, 1>> {
        Ok([[3.0 * (x[0] - 1f64).powi(2)]])
    }
}

#[test]
fn test_sine_root() -> Result<()> {
    let problem = Sine;
    let bisect = BisectionMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let newton = NewtonMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let false_pos = FalsePositionMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let root_bisect = bisect.find(&problem)?;
    let root_newton = newton.find(&problem)?;
    let root_false_pos = false_pos.find(&problem)?;

    let result_bisect = problem.eval(root_bisect)?[0];
    let result_newton = problem.eval(root_newton)?[0];
    let result_false_pos = problem.eval(root_false_pos)?[0];

    assert!(result_bisect.abs() < 1e-6);
    assert!(result_newton.abs() < 1e-6);
    assert!(result_false_pos.abs() < 1e-6);

    Ok(())
}

struct Sine;

impl Sine {
    fn eval(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        Ok([x[0].sin()])
    }
}

impl RootFindingProblem<1, 1, (f64, f64)> for Sine {
    fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        self.eval(x)
    }
    fn initial_guess(&self) -> (f64, f64) {
        (0.0, 2.0)
    }
}

impl RootFindingProblem<1, 1, f64> for Sine {
    fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        self.eval(x)
    }
    fn initial_guess(&self) -> f64 {
        1.0
    }
    fn derivative(&self, x: [f64; 1]) -> Result<Jaco<1, 1>> {
        Ok([[x[0].cos()]])
    }
}

#[test]
fn test_cosine_root() -> Result<()> {
    let problem = Cosine;
    let newton = NewtonMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let root_newton = match newton.find(&problem) {
        Ok(x) => x,
        Err(e) => {
            println!("{:?}", e);
            match e.downcast::<RootError<1>>() {
                Ok(RootError::ZeroDerivative(x)) => x,
                Ok(e) => panic!("ok but {:?}", e),
                Err(e) => panic!("err {:?}", e),
            }
        }
    };
    assert_eq!(root_newton[0], 0.0);

    Ok(())
}

struct Cosine;

impl Cosine {
    fn eval(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        Ok([x[0].cos()])
    }
}

impl RootFindingProblem<1, 1, f64> for Cosine {
    fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
        self.eval(x)
    }
    fn initial_guess(&self) -> f64 {
        0.0 // should fail in newton (derivative is 0)
    }
    fn derivative(&self, x: [f64; 1]) -> Result<Jaco<1, 1>> {
        Ok([[-x[0].sin()]])
    }
}
