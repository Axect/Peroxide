//! Plotting module for peroxide
//!
//! For Rust, there are some plot libraries but, still difficult to use.
//! Practically, using python is best choice to plot. And there is awesome crate - [pyo3](https://crates.io/crates/pyo3).
//!
//! Let's see next ordinary code file.
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::*;
//!
//! fn main() {
//!     let init_state = State::<f64>::new(0f64, c!(1), c!(0));
//!
//!     let mut ode_solver = ExplicitODE::new(test_fn);
//!
//!     ode_solver
//!         .set_method(ExMethod::RK4)
//!         .set_initial_condition(init_state)
//!         .set_step_size(0.01)
//!         .set_times(1000);
//!
//!     let result = ode_solver.integrate();
//!
//!     let mut st = SimpleWriter::new();
//!     st.set_path("example_data/rk4_test.pickle")
//!         .insert_matrix(result)
//!         .write_pickle();
//! }
//!
//! fn test_fn(st: &mut State<f64>) {
//!     let x = st.param;
//!     let y = &st.value;
//!     let dy = &mut st.deriv;
//!     dy[0] = (5f64*x.powi(2) - y[0]) / (x + y[0]).exp();
//! }
//! ```
//!
//! Now, let's modify this code to below. Then it works surprisingly!
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::*;
//!
//! fn main() {
//!     let init_state = State::<f64>::new(0f64, c!(1), c!(0));
//!
//!     let mut ode_solver = ExplicitODE::new(test_fn);
//!
//!     ode_solver
//!         .set_method(ExMethod::RK4)
//!         .set_initial_condition(init_state)
//!         .set_step_size(0.01)
//!         .set_times(1000);
//!
//!     let result = ode_solver.integrate();
//!
//!     let x = result.col(0);
//!     let y = result.col(1);
//!
//!     let mut plt = Plot2D::new();
//!     plt.set_domain(x)
//!         .insert_image(y)
//!         .set_title("Test Figure")
//!         .set_fig_size((10, 6))
//!         .set_dpi(300)
//!         .set_legends(vec!["RK4"])
//!         .set_path("example_data/test_plot.png");
//!
//!     plt.savefig();
//! }
//!
//! fn test_fn(st: &mut State<f64>) {
//!     let x = st.param;
//!     let y = &st.value;
//!     let dy = &mut st.deriv;
//!     dy[0] = (5f64 * x.powi(2) - y[0]) / (x + y[0]).exp();
//! }
//! ```
//!
//! It draws next image
//!
//! ![test_plot](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/test_plot.png)


extern crate pyo3;
use self::pyo3::types::IntoPyDict;
use self::pyo3::{PyResult, Python};
use std::collections::HashMap;
pub use Grid::*;
use PlotOptions::{Domain, Images, Legends, Path};

type Vector = Vec<f64>;

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum PlotOptions {
    Domain,
    Images,
    Legends,
    Path,
}

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum Grid {
    On,
    Off,
}

pub trait Plot {
    fn set_domain(&mut self, x: Vec<f64>) -> &mut Self;
    fn insert_image(&mut self, y: Vec<f64>) -> &mut Self;
    fn set_title(&mut self, title: &str) -> &mut Self;
    fn set_xlabel(&mut self, xlabel: &str) -> &mut Self;
    fn set_ylabel(&mut self, ylabel: &str) -> &mut Self;
    fn set_zlabel(&mut self, zlabel: &str) -> &mut Self;
    fn set_legends(&mut self, legends: Vec<&str>) -> &mut Self;
    fn set_path(&mut self, path: &str) -> &mut Self;
    fn set_fig_size(&mut self, fig_size: (usize, usize)) -> &mut Self;
    fn set_dpi(&mut self, dpi: usize) -> &mut Self;
    fn grid(&mut self, grid: Grid) -> &mut Self;
    fn savefig(&self) -> PyResult<()>;
}

#[derive(Debug)]
pub struct Plot2D {
    domain: Vector,
    images: Vec<Vector>,
    title: String,
    xlabel: String,
    ylabel: String,
    legends: Vec<String>,
    path: String,
    fig_size: (usize, usize),
    dpi: usize,
    grid: Grid,
    options: HashMap<PlotOptions, bool>,
}

impl Plot2D {
    pub fn new() -> Self {
        let mut default_options: HashMap<PlotOptions, bool> = HashMap::new();
        default_options.insert(Domain, false);
        default_options.insert(Images, false);
        default_options.insert(Legends, false);
        default_options.insert(Path, false);

        Plot2D {
            domain: vec![],
            images: vec![],
            title: "Title".to_string(),
            xlabel: "$x$".to_string(),
            ylabel: "$y$".to_string(),
            legends: vec![],
            path: "".to_string(),
            fig_size: (10, 6),
            dpi: 300,
            grid: On,
            options: default_options,
        }
    }
}

impl Plot for Plot2D {
    fn set_domain(&mut self, x: Vec<f64>) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Domain) {
            *x = true
        }
        self.domain = x;
        self
    }

    fn insert_image(&mut self, y: Vec<f64>) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Images) {
            *x = true
        }
        self.images.push(y);
        self
    }

    fn set_title(&mut self, title: &str) -> &mut Self {
        self.title = title.to_owned();
        self
    }

    fn set_xlabel(&mut self, xlabel: &str) -> &mut Self {
        self.xlabel = xlabel.to_owned();
        self
    }

    fn set_ylabel(&mut self, ylabel: &str) -> &mut Self {
        self.ylabel = ylabel.to_owned();
        self
    }

    fn set_zlabel(&mut self, _zlabel: &str) -> &mut Self {
        unimplemented!()
    }

    fn set_legends(&mut self, legends: Vec<&str>) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Legends) {
            *x = true
        }
        self.legends = legends.into_iter().map(|x| x.to_owned()).collect::<Vec<String>>();
        self
    }

    fn set_path(&mut self, path: &str) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Path) {
            *x = true
        }
        self.path = path.to_owned();
        self
    }

    fn set_fig_size(&mut self, fig_size: (usize, usize)) -> &mut Self {
        self.fig_size = fig_size;
        self
    }

    fn set_dpi(&mut self, dpi: usize) -> &mut Self {
        self.dpi = dpi;
        self
    }

    fn grid(&mut self, grid: Grid) -> &mut Self {
        self.grid = grid;
        self
    }

    fn savefig(&self) -> PyResult<()> {
        // Check domain
        match self.options.get(&Domain) {
            Some(x) => {
                assert!(*x, "Domain is not defined");
            }
            None => {
                assert!(false, "Domain is None");
            }
        }

        // Check images
        match self.options.get(&Images) {
            Some(x) => {
                assert!(*x, "Images are not defined");
            }
            None => {
                assert!(false, "Images are None");
            }
        }

        // Check legends
        match self.options.get(&Legends) {
            Some(x) => {
                assert!(*x, "Legends are not defined");
                assert!(
                    self.images.len() == self.legends.len(),
                    "Legends are not matched with images"
                );
            }
            None => {
                assert!(false, "Legends are None");
            }
        }

        // Plot
        let gil = Python::acquire_gil();
        let py = gil.python();

        // Input data
        let x = self.domain.clone();
        let ys = self.images.clone();
        let y_length = ys.len();
        let title = self.title.clone();
        let fig_size = self.fig_size;
        let dpi = self.dpi;
        let grid = match self.grid {
            On => true,
            Off => false,
        };
        let xlabel = self.xlabel.clone();
        let ylabel = self.ylabel.clone();
        let legends = self.legends.clone();
        let path = self.path.clone();

        // Global variables to plot
        let globals = vec![("plt", py.import("pylab")?)].into_py_dict(py);
        globals.set_item("x", x)?;
        globals.set_item("y", ys)?;
        globals.set_item("n", y_length)?;
        globals.set_item("fs", fig_size)?;
        globals.set_item("dp", dpi)?;
        globals.set_item("gr", grid)?;
        globals.set_item("pa", path)?;

        // Plot Code
        let mut plot_string = format!(
            "\
             plt.rc(\"text\", usetex=True)\n\
             plt.rc(\"font\", family=\"serif\")\n\
             plt.figure(figsize=fs, dpi=dp)\n\
             plt.title(r\"{}\", fontsize=16)\n\
             plt.xlabel(r\"{}\", fontsize=14)\n\
             plt.ylabel(r\"{}\", fontsize=14)\n\
             if gr:\n\
             \tplt.grid()\n",
            title, xlabel, ylabel
        );

        for i in 0..y_length {
            plot_string.push_str(&format!("plt.plot(x,y[{}],label=r\"{}\")\n", i, legends[i])[..])
        }

        plot_string.push_str("plt.legend(fontsize=12)\nplt.savefig(pa)");

        py.run(&plot_string[..], Some(&globals), None)?;
        Ok(())
    }
}
