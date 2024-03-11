//! Plotting module for peroxide
//!
//! For Rust, there are some plot libraries but, still difficult to use.
//! Practically, using python is best choice to plot. And there is awesome crate - [pyo3](https://crates.io/crates/pyo3).
//!
//! Let's see next ordinary code file.
//!
//! ```no-run
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
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
//!     result.write("example_data/test.csv");
//! }
//!
//! fn test_fn(st: &mut State<f64>, _: &NoEnv) {
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
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
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
//!     // Remove below comments to execute
//!     //let mut plt = Plot2D::new();
//!     //plt.set_domain(x)
//!     //    .insert_image(y)
//!     //    .set_title("Test Figure")
//!     //    .set_fig_size((10, 6))
//!     //    .set_dpi(300)
//!     //    .set_legend(vec!["RK4"])
//!     //    .set_path("example_data/test_plot.png");
//!
//!     //plt.savefig();
//! }
//!
//! fn test_fn(st: &mut State<f64>, _: &NoEnv) {
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
//!
//! But now, the recommended way is exporting `netcdf` files. Refer to [dataframe](../../structure/dataframe/index.html)

extern crate pyo3;
use self::pyo3::types::IntoPyDict;
use self::pyo3::{PyResult, Python};
pub use self::Grid::{Off, On};
pub use self::Markers::{Circle, Line, Point};
use self::PlotOptions::{Domain, Images, Legends, Pairs, Path};
use std::collections::HashMap;
use std::fmt::Display;

type Vector = Vec<f64>;

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum PlotOptions {
    Domain,
    Images,
    Pairs,
    Legends,
    Path,
}

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum Markers {
    Point,
    Line,
    Circle,
    Pixel,
    TriangleDown,
    TriangleUp,
    TriangleLeft,
    TriangleRight,
    Square,
    Pentagon,
    Star,
    Hexagon1,
    Hexagon2,
    Plus,
    X,
    Diamond,
    ThinDiamond,
    VLine,
    HLine,
}

impl Display for Markers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Markers::Point => ".".to_string(),
            Markers::Line => "-".to_string(),
            Markers::Circle => "o".to_string(),
            Markers::Pixel => ",".to_string(),
            Markers::TriangleDown => "v".to_string(),
            Markers::TriangleUp => "^".to_string(),
            Markers::TriangleLeft => "<".to_string(),
            Markers::TriangleRight => ">".to_string(),
            Markers::Square => "s".to_string(),
            Markers::Pentagon => "p".to_string(),
            Markers::Star => "*".to_string(),
            Markers::Hexagon1 => "h".to_string(),
            Markers::Hexagon2 => "H".to_string(),
            Markers::Plus => "+".to_string(),
            Markers::X => "x".to_string(),
            Markers::Diamond => "D".to_string(),
            Markers::ThinDiamond => "d".to_string(),
            Markers::VLine => "|".to_string(),
            Markers::HLine => "_".to_string(),
        };
        write!(f, "{}", str)
    }
}

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum LineStyle {
    Solid,
    Dashed,
    Dotted,
    DashDot,
}

impl Display for LineStyle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            LineStyle::Solid => "solid".to_string(),
            LineStyle::Dashed => "dashed".to_string(),
            LineStyle::Dotted => "dotted".to_string(),
            LineStyle::DashDot => "dashdot".to_string(),
        };
        write!(f, "{}", str)
    }
}

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum Grid {
    On,
    Off,
}

/// Plot Style (`scienceplots` should be installed)
///
/// * Nature
/// * IEEE
/// * Default (Matplotlib default style)
/// * Science
#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum PlotStyle {
    Nature,
    IEEE,
    Default,
    Science,
}

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum PlotScale {
    Linear,
    Log,
}

pub trait Plot {
    fn set_domain(&mut self, x: Vec<f64>) -> &mut Self;
    fn insert_image(&mut self, y: Vec<f64>) -> &mut Self;
    fn insert_pair(&mut self, xy: (Vec<f64>, Vec<f64>)) -> &mut Self;
    fn set_title(&mut self, title: &str) -> &mut Self;
    fn set_xlabel(&mut self, xlabel: &str) -> &mut Self;
    fn set_ylabel(&mut self, ylabel: &str) -> &mut Self;
    fn set_zlabel(&mut self, zlabel: &str) -> &mut Self;
    fn set_xscale(&mut self, xscale: PlotScale) -> &mut Self;
    fn set_yscale(&mut self, yscale: PlotScale) -> &mut Self;
    fn set_xlim(&mut self, xlim: (f64, f64)) -> &mut Self;
    fn set_ylim(&mut self, ylim: (f64, f64)) -> &mut Self;
    fn set_legend(&mut self, legends: Vec<&str>) -> &mut Self;
    fn set_path(&mut self, path: &str) -> &mut Self;
    fn set_fig_size(&mut self, fig_size: (usize, usize)) -> &mut Self;
    fn set_dpi(&mut self, dpi: usize) -> &mut Self;
    fn grid(&mut self, grid: Grid) -> &mut Self;
    fn set_marker(&mut self, styles: Vec<Markers>) -> &mut Self;
    fn set_style(&mut self, style: PlotStyle) -> &mut Self;
    fn tight_layout(&mut self) -> &mut Self;
    fn set_line_style(&mut self, style: Vec<LineStyle>) -> &mut Self;
    fn set_color(&mut self, color: Vec<&str>) -> &mut Self;
    fn set_alpha(&mut self, alpha: Vec<f64>) -> &mut Self;
    fn savefig(&self) -> PyResult<()>;
}

#[derive(Debug)]
pub struct Plot2D {
    domain: Vector,
    images: Vec<Vector>,
    pairs: Vec<(Vector, Vector)>,
    title: Option<String>,
    xlabel: Option<String>,
    ylabel: Option<String>,
    xscale: PlotScale,
    yscale: PlotScale,
    xlim: Option<(f64, f64)>,
    ylim: Option<(f64, f64)>,
    legends: Vec<String>,
    markers: Vec<Markers>,
    line_style: Vec<LineStyle>,
    color: Vec<String>,
    alpha: Vec<f64>,
    path: String,
    fig_size: Option<(usize, usize)>,
    dpi: usize,
    grid: Grid,
    style: PlotStyle,
    tight: bool,
    options: HashMap<PlotOptions, bool>,
}

impl Plot2D {
    pub fn new() -> Self {
        let mut default_options: HashMap<PlotOptions, bool> = HashMap::new();
        default_options.insert(Domain, false);
        default_options.insert(Images, false);
        default_options.insert(Pairs, false);
        default_options.insert(Path, false);

        Plot2D {
            domain: vec![],
            images: vec![],
            pairs: vec![],
            title: None,
            xlabel: None,
            ylabel: None,
            xscale: PlotScale::Linear,
            yscale: PlotScale::Linear,
            xlim: None,
            ylim: None,
            legends: vec![],
            markers: vec![],
            line_style: vec![],
            color: vec![],
            alpha: vec![],
            path: "".to_string(),
            fig_size: None,
            dpi: 300,
            grid: On,
            style: PlotStyle::Default,
            tight: false,
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

    fn insert_pair(&mut self, xy: (Vec<f64>, Vec<f64>)) -> &mut Self {
        if let Some(t) = self.options.get_mut(&Pairs) {
            *t = true
        }
        self.pairs.push(xy);
        self
    }

    fn set_title(&mut self, title: &str) -> &mut Self {
        self.title = Some(title.to_owned());
        self
    }

    fn set_xlabel(&mut self, xlabel: &str) -> &mut Self {
        self.xlabel = Some(xlabel.to_owned());
        self
    }

    fn set_ylabel(&mut self, ylabel: &str) -> &mut Self {
        self.ylabel = Some(ylabel.to_owned());
        self
    }

    fn set_zlabel(&mut self, _zlabel: &str) -> &mut Self {
        unimplemented!()
    }

    fn set_xscale(&mut self, xscale: PlotScale) -> &mut Self {
        self.xscale = xscale;
        self
    }

    fn set_yscale(&mut self, yscale: PlotScale) -> &mut Self {
        self.yscale = yscale;
        self
    }

    fn set_xlim(&mut self, xlim: (f64, f64)) -> &mut Self {
        self.xlim = Some(xlim);
        self
    }

    fn set_ylim(&mut self, ylim: (f64, f64)) -> &mut Self {
        self.ylim = Some(ylim);
        self
    }

    fn set_legend(&mut self, legends: Vec<&str>) -> &mut Self {
        self.legends = legends
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
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
        self.fig_size = Some(fig_size);
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

    fn set_marker(&mut self, styles: Vec<Markers>) -> &mut Self {
        self.markers = styles;
        self
    }

    fn set_style(&mut self, style: PlotStyle) -> &mut Self {
        self.style = style;
        self
    }

    fn tight_layout(&mut self) -> &mut Self {
        self.tight = true;
        self
    }

    fn set_line_style(&mut self, style: Vec<LineStyle>) -> &mut Self {
        self.line_style = style;
        self
    }

    fn set_color(&mut self, color: Vec<&str>) -> &mut Self {
        self.color = color.into_iter().map(|x| x.to_owned()).collect();
        self
    }

    fn set_alpha(&mut self, alpha: Vec<f64>) -> &mut Self {
        self.alpha = alpha;
        self
    }

    fn savefig(&self) -> PyResult<()> {
        // Check domain
        match self.options.get(&Domain) {
            Some(x) if !*x => match self.options.get(&Pairs) {
                Some(xy) if !*xy => {
                    panic!("There are no data to plot");
                }
                None => {
                    panic!("There are some serious problems in plot system");
                }
                _ => (),
            },
            None => {
                panic!("There are some serious problems in plot system");
            }
            _ => (),
        }

        // Check images
        match self.options.get(&Images) {
            Some(x) if !*x => match self.options.get(&Pairs) {
                Some(xy) if !*xy => {
                    panic!("there are no data to plot");
                }
                None => {
                    panic!("There are some serious problems in plot system");
                }
                _ => (),
            },
            None => {
                panic!("There are some serious problems in plot system");
            }
            _ => (),
        }

        // Plot
        Python::with_gil(|py| {
            // Input data
            let x = self.domain.clone();
            let ys = self.images.clone();
            let pairs = self.pairs.clone();
            let y_length = ys.len();
            let pair_length = pairs.len();
            let title = self.title.clone();
            let fig_size = self.fig_size;
            let dpi = self.dpi;
            let grid = match self.grid {
                On => true,
                Off => false,
            };
            let style = match self.style {
                PlotStyle::Nature => "nature",
                PlotStyle::IEEE => "ieee",
                PlotStyle::Default => "default",
                PlotStyle::Science => "science",
            };
            let xlabel = self.xlabel.clone();
            let ylabel = self.ylabel.clone();
            let legends = self.legends.clone();
            let path = self.path.clone();
            let markers = self.markers.iter().map(|x| format!("{}", x)).collect::<Vec<String>>();
            let line_style = self.line_style.iter().map(|x| format!("{}", x)).collect::<Vec<String>>();
            let color = self.color.clone();
            let alpha = self.alpha.clone();

            // Global variables to plot
            let globals = vec![("plt", py.import("matplotlib.pyplot")?)].into_py_dict(py);
            globals.set_item("x", x)?;
            globals.set_item("y", ys)?;
            globals.set_item("pair", pairs)?;
            globals.set_item("n", y_length)?;
            globals.set_item("p", pair_length)?;
            if let Some(fs) = fig_size {
                globals.set_item("fs", fs)?;
            }
            globals.set_item("dp", dpi)?;
            globals.set_item("gr", grid)?;
            globals.set_item("pa", path)?;
            if let Some(xl) = self.xlim {
                globals.set_item("xl", xl)?;
            }
            if let Some(yl) = self.ylim {
                globals.set_item("yl", yl)?;
            }

            // Plot Code
            let mut plot_string = match self.style {
                PlotStyle::Default => {
                    "\
                    plt.rc(\"text\", usetex=True)\n\
                    plt.rc(\"font\", family=\"serif\")\n".to_string()
                }
                PlotStyle::Science => {
                    "\
                    import scienceplots\n\
                    plt.style.use(\"science\")\n".to_string()
                }
                _ => format!(
                    "\
                    import scienceplots\n\
                    plt.style.use([\"science\", \"{}\"])\n",
                    style
                ),
            };
            if let Some(fs) = fig_size {
                plot_string.push_str(&format!("plt.figure(figsize=fs, dpi=dp)\n")[..]);
            } else {
                plot_string.push_str(&format!("plt.figure()\n")[..]);
            }
            if self.tight {
                plot_string.push_str(&format!("plt.autoscale(tight=True)\n")[..]);
            }
            if let Some(t) = title {
                plot_string.push_str(&format!("plt.title(r\"{}\")\n", t)[..]);
            }
            if let Some(x) = xlabel {
                plot_string.push_str(&format!("plt.xlabel(r\"{}\")\n", x)[..]);
            }
            if let Some(y) = ylabel {
                plot_string.push_str(&format!("plt.ylabel(r\"{}\")\n", y)[..]);
            }
            match self.xscale {
                PlotScale::Linear => plot_string.push_str(&"plt.xscale(\"linear\")\n".to_string()[..]),
                PlotScale::Log => plot_string.push_str(&"plt.xscale(\"log\")\n".to_string()[..]),
            }
            match self.yscale {
                PlotScale::Linear => plot_string.push_str(&"plt.yscale(\"linear\")\n".to_string()[..]),
                PlotScale::Log => plot_string.push_str(&"plt.yscale(\"log\")\n".to_string()[..]),
            }
            if let Some(xl) = self.xlim {
                plot_string.push_str(&"plt.xlim(xl)\n".to_string()[..]);
            }
            if let Some(yl) = self.ylim {
                plot_string.push_str(&"plt.ylim(yl)\n".to_string()[..]);
            }

            for i in 0..y_length {
                let mut inner_string = format!("x,y[{}]", i);
                if !markers.is_empty() {
                    inner_string.push_str(&format!(",marker=\"{}\"", markers[i])[..]);
                }
                if !line_style.is_empty() {
                    inner_string.push_str(&format!(",linestyle=\"{}\"", line_style[i])[..]);
                }
                if !color.is_empty() {
                    inner_string.push_str(&format!(",color=\"{}\"", color[i])[..]);
                }
                if !legends.is_empty() {
                    inner_string.push_str(&format!(",label=r\"{}\"", legends[i])[..]);
                }
                if !alpha.is_empty() {
                    inner_string.push_str(&format!(",alpha={}", alpha[i])[..]);
                }
                plot_string.push_str(&format!("plt.plot({})\n", inner_string)[..]);
            }
            for i in 0..pair_length {
                let mut inner_string = format!("pair[{}][0],pair[{}][1]", i, i);
                if !markers.is_empty() {
                    inner_string.push_str(&format!(",marker=\"{}\"", markers[i + y_length])[..]);
                }
                if !line_style.is_empty() {
                    inner_string.push_str(&format!(",linestyle=\"{}\"", line_style[i + y_length])[..]);
                }
                if !color.is_empty() {
                    inner_string.push_str(&format!(",color=\"{}\"", color[i + y_length])[..]);
                }
                if !legends.is_empty() {
                    inner_string.push_str(&format!(",label=r\"{}\"", legends[i + y_length])[..]);
                }
                if !alpha.is_empty() {
                    inner_string.push_str(&format!(",alpha={}", alpha[i + y_length])[..]);
                }
                plot_string.push_str(&format!("plt.plot({})\n", inner_string)[..]);
            }

            if self.tight {
                plot_string.push_str(&format!("plt.legend()\nplt.savefig(pa, dpi={}, bbox_inches='tight')", dpi)[..]);
            } else {
                plot_string.push_str(&format!("plt.legend()\nplt.savefig(pa, dpi={})", dpi)[..]);
            }

            py.run(&plot_string[..], Some(&globals), None)?;
            Ok(())
        })
    }
}
