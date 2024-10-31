//! Plotting module for peroxide
//!
//! For Rust, there are some plot libraries but, still difficult to use.
//! Practically, using python is best choice to plot. And there is awesome crate - [pyo3](https://crates.io/crates/pyo3).
//!
//! # Prerequisite
//!
//! - python 3
//! - matplotlib
//! - scienceplots (Optional)
//!
//! # Usage
//!
//! To use this module, you should enable `plot` feature in `Cargo.toml`
//!
//! ```
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let x = linspace(0, 1, 100);
//!     let y1 = x.fmap(|t| t.powi(2));
//!     let y2 = x.fmap(|t| t.powi(3));
//!
//!     let mut rng = SmallRng::seed_from_u64(42);
//!     let normal = Normal(0f64, 0.1);
//!     let eps = normal.sample_with_rng(&mut rng, x.len());
//!     let y3 = y2.add_v(&eps);
//!
//!     let mut plt = Plot2D::new();
//!     plt.set_domain(x)
//!         .insert_image(y1)
//!         .insert_image(y2)
//!         .insert_image(y3)
//!         .set_legend(vec![r"$y=x^2$", r"$y=x^3$", r"$y=x^2 + \epsilon$"])
//!         .set_line_style(vec![(0, LineStyle::Dashed), (1, LineStyle::Dotted)])
//!         .set_plot_type(vec![(2, PlotType::Scatter)])
//!         .set_marker(vec![(2, Markers::Point)])
//!         .set_color(vec![(0, "red"), (1, "darkblue"), (2, "olive")])
//!         .set_xlabel(r"$x$")
//!         .set_ylabel(r"$y$")
//!         .set_style(PlotStyle::Nature) // if you want to use scienceplots
//!         .set_dpi(600)
//!         .tight_layout()
//!         .set_path("example_data/test_plot.png")
//!         .savefig().unwrap();
//! }
//! ```
//!
//! This code will generate below plot
//!
//! ![test_plot](https://github.com/Axect/Peroxide/blob/master/example_data/test_plot.png?raw=true)
//!
//! # Available Plot Options
//! - `set_domain` : Set x data
//! - `insert_image` : Insert y data
//! - `insert_pair` : Insert (x, y) data
//! - `set_title` : Set title of plot (optional)
//! - `set_xlabel` : Set xlabel of plot (optional)
//! - `set_ylabel` : Set ylabel of plot (optional)
//! - `set_zlabel` : Set zlabel of plot (optional; for 3D plot)
//! - `set_xscale` : Set xscale of plot (optional; `PlotScale::Linear` or `PlotScale::Log`)
//! - `set_yscale` : Set yscale of plot (optional; `PlotScale::Linear` or `PlotScale::Log`)
//! - `set_xlim` : Set xlim of plot (optional)
//! - `set_ylim` : Set ylim of plot (optional)
//! - `set_legend` : Set legend of plot (optional)
//! - `set_path` : Set path of plot (with filename - e.g. "example_data/test_plot.png")
//! - `set_fig_size` : Set figure size of plot (optional)
//! - `set_dpi` : Set dpi of plot (optional)
//! - `grid` : Set grid of plot (Grid::On, Grid::Off (default))
//! - `set_marker` : Set marker of plot (optional; `Markers::{Point, Line, Circle, TriangleUp, ...}`)
//! - `set_style` : Set style of plot (`PlotStyle::Nature`, `PlotStyle::IEEE`, `PlotStyle::Default` (default), `PlotStyle::Science`)
//! - `tight_layout` : Set tight layout of plot (optional)
//! - `set_line_style` : Set line style of plot (optional; `LineStyle::{Solid, Dashed, Dotted, DashDot}`)
//! - `set_color` : Set color of plot (optional; Vec<(usize, &str)>)
//! - `set_alpha` : Set alpha of plot (optional; Vec<(usize, f64)>)
//! - `set_plot_type` : Set plot type of plot (optional; `PlotType::{Scatter, Line, Bar}`)
//! - `savefig` : Save plot with given path

extern crate pyo3;
use self::pyo3::types::IntoPyDict;
use self::pyo3::{PyResult, Python};
pub use self::Grid::{Off, On};
use self::PlotOptions::{Domain, Images, Pairs, Path};
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

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum PlotType {
    Scatter,
    Line,
    Bar,
}

impl Display for PlotType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            PlotType::Scatter => "scatter".to_string(),
            PlotType::Line => "line".to_string(),
            PlotType::Bar => "bar".to_string(),
        };
        write!(f, "{}", str)
    }
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
    fn set_marker(&mut self, styles: Vec<(usize, Markers)>) -> &mut Self;
    fn set_style(&mut self, style: PlotStyle) -> &mut Self;
    fn tight_layout(&mut self) -> &mut Self;
    fn set_line_style(&mut self, style: Vec<(usize, LineStyle)>) -> &mut Self;
    fn set_color(&mut self, color: Vec<(usize, &str)>) -> &mut Self;
    fn set_alpha(&mut self, alpha: Vec<(usize, f64)>) -> &mut Self;
    fn set_plot_type(&mut self, plot_type: Vec<(usize, PlotType)>) -> &mut Self;
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
    markers: Vec<(usize, Markers)>,
    line_style: Vec<(usize, LineStyle)>,
    color: Vec<(usize, String)>,
    alpha: Vec<(usize, f64)>,
    path: String,
    fig_size: Option<(usize, usize)>,
    dpi: usize,
    grid: Grid,
    style: PlotStyle,
    tight: bool,
    plot_type: Vec<(usize, PlotType)>,
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
            plot_type: vec![],
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

    fn set_marker(&mut self, styles: Vec<(usize, Markers)>) -> &mut Self {
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

    fn set_line_style(&mut self, style: Vec<(usize, LineStyle)>) -> &mut Self {
        self.line_style = style;
        self
    }

    fn set_color(&mut self, color: Vec<(usize, &str)>) -> &mut Self {
        self.color = color.into_iter().map(|(i, x)| (i, x.to_owned())).collect();
        self
    }

    fn set_alpha(&mut self, alpha: Vec<(usize, f64)>) -> &mut Self {
        self.alpha = alpha;
        self
    }

    fn set_plot_type(&mut self, plot_type: Vec<(usize, PlotType)>) -> &mut Self {
        self.plot_type = plot_type;
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
            let markers = self
                .markers
                .iter()
                .map(|(i, x)| (i, format!("{}", x)))
                .collect::<Vec<_>>();
            let line_style = self
                .line_style
                .iter()
                .map(|(i, x)| (i, format!("{}", x)))
                .collect::<Vec<_>>();
            let color = self.color.clone();
            let alpha = self.alpha.clone();
            let plot_type = self.plot_type.clone();

            // Global variables to plot
            let globals =
                vec![("plt", py.import_bound("matplotlib.pyplot")?)].into_py_dict_bound(py);
            globals.as_gil_ref().set_item("x", x)?;
            globals.as_gil_ref().set_item("y", ys)?;
            globals.as_gil_ref().set_item("pair", pairs)?;
            globals.as_gil_ref().set_item("n", y_length)?;
            globals.as_gil_ref().set_item("p", pair_length)?;
            if let Some(fs) = fig_size {
                globals.as_gil_ref().set_item("fs", fs)?;
            }
            globals.as_gil_ref().set_item("dp", dpi)?;
            globals.as_gil_ref().set_item("gr", grid)?;
            globals.as_gil_ref().set_item("pa", path)?;
            if let Some(xl) = self.xlim {
                globals.as_gil_ref().set_item("xl", xl)?;
            }
            if let Some(yl) = self.ylim {
                globals.as_gil_ref().set_item("yl", yl)?;
            }

            // Plot Code
            let mut plot_string = match self.style {
                PlotStyle::Default => "\
                    plt.rc(\"text\", usetex=True)\n\
                    plt.rc(\"font\", family=\"serif\")\n"
                    .to_string(),
                PlotStyle::Science => "\
                    import scienceplots\n\
                    plt.style.use(\"science\")\n"
                    .to_string(),
                _ => format!(
                    "\
                    import scienceplots\n\
                    plt.style.use([\"science\", \"{}\"])\n",
                    style
                ),
            };
            if fig_size.is_some() {
                plot_string.push_str(&"plt.figure(figsize=fs, dpi=dp)\n".to_string()[..]);
            } else {
                plot_string.push_str(&"plt.figure()\n".to_string()[..]);
            }
            if self.tight {
                plot_string.push_str(&"plt.autoscale(tight=True)\n".to_string()[..]);
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
                PlotScale::Linear => {
                    plot_string.push_str(&"plt.xscale(\"linear\")\n".to_string()[..])
                }
                PlotScale::Log => plot_string.push_str(&"plt.xscale(\"log\")\n".to_string()[..]),
            }
            match self.yscale {
                PlotScale::Linear => {
                    plot_string.push_str(&"plt.yscale(\"linear\")\n".to_string()[..])
                }
                PlotScale::Log => plot_string.push_str(&"plt.yscale(\"log\")\n".to_string()[..]),
            }
            if self.xlim.is_some() {
                plot_string.push_str(&"plt.xlim(xl)\n".to_string()[..]);
            }
            if self.ylim.is_some() {
                plot_string.push_str(&"plt.ylim(yl)\n".to_string()[..]);
            }

            for i in 0..y_length {
                let mut inner_string = format!("x,y[{}]", i);
                let is_corresponding_marker =
                    !markers.is_empty() && (markers.iter().any(|(&j, _)| j == i));
                if is_corresponding_marker {
                    let marker = markers.iter().find(|(&j, _)| j == i).unwrap().1.as_str();
                    inner_string.push_str(&format!(",marker=\"{}\"", marker)[..]);
                }
                let is_corresponding_line_style =
                    !line_style.is_empty() && (line_style.iter().any(|(&j, _)| j == i));
                if is_corresponding_line_style {
                    let style = line_style.iter().find(|(&j, _)| j == i).unwrap().1.as_str();
                    inner_string.push_str(&format!(",linestyle=\"{}\"", style)[..]);
                }
                let is_corresponding_color =
                    !color.is_empty() && (color.iter().any(|(j, _)| j == &i));
                if is_corresponding_color {
                    let color = color.iter().find(|(j, _)| j == &i).unwrap().1.as_str();
                    inner_string.push_str(&format!(",color=\"{}\"", color)[..]);
                }
                if !legends.is_empty() {
                    inner_string.push_str(&format!(",label=r\"{}\"", legends[i])[..]);
                }
                let is_corresponding_alpha =
                    !alpha.is_empty() && (alpha.iter().any(|(j, _)| j == &i));
                if is_corresponding_alpha {
                    let alpha = alpha.iter().find(|(j, _)| j == &i).unwrap().1;
                    inner_string.push_str(&format!(",alpha={}", alpha)[..]);
                }
                let is_corresponding_plot_type =
                    !plot_type.is_empty() && (plot_type.iter().any(|(j, _)| j == &i));
                if is_corresponding_plot_type {
                    let plot_type = plot_type.iter().find(|(j, _)| j == &i).unwrap().1;
                    match plot_type {
                        PlotType::Scatter => {
                            plot_string.push_str(&format!("plt.scatter({})\n", inner_string)[..]);
                        }
                        PlotType::Line => {
                            plot_string.push_str(&format!("plt.plot({})\n", inner_string)[..]);
                        }
                        PlotType::Bar => {
                            plot_string.push_str(&format!("plt.bar({})\n", inner_string)[..]);
                        }
                    }
                } else {
                    plot_string.push_str(&format!("plt.plot({})\n", inner_string)[..]);
                }
            }
            for i in 0..pair_length {
                let mut inner_string = format!("pair[{}][0],pair[{}][1]", i, i);
                let is_corresponding_marker =
                    !markers.is_empty() && (markers.iter().any(|(&j, _)| j == (i + y_length)));
                if is_corresponding_marker {
                    let marker = markers
                        .iter()
                        .find(|(&j, _)| j == (i + y_length))
                        .unwrap()
                        .1
                        .as_str();
                    inner_string.push_str(&format!(",marker=\"{}\"", marker)[..]);
                }
                let is_corresponding_line_style = !line_style.is_empty()
                    && (line_style.iter().any(|(&j, _)| j == (i + y_length)));
                if is_corresponding_line_style {
                    let style = line_style
                        .iter()
                        .find(|(&j, _)| j == (i + y_length))
                        .unwrap()
                        .1
                        .as_str();
                    inner_string.push_str(&format!(",linestyle=\"{}\"", style)[..]);
                }
                let is_corresponding_color =
                    !color.is_empty() && (color.iter().any(|(j, _)| j == &(i + y_length)));
                if is_corresponding_color {
                    let color = color
                        .iter()
                        .find(|(j, _)| j == &(i + y_length))
                        .unwrap()
                        .1
                        .as_str();
                    inner_string.push_str(&format!(",color=\"{}\"", color)[..]);
                }
                if !legends.is_empty() {
                    inner_string.push_str(&format!(",label=r\"{}\"", legends[i + y_length])[..]);
                }
                let is_corresponding_alpha =
                    !alpha.is_empty() && (alpha.iter().any(|(j, _)| j == &(i + y_length)));
                if is_corresponding_alpha {
                    let alpha = alpha.iter().find(|(j, _)| j == &(i + y_length)).unwrap().1;
                    inner_string.push_str(&format!(",alpha={}", alpha)[..]);
                }
                let is_corresponding_plot_type =
                    !plot_type.is_empty() && (plot_type.iter().any(|(j, _)| j == &(i + y_length)));
                if is_corresponding_plot_type {
                    let plot_type = plot_type
                        .iter()
                        .find(|(j, _)| j == &(i + y_length))
                        .unwrap()
                        .1;
                    match plot_type {
                        PlotType::Scatter => {
                            plot_string.push_str(&format!("plt.scatter({})\n", inner_string)[..]);
                        }
                        PlotType::Line => {
                            plot_string.push_str(&format!("plt.plot({})\n", inner_string)[..]);
                        }
                        PlotType::Bar => {
                            plot_string.push_str(&format!("plt.bar({})\n", inner_string)[..]);
                        }
                    }
                } else {
                    plot_string.push_str(&format!("plt.plot({})\n", inner_string)[..]);
                }
            }

            if !legends.is_empty() {
                plot_string.push_str("plt.legend()\n");
            }

            if self.tight {
                plot_string
                    .push_str(&format!("plt.savefig(pa, dpi={}, bbox_inches='tight')", dpi)[..]);
            } else {
                plot_string.push_str(&format!("plt.savefig(pa, dpi={})", dpi)[..]);
            }

            py.run_bound(&plot_string[..], Some(&globals), None)?;
            Ok(())
        })
    }
}
