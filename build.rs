use std::env;

macro_rules! feature(($name:expr) => (env::var(concat!("CARGO_FEATURE_", $name)).is_ok()));

fn main() {
    if feature!("NATIVE") {
        println!("cargo:rustc-link-search=/opt/OpenBLAS/lib");
        println!("cargo:rustc-link-lib=dylib=gfortran");
        println!("cargo:rustc-link-lib=dylib=openblas");
    } else {

    }
}