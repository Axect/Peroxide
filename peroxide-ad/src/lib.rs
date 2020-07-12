extern crate proc_macro;
use proc_macro::TokenStream;

/// Maximum order for Taylor mode AD
const N: usize = 10;

#[proc_macro]
pub fn ad_struct_def(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. N {
        let mut body = "pub d0: f64,\n".to_string();
        for j in 1 .. i+1 {
            body.push_str(&format!("pub d{}: f64,\n", j));
        }
        let one = format!("#[derive(Debug, Copy, Clone, PartialEq)]
        pub struct AD{} {{
            {}
        }}\n", i, body);
        total.push_str(&one);
    }
    total.parse().unwrap()
}

#[proc_macro]
pub fn ad_display(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. N {
        let mut body = "".to_string();
        for j in 1 .. i+1 {
            body.push_str(&format!("s.push_str(&format!(\"  d{}: {{}}\n\", self.d{}));", j, j));
        }
        let one = format!("impl std::fmt::Display for AD{} {{
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {{
                let mut s = format!(\"AD\n  d0: {{}}\n\", self.d0);
                {}
                write!(f, \"{{}}\", s)
            }}
        }}\n", i, body);
        total.push_str(&one);
    }
    total.parse().unwrap()
}

#[proc_macro]
pub fn ad_impl(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. N {
        let mut arg = "d0: f64, ".to_string();
        let mut body = "d0,\n".to_string();
        for j in 1 .. i {
            arg.push_str(&format!("d{}: f64, ", j));
            body.push_str(&format!("d{},\n", j));
        }
        arg.push_str(&format!("d{}: f64", i));
        body.push_str(&format!("d{}", i));

        let one = format!("impl AD{} {{
            pub fn new({}) -> Self {{
                Self {{
                    {}
                }}
            }}

            pub fn print(&self) {{
                println!(\"{{}}\", self);
            }}
        }}", i, arg, body);
        total.push_str(&one);
    }

    total.parse().unwrap()
}

#[proc_macro]
pub fn ad_impl_neg(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. N {
        let mut body = "-self.d0, ".to_string();
        for j in 1 .. i {
            body.push_str(&format!("-self.d{}, ", j));
        }
        body.push_str(&format!("-self.d{}", i));
        let one = format!("impl Neg for AD{} {{
            type Output = Self;
            fn neg(self) -> Self::Output {{
                AD{}::new({})
            }}
        }}\n", i, i, body);
        total.push_str(&one);
    }
    total.parse().unwrap()
}

#[proc_macro]
pub fn ad_impl_add(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. N {
        for j in 1 .. i {
            let mut body = "".to_string();
            for k in 0 .. j+1 {
                body.push_str(&format!("z.d{} += rhs.d{};\n", k, k));
            }
            let one = format!("impl Add<AD{}> for AD{} {{
                type Output = AD{};
                
                fn add(self, rhs: AD{}) -> Self::Output {{
                    let mut z = self.clone();
                    {}
                    z
                }}
            }}\n", j, i, i, j, body);
            total.push_str(&one);
        }
        for j in i .. N {
            let mut body = "".to_string();
            for k in 0 .. i+1 {
                body.push_str(&format!("z.d{} += self.d{};\n", k, k));
            }
            let one = format!("impl Add<AD{}> for AD{} {{
                type Output = AD{};

                fn add(self, rhs: AD{}) -> Self::Output {{
                    let mut z = rhs.clone();
                    {}
                    z
                }}
            }}\n", j, i, j, j, body);
            total.push_str(&one);
        }
    }
    total.parse().unwrap()
}

#[proc_macro]
pub fn ad_impl_sub(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. N {
        for j in 1 .. i {
            let mut body = "".to_string();
            for k in 0 .. j+1 {
                body.push_str(&format!("z.d{} -= rhs.d{};\n", k, k));
            }
            let one = format!("impl Sub<AD{}> for AD{} {{
                type Output = AD{};
                
                fn sub(self, rhs: AD{}) -> Self::Output {{
                    let mut z = self.clone();
                    {}
                    z
                }}
            }}\n", j, i, i, j, body);
            total.push_str(&one);
        }
        for j in i .. N {
            let mut body = "".to_string();
            for k in 0 .. i+1 {
                body.push_str(&format!("z.d{} += self.d{};\n", k, k));
            }
            let one = format!("impl Sub<AD{}> for AD{} {{
                type Output = AD{};

                fn sub(self, rhs: AD{}) -> Self::Output {{
                    let mut z = -rhs.clone();
                    {}
                    z
                }}
            }}\n", j, i, j, j, body);
            total.push_str(&one);
        }
    }
    total.parse().unwrap()
}

fn factorial(n: usize) -> usize {
    let mut p = 1usize;
    for i in 1..(n + 1) {
        p *= i;
    }
    p
}

#[allow(non_snake_case)]
fn P(n: usize, r: usize) -> usize {
    let mut p = 1usize;
    for i in 0..r {
        p *= n - i;
    }
    p
}

#[allow(non_snake_case)]
fn C(n: usize, r: usize) -> usize {
    if r > n / 2 {
        return C(n, n - r);
    }

    P(n, r) / factorial(r)
}
