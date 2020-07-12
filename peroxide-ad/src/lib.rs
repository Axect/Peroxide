extern crate proc_macro;
use proc_macro::TokenStream;

#[proc_macro]
pub fn ad_struct_def(_item: TokenStream) -> TokenStream {
    let mut total = "".to_string();
    for i in 1 .. 33 {
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
    for i in 1 .. 33 {
        let mut body = "".to_string();
        for j in 1 .. i+1 {
            body.push_str(&format!("s.push_str(&format!(\"  d{}: {{}}\", self.d{}));", j, j));
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
    for i in 1 .. 33 {
        let mut body1 = "d0: 0f64,\n".to_string();
        for j in 1 .. i+1 {
            body1.push_str(&format!("d{}: 0f64,\n", j));
        }

        let one = format!("impl AD{} {{
            pub fn new() -> Self {{
                Self {{
                    {}
                }}
            }}

            pub fn print(&self) {{
                println!(\"{{}}\", self);
            }}
        }}", i, body1);
        total.push_str(&one);
    }

    total.parse().unwrap()
}
