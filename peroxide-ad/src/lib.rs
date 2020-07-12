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
