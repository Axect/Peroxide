use proc_macro::TokenStream;
use quote::{quote, format_ident};
use syn::{Ident, Type, Token, parenthesized, Block}; 
use syn::parse::{Parse, ParseStream};
use syn::Result;

type Function = Ident;
type Variable = Ident;
type Input = Type;
type Return = Type;

struct Signature {
    function: Function,
    var: Variable,
    var_type: Input,
    return_t: Return,
    block: Block,
}

struct Closure {
    var: Variable,
    block: Block,
}

struct Syntax {
    _fn_token: Token!(fn),
    ident: Function,
    _paren_token: syn::token::Paren,
    paren_var: Ident,
    _colon_token: syn::token::Colon,
    paren_type: Type,
    _rarrow_token: Token!(->),
    return_t: Return,
    block: Block,
}

struct ClosureSyntax {
    _lb_token: Token!(|),
    var: Ident,
    _colon_token: syn::token::Colon,
    _bar_type: Type,
    _rb_token: Token!(|),
    block: Block,
}

impl Parse for Signature {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.is_empty() {
            panic!("Empty signature!");
        }

        let content1;
        let _fn_token = input.parse().unwrap();
        let ident = input.parse().unwrap();
        let _paren_token = parenthesized!(content1 in input);
        let paren_var = content1.parse().unwrap();
        let _colon_token = content1.parse().unwrap();
        let paren_type = content1.parse().unwrap();
        let _rarrow_token = input.parse().unwrap();
        let return_t = input.parse().unwrap();
        let block = input.parse().unwrap();
        let syntax = Syntax {
            _fn_token,
            ident,
            _paren_token,
            paren_var,
            _colon_token,
            paren_type,
            _rarrow_token,
            return_t,
            block,
        };

        Ok(Signature {
            function: syntax.ident,
            var: syntax.paren_var,
            var_type: syntax.paren_type,
            return_t: syntax.return_t,
            block: syntax.block,
        })
    }
}

impl Parse for Closure {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.is_empty() {
            panic!("Empty signature!");
        }

        let _lb_token = input.parse()?;
        let var = input.parse()?;
        let _colon_token = input.parse()?;
        let _bar_type = input.parse()?;
        let _rb_token = input.parse()?;
        let block = input.parse()?;

        let syntax = ClosureSyntax {
            _lb_token,
            var,
            _colon_token,
            _bar_type,
            _rb_token,
            block,
        };

        Ok(Closure {
            var: syntax.var,
            block: syntax.block,
        })
    }
}

#[proc_macro_attribute]
pub fn ad_function(_attr: TokenStream, input: TokenStream) -> TokenStream {
    let signature = syn::parse_macro_input!(input as Signature);
    let function = signature.function;
    let var = signature.var;
    let var_type = signature.var_type;
    let return_t = signature.return_t;
    let block = signature.block;

    let ad_fn = format_ident!("{}_ad", function);
    let grad_fn = format_ident!("{}_grad", function);
    let hess_fn = format_ident!("{}_hess", function);

    let tokens = quote!{
        fn #function(#var: #var_type) -> #return_t {
            #block
        }

        fn #ad_fn(#var: AD) -> AD {
            #block
        }

        fn #grad_fn(#var: f64) -> f64 {
            let ad = #ad_fn(AD1(#var, 1f64));
            ad.dx()
        }

        fn #hess_fn(#var: f64) -> f64 {
            let ad = #ad_fn(AD2(#var, 1f64, 0f64));
            ad.ddx()
        }
    };
    tokens.into()
}

#[proc_macro]
pub fn ad_closure(input: TokenStream) -> TokenStream {
    let closure = syn::parse_macro_input!(input as Closure);
    let var = closure.var;
    let block = closure.block;

    let tokens = quote!{
        |#var: AD| #block
    };
    tokens.into()
}

///// Maximum order for Taylor mode AD
//const N: usize = 5;
//
//#[proc_macro]
//pub fn ad_struct_def(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut body = "pub d0: f64,\n".to_string();
//        for j in 1 .. i+1 {
//            body.push_str(&format!("pub d{}: f64,\n", j));
//        }
//        let one = format!("#[derive(Debug, Copy, Clone, PartialEq, Default)]
//        pub struct AD{} {{
//            {}
//        }}\n", i, body);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_display(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut body = "".to_string();
//        for j in 1 .. i+1 {
//            body.push_str(&format!("s.push_str(&format!(\"  d{}: {{}}\n\", self.d{}));", j, j));
//        }
//        let one = format!("impl std::fmt::Display for AD{} {{
//            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {{
//                let mut s = format!(\"AD\n  d0: {{}}\n\", self.d0);
//                {}
//                write!(f, \"{{}}\", s)
//            }}
//        }}\n", i, body);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut arg = "d0: f64, ".to_string();
//        let mut body = "d0,\n".to_string();
//        for j in 1 .. i {
//            arg.push_str(&format!("d{}: f64, ", j));
//            body.push_str(&format!("d{},\n", j));
//        }
//        arg.push_str(&format!("d{}: f64", i));
//        body.push_str(&format!("d{}", i));
//
//        let one = format!("impl AD{} {{
//            pub fn new({}) -> Self {{
//                Self {{
//                    {}
//                }}
//            }}
//
//            pub fn print(&self) {{
//                println!(\"{{}}\", self);
//            }}
//
//            pub fn iter(&self) -> ADIter{} {{
//                self.into_iter()
//            }}
//
//            pub fn iter_mut(&mut self) -> ADIterMut{} {{
//                self.into_iter()
//            }}
//
//            pub fn len(&self) -> usize {{
//                {}
//            }}
//        }}", i, arg, body, i, i, i+1);
//        total.push_str(&one);
//    }
//
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_from(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        for j in 1 .. i {
//            let mut body = "high.d0, ".to_string();
//            for k in 1 .. j {
//                body.push_str(&format!("high.d{}, ", k));
//            }
//            body.push_str(&format!("high.d{}", j));
//            // i -> j (i > j)
//            let one = format!("impl From<AD{}> for AD{} {{
//                fn from(high: AD{}) -> Self {{
//                    Self::new({})
//                }}
//            }}\n", i, j, i, body);
//            let two = format!("impl<'a> From<&'a AD{}> for AD{} {{
//                fn from(high: &'a AD{}) -> Self {{
//                    Self::new({})
//                }}
//            }}\n", i, j, i, body);
//            total.push_str(&one);
//            total.push_str(&two);
//        }
//        // i vs i (&AD -> AD)
//        {
//            let mut body = "high.d0, ".to_string();
//            for k in 1 .. i {
//                body.push_str(&format!("high.d{}, ", k));
//            }
//            body.push_str(&format!("high.d{}", i));
//            let two = format!("impl<'a> From<&'a AD{}> for AD{} {{
//                fn from(high: &'a AD{}) -> Self {{
//                    Self::new({})
//                }}
//            }}\n", i, i, i, body);
//            total.push_str(&two);
//        }
//        for j in i+1 .. N+1 {
//            let mut body = "low.d0, ".to_string();
//            for k in 1 .. i+1 {
//                body.push_str(&format!("low.d{}, ", k));
//            }
//            for _k in i+1 .. j {
//                body.push_str("0f64, ");
//            }
//            body.push_str("0f64");
//            // i -> j (i <= j)
//            let one = format!("impl From<AD{}> for AD{} {{
//                fn from(low: AD{}) -> Self {{
//                    Self::new({})
//                }}
//            }}\n", i, j, i, body);
//            let two = format!("impl<'a> From<&'a AD{}> for AD{} {{
//                fn from(low: &'a AD{}) -> Self {{
//                    Self::new({})
//                }}
//            }}\n", i, j, i, body);
//            total.push_str(&one);
//            total.push_str(&two);
//        }
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_iter_def(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("#[derive(Debug)]
//        pub struct ADIntoIter{} {{
//            ad: AD{},
//            index: usize,
//            r_index: usize,
//        }}\n", i, i);
//        let two = format!("#[derive(Debug)]
//        pub struct ADIter{}<'a> {{
//            ad: &'a AD{},
//            index: usize,
//            r_index: usize,
//        }}\n", i, i);
//        let three = format!("#[derive(Debug)]
//        pub struct ADIterMut{}<'a> {{
//            ad: &'a mut AD{},
//            index: usize,
//            r_index: usize,
//        }}\n", i, i);
//        total.push_str(&one);
//        total.push_str(&two);
//        total.push_str(&three);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_into_iter(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl IntoIterator for AD{} {{
//            type Item = f64;
//            type IntoIter = ADIntoIter{};
//
//            fn into_iter(self) -> Self::IntoIter {{
//                ADIntoIter{} {{
//                    ad: self,
//                    index: 0,
//                    r_index: 0,
//                }}
//            }}
//        }}\n", i, i, i);
//        let two = format!("impl<'a> IntoIterator for &'a AD{} {{
//            type Item = f64;
//            type IntoIter = ADIter{}<'a>;
//
//            fn into_iter(self) -> Self::IntoIter {{
//                ADIter{} {{
//                    ad: self,
//                    index: 0,
//                    r_index: 0,
//                }}
//            }}
//        }}\n", i, i, i);
//        let three = format!("impl<'a> IntoIterator for &'a mut AD{} {{
//            type Item = f64;
//            type IntoIter = ADIterMut{}<'a>;
//
//            fn into_iter(self) -> Self::IntoIter {{
//                ADIterMut{} {{
//                    ad: self,
//                    index: 0,
//                    r_index: 0,
//                }}
//            }}
//        }}\n", i, i, i);
//        total.push_str(&one);
//        total.push_str(&two);
//        total.push_str(&three);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_iter(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut body = "0 => self.ad.d0,\n".to_string();
//        for j in 1 .. i+1 {
//            body.push_str(&format!("{} => self.ad.d{},\n", j, j));
//        }
//        body.push_str("_ => return None,");
//        let one = format!("impl Iterator for ADIntoIter{} {{
//            type Item = f64;
//            fn next(&mut self) -> Option<Self::Item> {{
//                if self.index + self.r_index < {} {{
//                    let result = match self.index {{
//                        {}
//                    }};
//                    self.index += 1;
//                    Some(result)
//                }} else {{
//                    None
//                }}
//            }}
//
//            fn size_hint(&self) -> (usize, Option<usize>) {{
//                let lower = self.ad.len() - (self.index + self.r_index);
//                let upper = self.ad.len() - (self.index + self.r_index);
//                (lower, Some(upper))
//            }}
//        }}\n", i, i+1, body);
//        let two = format!("impl<'a> Iterator for ADIter{}<'a> {{
//            type Item = f64;
//            fn next(&mut self) -> Option<Self::Item> {{
//                if self.index + self.r_index < {} {{
//                    let result = match self.index {{
//                        {}
//                    }};
//                    self.index += 1;
//                    Some(result)
//                }} else {{
//                    None
//                }}
//            }}
//            fn size_hint(&self) -> (usize, Option<usize>) {{
//                let lower = self.ad.len() - (self.index + self.r_index);
//                let upper = self.ad.len() - (self.index + self.r_index);
//                (lower, Some(upper))
//            }}
//        }}\n", i, i+1, body);
//        let three = format!("impl<'a> Iterator for ADIterMut{}<'a> {{
//            type Item = f64;
//            fn next(&mut self) -> Option<Self::Item> {{
//                if self.index + self.r_index < {} {{
//                    let result = match self.index {{
//                        {}
//                    }};
//                    self.index += 1;
//                    Some(result)
//                }} else {{
//                    None
//                }}
//            }}
//            fn size_hint(&self) -> (usize, Option<usize>) {{
//                let lower = self.ad.len() - (self.index + self.r_index);
//                let upper = self.ad.len() - (self.index + self.r_index);
//                (lower, Some(upper))
//            }}
//        }}\n", i, i+1, body);
//        total.push_str(&one);
//        total.push_str(&two);
//        total.push_str(&three);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_from_iter(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl FromIterator<f64> for AD{} {{
//            fn from_iter<I: IntoIterator<Item=f64>>(iter: I) -> Self {{
//                let mut z = Self::default();
//                for (i, elem) in iter.into_iter().enumerate() {{
//                    z[i] = elem;
//                }}
//                z
//            }}
//        }}\n", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_double_ended_iter(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut body = format!("0 => self.ad.d{},\n", i);
//        for j in 1 .. i+1 {
//            body.push_str(&format!("{} => self.ad.d{},\n", j, i-j));
//        }
//        body.push_str("_ => return None,");
//        let one = format!("impl DoubleEndedIterator for ADIntoIter{} {{
//            fn next_back(&mut self) -> Option<Self::Item> {{
//                if self.index + self.r_index == {} {{
//                    return None;
//                }}
//                let result = match self.r_index {{
//                    {}
//                }};
//                self.r_index += 1;
//                Some(result)
//            }}
//        }}\n", i, i+1, body);
//        let two = format!("impl<'a> DoubleEndedIterator for ADIter{}<'a> {{
//            fn next_back(&mut self) -> Option<Self::Item> {{
//                if self.index + self.r_index == {} {{
//                    return None;
//                }}
//                let result = match self.r_index {{
//                    {}
//                }};
//                self.r_index += 1;
//                Some(result)
//            }}
//        }}\n", i, i+1, body);
//        let three = format!("impl<'a> DoubleEndedIterator for ADIterMut{}<'a> {{
//            fn next_back(&mut self) -> Option<Self::Item> {{
//                if self.index + self.r_index == {} {{
//                    return None;
//                }}
//                let result = match self.r_index {{
//                    {}
//                }};
//                self.r_index += 1;
//                Some(result)
//            }}
//        }}\n", i, i+1, body);
//        total.push_str(&one);
//        total.push_str(&two);
//        total.push_str(&three);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_exact_size_iter(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl ExactSizeIterator for ADIntoIter{} {{
//            fn len(&self) -> usize {{
//                self.ad.len() - (self.index + self.r_index)
//            }}
//        }}\n", i);
//        let two = format!("impl<'a> ExactSizeIterator for ADIter{}<'a> {{
//            fn len(&self) -> usize {{
//                self.ad.len() - (self.index + self.r_index)
//            }}
//        }}\n", i);
//        let three = format!("impl<'a> ExactSizeIterator for ADIterMut{}<'a> {{
//            fn len(&self) -> usize {{
//                self.ad.len() - (self.index + self.r_index)
//            }}
//        }}\n", i);
//        total.push_str(&one);
//        total.push_str(&two);
//        total.push_str(&three);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_index(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut body1 = "0 => &self.d0,\n".to_string();
//        let mut body2 = "0 => &mut self.d0,\n".to_string();
//        for j in 1 .. i+1 {
//            body1.push_str(&format!("{} => &self.d{},\n", j, j));
//            body2.push_str(&format!("{} => &mut self.d{},\n", j, j));
//        }
//        let one = format!("impl Index<usize> for AD{} {{
//            type Output = f64;
//
//            fn index(&self, n: usize) -> &Self::Output {{
//                match n {{
//                    {}
//                    _ => panic!(\"{{}} exceed order of AD{}\", n),
//                }}
//            }}
//        }}\n", i, body1, i);
//        let two = format!("impl IndexMut<usize> for AD{} {{
//            fn index_mut(&mut self, n: usize) -> &mut Self::Output {{
//                match n {{
//                    {}
//                    _ => panic!(\"{{}} exceed order of AD{}\", n),
//                }}
//            }}
//        }}\n", i, body2, i);
//        total.push_str(&one);
//        total.push_str(&two);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_neg(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let mut body = "-self.d0, ".to_string();
//        for j in 1 .. i {
//            body.push_str(&format!("-self.d{}, ", j));
//        }
//        body.push_str(&format!("-self.d{}", i));
//        let one = format!("impl Neg for AD{} {{
//            type Output = Self;
//            fn neg(self) -> Self::Output {{
//                AD{}::new({})
//            }}
//        }}\n", i, i, body);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_add(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        for j in 1 .. i {
//            let mut body = "".to_string();
//            for k in 0 .. j+1 {
//                body.push_str(&format!("z.d{} += rhs.d{};\n", k, k));
//            }
//            let one = format!("impl Add<AD{}> for AD{} {{
//                type Output = AD{};
//                
//                fn add(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = self.clone();
//                    {}
//                    z
//                }}
//            }}\n", j, i, i, j, body);
//            total.push_str(&one);
//        }
//        for j in i .. N+1 {
//            let mut body = "".to_string();
//            for k in 0 .. i+1 {
//                body.push_str(&format!("z.d{} += self.d{};\n", k, k));
//            }
//            let one = format!("impl Add<AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn add(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = rhs.clone();
//                    {}
//                    z
//                }}
//            }}\n", j, i, j, j, body);
//            total.push_str(&one);
//        }
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_sub(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        for j in 1 .. i {
//            let mut body = "".to_string();
//            for k in 0 .. j+1 {
//                body.push_str(&format!("z.d{} -= rhs.d{};\n", k, k));
//            }
//            let one = format!("impl Sub<AD{}> for AD{} {{
//                type Output = AD{};
//                
//                fn sub(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = self.clone();
//                    {}
//                    z
//                }}
//            }}\n", j, i, i, j, body);
//            total.push_str(&one);
//        }
//        for j in i .. N+1 {
//            let mut body = "".to_string();
//            for k in 0 .. i+1 {
//                body.push_str(&format!("z.d{} += self.d{};\n", k, k));
//            }
//            let one = format!("impl Sub<AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn sub(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = -rhs.clone();
//                    {}
//                    z
//                }}
//            }}\n", j, i, j, j, body);
//            total.push_str(&one);
//        }
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_mul(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        for j in 1 .. i+1 {
//            let one = format!("impl Mul<AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn mul(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = self.clone();
//                    let y = Self::Output::from(rhs);
//                    for t in 0 .. z.len() {{
//                        z[t] = self.iter()
//                                .take(t + 1)
//                                .zip(y.iter().take(t+1).rev())
//                                .enumerate()
//                                .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
//                    }}
//                    z
//                }}
//            }}\n", j, i, i, j);
//            let two = format!("impl<'a> Mul<&'a AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn mul(self, rhs: &'a AD{}) -> Self::Output {{
//                    let mut z = self.clone();
//                    let y = Self::Output::from(rhs);
//                    for t in 0 .. z.len() {{
//                        z[t] = self.iter()
//                                .take(t + 1)
//                                .zip(y.iter().take(t+1).rev())
//                                .enumerate()
//                                .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
//                    }}
//                    z
//                }}
//            }}\n", j, i, i, j);
//            total.push_str(&one);
//            total.push_str(&two);
//        }
//        for j in i+1 .. N+1 {
//            let one = format!("impl Mul<AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn mul(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = rhs.clone();
//                    let x = Self::Output::from(self);
//                    for t in 0 .. z.len() {{
//                        z[t] = x.iter()
//                                .take(t + 1)
//                                .zip(rhs.iter().take(t+1).rev())
//                                .enumerate()
//                                .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
//                    }}
//                    z
//                }}
//            }}\n", j, i, j, j);
//            let two = format!("impl<'a> Mul<&'a AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn mul(self, rhs: &'a AD{}) -> Self::Output {{
//                    let mut z = rhs.clone();
//                    let x = Self::Output::from(self);
//                    let y = Self::Output::from(rhs);
//                    for t in 0 .. z.len() {{
//                        z[t] = x.iter()
//                                .take(t + 1)
//                                .zip(y.iter().take(t+1).rev())
//                                .enumerate()
//                                .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
//                    }}
//                    z
//                }}
//            }}\n", j, i, j, j);
//            total.push_str(&one);
//            total.push_str(&two);
//        }
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_div(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        // i >= j (self >= rhs)
//        for j in 1 .. i+1 {
//            let one = format!("impl Div<AD{}> for AD{} {{
//                type Output = AD{};
//                
//                fn div(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = Self::Output::default();
//                    let y = Self::Output::from(rhs);
//                    z[0] = self[0] / y[0];
//                    let y0 = 1f64 / y[0];
//                    for i in 1 .. z.len() {{
//                        let mut s = 0f64;
//                        for (j, (y1, z1)) in y.iter().skip(1).take(i).zip(z.iter().take(i).rev()).enumerate() {{
//                            s += (C(i, j+1) as f64) * y1 * z1;
//                        }}
//                        z[i] = y0 * (self[i] - s);
//                    }}
//                    z
//                }}
//            }}", j, i, i, j);
//            total.push_str(&one);
//        }
//        // i < j (self < rhs)
//        for j in i+1 .. N+1 {
//            let one = format!("impl Div<AD{}> for AD{} {{
//                type Output = AD{};
//
//                fn div(self, rhs: AD{}) -> Self::Output {{
//                    let mut z = Self::Output::default();
//                    let x = Self::Output::from(self);
//                    z[0] = x[0] / rhs[0];
//                    let y0 = 1f64 / rhs[0];
//                    for i in 1 .. z.len() {{
//                        let mut s = 0f64;
//                        for (j, (y1, z1)) in rhs.iter().skip(1).take(i).zip(z.iter().take(i).rev()).enumerate() {{
//                            s += (C(i, j+1) as f64) * y1 * z1;
//                        }}
//                        z[i] = y0 * (x[i] - s);
//                    }}
//                    z
//                }}
//            }}", j, i, j, j);
//            total.push_str(&one);
//        }
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_explogops(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl ExpLogOps for AD{} {{
//            fn exp(&self) -> Self {{
//                let mut z = Self::default();
//                z[0] = self[0].exp();
//                for i in 1 .. z.len() {{
//                    z[i] = z.iter()
//                        .take(i)
//                        .zip(self.iter().skip(1).take(i).rev())
//                        .enumerate()
//                        .fold(0f64, |x, (k, (z1, x1))| x + (C(i-1, k) as f64) * x1 * z1);
//                }}
//                z
//            }}
//
//            fn ln(&self) -> Self {{
//                let mut z = Self::default();
//                z[0] = self[0].ln();
//                let x0 = 1f64 / self[0];
//                for i in 1 .. z.len() {{
//                    let mut s = 0f64;
//                    for (k, (z1, x1)) in z.iter().skip(1).take(i-1).zip(self.iter().skip(1).take(i-1).rev()).enumerate() {{
//                        s += (C(i-1, k+1) as f64) * z1 * x1;
//                    }}
//                    z[i] = x0 * (self[i] - s);
//                }}
//                z
//            }}
//
//            fn log(&self, base: f64) -> Self {{
//                self.ln().iter().map(|x| x / base.ln()).collect()
//            }}
//        }}", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_powops(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl PowOps for AD{} {{
//            fn powi(&self, n: i32) -> Self {{
//                let mut z = self.clone();
//                for _i in 1 .. n {{
//                    z = z * self;
//                }}
//                z
//            }}
//
//            fn powf(&self, f: f64) -> Self {{
//                let ln_x = self.ln();
//                let mut z = Self::default();
//                z[0] = self.d0.powf(f);
//                for i in 1 .. z.len() {{
//                    let mut s = 0f64;
//                    for (j, (z1, ln_x1)) in z.iter().skip(1).take(i-1).zip(ln_x.iter().skip(1).take(i-1).rev()).enumerate() {{
//                        s += (C(i-1, j+1) as f64) * z1 * ln_x1;
//                    }}
//                    z[i] = f * (z[0] * ln_x[i] + s);
//                }}
//                z
//            }}
//
//            fn pow(&self, _f: Self) -> Self {{
//                unimplemented!()
//            }}
//        }}", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_trigops(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl TrigOps for AD{} {{
//            fn sin_cos(&self) -> (Self, Self) {{
//                let mut u = Self::default();
//                let mut v = Self::default();
//                u[0] = self[0].sin();
//                v[0] = self[0].cos();
//                for i in 1 .. u.len() {{
//                    u[i] = v.iter()
//                        .take(i)
//                        .zip(self.iter().skip(1).take(i).rev())
//                        .enumerate()
//                        .fold(0f64, |x, (k, (v1, x1))| x + (C(i-1, k) as f64) * x1 * v1);
//                    v[i] = -u.iter()
//                        .take(i)
//                        .zip(self.iter().skip(1).take(i).rev())
//                        .enumerate()
//                        .fold(0f64, |x, (k, (u1, x1))| x + (C(i-1, k) as f64) * x1 * u1);
//                }}
//                (u, v)
//            }}
//
//            fn sinh_cosh(&self) -> (Self, Self) {{
//                let mut u = Self::default();
//                let mut v = Self::default();
//                u[0] = self[0].sinh();
//                v[0] = self[0].cosh();
//                for i in 1 .. u.len() {{
//                    u[i] = v.iter()
//                        .take(i)
//                        .zip(self.iter().skip(1).take(i).rev())
//                        .enumerate()
//                        .fold(0f64, |x, (k, (v1, x1))| x + (C(i-1, k) as f64) * x1 * v1);
//                    v[i] = u.iter()
//                        .take(i)
//                        .zip(self.iter().skip(1).take(i).rev())
//                        .enumerate()
//                        .fold(0f64, |x, (k, (u1, x1))| x + (C(i-1, k) as f64) * x1 * u1);
//                }}
//                (u, v)
//                
//            }}
//
//            fn asin(&self) -> Self {{
//                unimplemented!()
//            }}
//
//            fn acos(&self) -> Self {{
//                unimplemented!()
//            }}
//            
//            fn atan(&self) -> Self {{
//                unimplemented!()
//            }}
//
//            fn asinh(&self) -> Self {{
//                unimplemented!()
//            }}
//
//            fn acosh(&self) -> Self {{
//                unimplemented!()
//            }}
//
//            fn atanh(&self) -> Self {{
//                unimplemented!()
//            }}
//        }}\n", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn def_ad(_item: TokenStream) -> TokenStream {
//    let mut from_into = "".to_string();
//    let mut to_ad = "".to_string();
//    for i in 1 .. N+1 {
//        from_into.push_str(&format!("+ From<AD{}>\n", i));
//        from_into.push_str(&format!("+ Into<AD{}>\n", i));
//        to_ad.push_str(&format!("fn to_ad{}(self) -> AD{} {{ self.into() }}\n", i, i));
//    }
//    let total = format!("pub trait AD: 
//        std::fmt::Display
//        + Clone
//        + Copy
//        + PartialEq
//        + From<f64>
//        + Into<f64>
//        + Add<Output = Self>
//        + Sub<Output = Self>
//        + Mul<Output = Self>
//        + Div<Output = Self>
//        + Add<f64, Output = Self>
//        + Sub<f64, Output = Self>
//        + Mul<f64, Output = Self>
//        + Div<f64, Output = Self>
//        + PowOps
//        + ExpLogOps
//        + TrigOps
//        {}
//    {{
//         {}
//    }}", from_into, to_ad);
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_ad(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl AD for AD{} {{}}", i);
//        total.push_str(&one);
//    } 
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_from_type(item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl From<{}> for AD{} {{
//            fn from(other: {}) -> Self {{
//                let f = other as f64;
//                let mut z = Self::default();
//                z.d0 = f;
//                z
//            }}
//        }}\n", item, i, item);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_add_f64(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Add<f64> for AD{} {{
//            type Output = Self;
//
//            fn add(self, rhs: f64) -> Self::Output {{
//                let mut z = self;
//                z.d0 += rhs;
//                z
//            }}
//        }}\n", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_sub_f64(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Sub<f64> for AD{} {{
//            type Output = Self;
//
//            fn sub(self, rhs: f64) -> Self::Output {{
//                let mut z = self;
//                z.d0 -= rhs;
//                z
//            }}
//        }}\n", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_mul_f64(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Mul<f64> for AD{} {{
//            type Output = Self;
//
//            fn mul(self, rhs: f64) -> Self::Output {{
//                self.iter().map(|x| x * rhs).collect()
//            }}
//        }}\n", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_div_f64(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Div<f64> for AD{} {{
//            type Output = Self;
//
//            fn div(self, rhs: f64) -> Self::Output {{
//                self.iter().map(|x| x / rhs).collect()
//            }}
//        }}\n", i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn f64_impl_add_ad(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Add<AD{}> for f64 {{
//            type Output = AD{};
//
//            fn add(self, rhs: AD{}) -> Self::Output {{
//                let mut z = rhs;
//                z.d0 += self;
//                z
//            }}
//        }}\n", i, i, i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn f64_impl_sub_ad(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Sub<AD{}> for f64 {{
//            type Output = AD{};
//
//            fn sub(self, rhs: AD{}) -> Self::Output {{
//                let mut z = -rhs;
//                z.d0 += self;
//                z
//            }}
//        }}\n", i, i, i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn f64_impl_mul_ad(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Mul<AD{}> for f64 {{
//            type Output = AD{};
//
//            fn mul(self, rhs: AD{}) -> Self::Output {{
//                rhs.iter().map(|x| x * self).collect()
//            }}
//        }}\n", i, i, i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn f64_impl_div_ad(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl Div<AD{}> for f64 {{
//            type Output = AD{};
//
//            fn div(self, rhs: AD{}) -> Self::Output {{
//                let ad1 = AD1::from(self);
//                ad1 / rhs
//            }}
//        }}\n", i, i, i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//
//#[proc_macro]
//pub fn f64_impl_from_ad(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 1 .. N+1 {
//        let one = format!("impl From<AD{}> for f64 {{
//            fn from(ad: AD{}) -> Self {{
//                ad.d0
//            }}
//        }}\n", i, i);
//        total.push_str(&one);
//    }
//    total.parse().unwrap()
//}
//
//#[proc_macro]
//pub fn ad_impl_stable_fn(_item: TokenStream) -> TokenStream {
//    let mut total = "".to_string();
//    for i in 2 .. N+1 {
//        for j in 2 .. (i+1) {
//            let one = format!("impl<F:Fn(AD{}) -> AD{}> StableFn<AD{}> for ADLift<F, AD{}> {{
//                type Output = AD{};
//                fn call_stable(&self, target: Self::Output) -> Self::Output {{
//                    self.f(AD{}::from(target)).into()
//                }}
//            }}\n", i, i, j, i, j, i);
//            total.push_str(&one);
//        }
//    }
//    total.parse().unwrap()
//}
