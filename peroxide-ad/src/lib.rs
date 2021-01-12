use proc_macro::TokenStream;
use watt::WasmMacro;

static MACRO: WasmMacro = WasmMacro::new(WASM);
static WASM: &[u8] = include_bytes!("peroxide_ad.wasm");

#[proc_macro]
pub fn ad_struct_def(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_struct_def", item)
}

#[proc_macro]
pub fn ad_display(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_display", item)
}

#[proc_macro]
pub fn ad_impl(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl", item)
}

#[proc_macro]
pub fn ad_impl_from(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_from", item)
}

#[proc_macro]
pub fn ad_iter_def(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_iter_def", item)
}

#[proc_macro]
pub fn ad_impl_into_iter(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_into_iter", item)
}

#[proc_macro]
pub fn ad_impl_iter(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_iter", item)
}

#[proc_macro]
pub fn ad_impl_from_iter(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_from_iter", item)
}

#[proc_macro]
pub fn ad_impl_double_ended_iter(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_double_ended_iter", item)
}

#[proc_macro]
pub fn ad_impl_exact_size_iter(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_exact_size_iter", item)
}

#[proc_macro]
pub fn ad_impl_index(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_index", item)
}

#[proc_macro]
pub fn ad_impl_neg(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_neg", item)
}

#[proc_macro]
pub fn ad_impl_add(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_add", item)
}

#[proc_macro]
pub fn ad_impl_sub(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_sub", item)
}

#[proc_macro]
pub fn ad_impl_mul(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_mul", item)
}

#[proc_macro]
pub fn ad_impl_div(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_div", item)
}

#[proc_macro]
pub fn ad_impl_explogops(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_explogops", item)
}

#[proc_macro]
pub fn ad_impl_powops(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_powops", item)
}

#[proc_macro]
pub fn ad_impl_trigops(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_trigops", item)
}

#[proc_macro]
pub fn def_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("def_ad", item)
}

#[proc_macro]
pub fn ad_impl_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_ad", item)
}

#[proc_macro]
pub fn ad_impl_from_type(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_from_type", item)
}

#[proc_macro]
pub fn ad_impl_add_f64(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_add_f64", item)
}

#[proc_macro]
pub fn ad_impl_sub_f64(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_sub_f64", item)
}

#[proc_macro]
pub fn ad_impl_mul_f64(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_mul_f64", item)
}

#[proc_macro]
pub fn ad_impl_div_f64(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_div_f64", item)
}

#[proc_macro]
pub fn f64_impl_add_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("f64_impl_add_ad", item)
}

#[proc_macro]
pub fn f64_impl_sub_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("f64_impl_sub_ad", item)
}

#[proc_macro]
pub fn f64_impl_mul_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("f64_impl_mul_ad", item)
}

#[proc_macro]
pub fn f64_impl_div_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("f64_impl_div_ad", item)
}


#[proc_macro]
pub fn f64_impl_from_ad(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("f64_impl_from_ad", item)
}

#[proc_macro]
pub fn ad_impl_stable_fn(item: TokenStream) -> TokenStream {
    MACRO.proc_macro("ad_impl_stable_fn", item)
}
