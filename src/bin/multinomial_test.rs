extern crate peroxide;
use peroxide::*;

fn main() {
    let a = Multinomial::new(
        c!(1,2,3,4)
    );
    a.print();
    a.eval(&c!(-1, 2, 3, 1)).print();
}