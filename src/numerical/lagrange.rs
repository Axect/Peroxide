#[allow(unused_imports)]
use structure::polynomial::*;
use structure::vector::*;

pub fn lagrange_polynomial(node_x: Vector, node_y: Vector) -> Polynomial {
    assert_eq!(node_x.len(), node_y.len());
    let l = node_x.len();
    let mut result = Polynomial::new(vec![0f64; l]);

    for i in 0 .. l {
        let fixed_val = node_x[i];
        let prod = node_y[i];
        let mut id = poly(vec![1f64]);

        for j in 0 .. l {
            if j == i {
                continue;
            } else {
                let target_val = node_x[j];
                let denom = fixed_val - target_val;
                id = id * (poly(vec![1f64, - target_val]) / denom);
            }
        }
        result = result + (id * prod);
    }

    result
}
