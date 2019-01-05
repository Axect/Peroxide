use structure::matrix::*;
use structure::dual::*;
use structure::vector::*;
use util::non_macro::*;

pub fn rk4<F>(init: Vec<f64>, f: F, step: f64, num: usize) -> Matrix
where
    F: Fn(Vec<f64>) -> Vec<f64> + Copy,
{
    let mut xs = init.clone();
    let mut records = zeros(num+1, init.len());
    records.subs_row(0, xs.clone());
    for i in 1 .. num+1 {
        xs = one_step_rk(xs.clone(), f, step);
        records.subs_row(i, xs.clone());
    }
    records
}

pub fn one_step_rk<F>(xs: Vec<f64>, f: F, step: f64) -> Vec<f64>
where
    F: Fn(Vec<f64>) -> Vec<f64> + Copy,
{
    let t = xs[0];
    let ys = xs.clone().into_iter().skip(1).collect::<Vec<f64>>();
    let h = step;
    let h2 = h / 2f64;
    let h3 = h / 3f64;
    let h6 = h / 6f64;

    let k1 = f(xs.clone());

    let t1 = t + h2;
    let t2 = t + h;

    let k2_add = k1.fmap(|x| x * h);
    let k2 = f(concat(vec![t1], ys.add(&k2_add)));

    let k3_add = k2.fmap(|x| x * h2);
    let k3 = f(concat(vec![t1], ys.add(&k3_add)));

    let k4_add = k3.fmap(|x| x * h2);
    let k4 = f(concat(vec![t2], ys.add(&k4_add)));

    let total_add_part1 = k1.zip_with(|x, y| h6 * x + h3 * y, &k2);
    let total_add_part2 = k3.zip_with(|x, y| h3 * x + h6 * y, &k4);

    let total_add = total_add_part1.zip_with(|x, y| x + y, &total_add_part2);

    let new_t = t + h;
    let new_ys = ys.add(&total_add);

    concat(vec![new_t], new_ys)
}