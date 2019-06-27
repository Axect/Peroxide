use structure::matrix::*;
use structure::vector::*;
use util::non_macro::*;

pub fn rk4<F>(init_param: f64, init_value: Vec<f64>, f: F, step: f64, num: usize) -> Matrix
where
    F: Fn(f64, Vec<f64>) -> Vec<f64> + Copy,
{
    let mut t = init_param;
    let mut xs = init_value.clone();
    let mut records = zeros(num+1, init_value.len() + 1);
    records.subs_row(0, concat(vec![t], xs.clone()));
    for i in 1 .. num+1 {
        xs = one_step_rk(t, xs.clone(), f, step);
        t += step;
        records.subs_row(i, concat(vec![t], xs.clone()));
    }
    records
}

pub fn one_step_rk<F>(t: f64, xs: Vec<f64>, f: F, step: f64) -> Vec<f64>
where
    F: Fn(f64, Vec<f64>) -> Vec<f64> + Copy,
{
    let ys = xs.clone();
    let h = step;
    let h2 = h / 2f64;
    let h3 = h / 3f64;
    let h6 = h / 6f64;

    let k1 = f(t, xs.clone());

    let t1 = t + h2;
    let t2 = t + h;

    let k2_add = k1.fmap(|x| x * h);
    let k2 = f(t1, ys.add(&k2_add));

    let k3_add = k2.fmap(|x| x * h2);
    let k3 = f(t1, ys.add(&k3_add));

    let k4_add = k3.fmap(|x| x * h2);
    let k4 = f(t2, ys.add(&k4_add));

    let total_add_part1 = k1.zip_with(|x, y| h6 * x + h3 * y, &k2);
    let total_add_part2 = k3.zip_with(|x, y| h3 * x + h6 * y, &k4);

    let total_add = total_add_part1.zip_with(|x, y| x + y, &total_add_part2);

    let new_ys = ys.add(&total_add);

    new_ys
}