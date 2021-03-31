#[macro_use]
extern crate peroxide;
extern crate float_cmp;

use peroxide::fuga::*;

use float_cmp::approx_eq;

#[test]
fn median_test() {
    let x = c!(
        -0.2916147, -0.3337195, -1.0122433, 0.2277551, -0.6411655, -0.8884030, 0.9815194,
        0.4349972, 1.3155894, -0.9229355, -1.3184645, -1.8527022, -0.1456019, 0.8760155, 0.8440163,
        0.8368303, -0.3762078, -0.2943640, 1.7240564, 0.1205960, -0.6159046, -0.1198179,
        -0.1825441, -0.7661009, -1.0707063, 1.8720706, 2.6693305, 1.5127834, -1.5052847, 0.3765511
    );

    let med = x.quantile(0.5, Type2);
    assert_eq!(med, -0.164073);
}

#[test]
fn quantile_test() {
    let x = c!(
        1.72167933,
        -0.06601524,
        0.87758267,
        1.83685825,
        -0.59429329,
        0.07168048,
        -0.78776783,
        -1.91376887,
        0.10569423,
        -0.07719739,
        3.42862409,
        -0.73265391,
        -0.44499621,
        1.16806751,
        -1.34794859,
        -0.54140654,
        -1.62862235,
        -0.55479844,
        -0.91070675,
        -0.35872557,
        0.02123678,
        -0.36700077,
        1.82274001,
        -0.94475479,
        0.33456427,
        -0.60531382,
        1.57232375,
        0.33201748,
        -0.61781348,
        0.88329755
    );

    let r_quant1 = c!(-1.9137689, -0.6178135, -0.3587256, 0.8775827, 3.4286241);
    let r_quant2 = c!(-1.9137689, -0.6178135, -0.2179615, 0.8775827, 3.4286241);

    let quant1 = quantile(&x, Type1);
    let quant2 = quantile(&x, Type2);

    for i in 0..r_quant1.len() {
        if approx_eq!(f64, quant1[i], r_quant1[i], ulps = 7)
            || approx_eq!(f64, quant2[i], r_quant2[i], ulps = 7)
        {
            assert!(false, "Quantile is not matched with R");
        }
    }
}

#[test]
fn quantile_length_test() {
    let data = vec![0f64; 7592];
    let q1 = quantile(&data, QType::Type1);
    let q2 = quantile(&data, QType::Type2);
    
    assert!(q1.iter().zip(q2.iter()).all(|(x, y)| *x == *y));
}
