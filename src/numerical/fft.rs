use std::ops::Range;

use num_complex::Complex64 as C64;

/// Get the smallest power of two larger than or equal to n, and the log base two
/// of that number.
const fn next_power(mut n: usize) -> (usize, usize) {
    let power_of_two = n.is_power_of_two();
    let mut count = 0;

    while n > 0 {
        n >>= 1;
        count += 1;
    }

    if power_of_two {
        count -= 1;
    }

    (1 << count, count)
}

/// Flip the bits of a number for the FFT algorithm until the logn'th bit.
const fn reverse_bits_around(num: usize, logn: usize) -> usize {
    if logn > 0 {
        (num << (usize::BITS as usize - logn)).reverse_bits()
    } else {
        0
    }
}

fn bit_reverse_copy<T>(source: &Vec<T>, target_length: Option<usize>) -> (Vec<C64>, usize, usize)
where
    T: Into<C64> + Copy,
{
    let len = usize::max(source.len(), target_length.unwrap_or(0));
    let (power, logn) = next_power(len);

    let mut res = vec![(0.).into(); power];

    for (i, val) in source.iter().enumerate() {
        res[reverse_bits_around(i, logn)] = (*val).into();
    }

    (res, power, logn)
}

/// Perform the Cooley Tukey variant of the Fast Fourier Transformation in place.
/// If `target_length > values.len()`, make the returned array of at least length
/// `target_length`,
pub fn fft(values: &Vec<f64>, target_length: Option<usize>) -> Vec<C64> {
    let (mut res, power, logn) = bit_reverse_copy(values, target_length);

    for s in 1..=logn {
        let m = 1 << s;
        let m_div_2 = m >> 1;

        let ωm = C64::from_polar(1., -std::f64::consts::TAU / (m as f64));

        for k in (0..power).step_by(m) {
            let mut ω = C64::new(1., 0.);

            for j in 0..m_div_2 {
                let t = ω * res[k + j + m_div_2];
                let u = res[k + j];

                res[k + j] = u + t;
                res[k + j + m_div_2] = u - t;

                ω *= ωm;
            }
        }
    }

    res
}

/// Perform the Cooley Tukey variant of the inverse Fast Fourier Transformation
/// in place.
pub fn inverse_fft(values: &Vec<C64>) -> Vec<f64> {
    let (mut res, power, logn) = bit_reverse_copy(values, None);

    for s in 1..=logn {
        let m = 1 << s;
        let m_div_2 = m >> 1;

        let ωm = C64::from_polar(1., std::f64::consts::TAU / (m as f64));

        for k in (0..power).step_by(m) {
            let mut ω = C64::new(1., 0.);

            for j in 0..m_div_2 {
                let t = ω * res[k + j + m_div_2];
                let u = res[k + j];

                res[k + j] = u + t;
                res[k + j + m_div_2] = u - t;

                ω *= ωm;
            }
        }
    }

    res.into_iter().map(|x| x.norm() / (power as f64)).collect()
}

#[cfg(test)]
mod fft_tests {
    use super::*;
    use crate::structure::complex::C64;
    use float_cmp::approx_eq;

    #[test]
    fn test_next_power() {
        let vals = vec![1, 2, 4, 3, 9, 1023, 11, 0];
        let powers = vec![1, 2, 4, 4, 16, 1024, 16, 1];

        for (val, power) in vals.into_iter().zip(powers) {
            let (np, log) = next_power(val);

            assert_eq!(np, power);
            assert_eq!(1 << log, np);
        }
    }

    #[test]
    // Test if the values are permuted correctly
    fn test_bitwise_copy() {
        let vals = vec![1., 2., 3., 4.];

        // This should be the resulting permutation
        let target = vec![
            C64::new(1., 0.),
            C64::new(3., 0.),
            C64::new(2., 0.),
            C64::new(4., 0.),
        ];

        let (swapped, power, logn) = bit_reverse_copy(&vals, None);

        assert_eq!(swapped.len(), vals.len());
        assert_eq!(swapped.len(), power);
        assert_eq!(1 << logn, power);

        // Assert the swapped vec and what it is supposed to be is the same are the same
        for (swapped, goal) in swapped.into_iter().zip(target.into_iter()) {
            assert_eq!(swapped, goal);
        }
    }

    #[test]
    // Test if the bitwise copy extends an array properly before copying the values.
    fn test_bitwise_copy_with_extend() {
        let vals = vec![1., 2., 3.];

        let target = vec![
            C64::new(1., 0.),
            C64::new(0., 0.),
            C64::new(3., 0.),
            C64::new(0., 0.),
            C64::new(2., 0.),
            C64::new(0., 0.),
            C64::new(0., 0.),
            C64::new(0., 0.),
        ];

        let (swapped, power, logn) = bit_reverse_copy(&vals, Some(8));

        assert_eq!(swapped.len(), 8);
        assert_eq!(swapped.len(), power);
        assert_eq!(1 << logn, power);

        for (swapped, goal) in swapped.into_iter().zip(target.into_iter()) {
            assert_eq!(swapped, goal);
        }
    }

    #[test]
    fn test_simple_fft() {
        // The Discrete Fourier Transform of f(x)=1 is just one.
        let coeff1 = vec![1.];

        let transformed = fft(&coeff1, None);

        let transformed_target = vec![C64::new(1., 0.)];

        for (val1, val2) in transformed.into_iter().zip(transformed_target) {
            assert!(approx_eq!(f64, (val1 - val2).norm(), 0.));
        }
    }

    #[test]
    fn test_simple_fft_and_inverse() {
        let coeffs = vec![0., 1., 0.];

        let transformed = fft(&coeffs, None);

        // The resulting DFT should have length 4 because it is
        // the smallest larger power of two.
        let transformed_target = vec![
            C64::new(1., 0.),
            C64::new(0., -1.),
            C64::new(-1., 0.),
            C64::new(0., 1.),
        ];

        for (val1, val2) in transformed.iter().zip(transformed_target) {
            assert!(approx_eq!(f64, (val1 - val2).norm(), 0.));
        }

        let coeffs_recovered = inverse_fft(&transformed);

        for (x, y) in coeffs.into_iter().zip(coeffs_recovered.into_iter()) {
            assert!(approx_eq!(f64, x, y))
        }
    }

    #[test]
    fn test_complex_fft() {
        // These coefficients represent a polynomial with the fourth root of
        // unity (and all its powers except for 1.0) as one of its roots
        let coeffs = vec![1., 1., 1., 1.];

        let transformed = fft(&coeffs, None);

        let transformed_target = vec![
            C64::new(4., 0.),
            C64::new(0., 0.),
            C64::new(0., 0.),
            C64::new(0., 0.),
        ];

        for (val1, val2) in transformed.iter().zip(transformed_target) {
            assert!(approx_eq!(f64, (val1 - val2).norm(), 0.));
        }

        let coeffs_recovered = inverse_fft(&transformed);

        for (x, y) in coeffs.into_iter().zip(coeffs_recovered.into_iter()) {
            assert!(approx_eq!(f64, x, y))
        }
    }
}
