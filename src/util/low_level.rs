pub unsafe fn copy_vec_ptr(dst: &mut Vec<*mut f64>, src: &Vec<f64>) {
    assert_eq!(dst.len(), src.len(), "Should use same length vectors");
    for (&mut p, &s) in dst.iter_mut().zip(src) {
        *p = s;
    }
}

pub unsafe fn swap_vec_ptr(lhs: &mut Vec<*mut f64>, rhs: &mut Vec<*mut f64>) {
    assert_eq!(lhs.len(), rhs.len(), "Should use same length vectors");
    for (&mut l, &mut r) in lhs.iter_mut().zip(rhs.iter_mut()) {
        std::ptr::swap(l, r);
    }
}

pub unsafe fn ptr_to_vec(pv: &Vec<*const f64>) -> Vec<f64> {
    pv.iter().map(|&x| *x).collect()
}
