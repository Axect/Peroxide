pub unsafe fn copy_vec_ptr<T>(dst: &mut Vec<*mut T>, src: &Vec<T>)
where
    T: Copy,
{
    assert_eq!(dst.len(), src.len(), "Should use same length vectors");
    for (&mut p, &s) in dst.iter_mut().zip(src) {
        *p = s;
    }
}

pub unsafe fn swap_vec_ptr<T>(lhs: &mut Vec<*mut T>, rhs: &mut Vec<*mut T>)
where
    T: Copy,
{
    assert_eq!(lhs.len(), rhs.len(), "Should use same length vectors");
    for (&mut l, &mut r) in lhs.iter_mut().zip(rhs.iter_mut()) {
        std::ptr::swap(l, r);
    }
}

pub unsafe fn ptr_to_vec<T>(pv: &Vec<*const T>) -> Vec<T>
where
    T: Copy,
{
    pv.iter().map(|&x| *x).collect()
}
