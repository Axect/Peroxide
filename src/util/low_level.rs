/// Copy each value from `src` through the corresponding pointer in `dst`.
///
/// # Safety
///
/// * `dst` and `src` must have equal length (enforced by debug assert).
/// * Every pointer in `dst` must be valid for writes of `T` and properly
///   aligned. Aliased pointers within `dst` produce undefined behaviour.
/// * No concurrent reader may access the locations written through `dst`.
pub unsafe fn copy_vec_ptr<T>(dst: &mut [*mut T], src: &[T])
where
    T: Copy,
{
    assert_eq!(dst.len(), src.len(), "Should use same length vectors");
    for (&mut p, &s) in dst.iter_mut().zip(src) {
        *p = s;
    }
}

/// Pairwise-swap the values pointed to by `lhs` and `rhs`.
///
/// # Safety
///
/// * `lhs` and `rhs` must have equal length (enforced by debug assert).
/// * Each `(lhs[i], rhs[i])` pair must point to valid, properly-aligned
///   `T` locations and must not alias another live reference.
/// * `lhs` and `rhs` themselves must point to disjoint memory.
pub unsafe fn swap_vec_ptr<T>(lhs: &mut [*mut T], rhs: &mut [*mut T])
where
    T: Copy,
{
    assert_eq!(lhs.len(), rhs.len(), "Should use same length vectors");
    for (&mut l, &mut r) in lhs.iter_mut().zip(rhs.iter_mut()) {
        std::ptr::swap(l, r);
    }
}

/// Dereference every pointer in `pv` and collect the values into a `Vec`.
///
/// # Safety
///
/// Every pointer in `pv` must be valid for reads of `T`, properly
/// aligned, and must outlive the call. The pointed-to data must not be
/// mutated through another reference during the call.
pub unsafe fn ptr_to_vec<T>(pv: &[*const T]) -> Vec<T>
where
    T: Copy,
{
    pv.iter().map(|&x| *x).collect()
}
