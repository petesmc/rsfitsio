use std::mem::ManuallyDrop;

/// Decomposes a `Vec<T>` into its raw components: `(pointer, length, capacity)`.
///
/// Returns the raw pointer to the underlying data, the length of
/// the vector (in elements), and the allocated capacity of the
/// data (in elements). These are the same arguments in the same
/// order as the arguments to [`from_raw_parts`].
///
/// After calling this function, the caller is responsible for the
/// memory previously managed by the `Vec`. The only way to do
/// this is to convert the raw pointer, length, and capacity back
/// into a `Vec` with the [`from_raw_parts`] function, allowing
/// the destructor to perform the cleanup.
#[must_use = "losing the pointer will leak memory"]
pub(crate) fn vec_into_raw_parts<T>(v: Vec<T>) -> (*mut T, usize, usize) {
    let mut me = ManuallyDrop::new(v);
    (me.as_mut_ptr(), me.len(), me.capacity())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec_into_raw_parts() {
        let v = vec![1, 2, 3];
        let (ptr, len, cap) = vec_into_raw_parts(v);

        assert_eq!(len, 3);
        assert_eq!(cap, 3);
        assert!(!ptr.is_null());

        // Convert back to Vec to ensure it works
        let v_back = unsafe { Vec::from_raw_parts(ptr, len, cap) };
        assert_eq!(v_back, vec![1, 2, 3]);
    }
}
