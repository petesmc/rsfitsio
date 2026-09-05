//! Handing a `Vec`'s storage out as a raw pointer, and reclaiming it later.
//!
//! CFITSIO's structs and its API traffic in bare pointers: `FITSfile` stores a
//! `char *filename` and a `tcolumn *tableptr`, and routines like `ffgkls` and
//! `fits_hdr2str` hand the caller a heap block to release with `fffree`. In C
//! that costs nothing, because `free` needs only the address.
//!
//! Rust's allocator needs more: [`alloc::alloc::dealloc`] must be given the same
//! [`Layout`](core::alloc::Layout) the block was allocated with, and
//! `Vec::from_raw_parts` needs the capacity, neither of which can be recovered
//! from the pointer. This module is the side channel that carries them from the
//! hand-out site to the reclaim site.
//!
//! Reach for it only when nothing simpler fits. In order of preference:
//!
//! 1. The block never escapes: use a plain `Vec<T>` or `Box<[T]>` and let RAII
//!    free it.
//! 2. A raw pointer has to be visible, but something we own can hold the
//!    storage: keep the `Vec<T>` in a local or a struct field and take
//!    `as_mut_ptr()` at the point of use.
//! 3. The pointer really is the handle, and the length is recoverable at the
//!    free site (from a sibling field, or from a NUL terminator): use
//!    `Box::into_raw(v.into_boxed_slice())` and
//!    `Box::from_raw(ptr::slice_from_raw_parts_mut(p, len))`. A `Box<[T]>` has
//!    no spare capacity, so its layout is exactly `Layout::array::<T>(len)`.
//! 4. Otherwise, this module.
//!
//! There is deliberately no "assume the length" fallback. Every reclaim site
//! used to end in `else { Vec::from_raw_parts(p, <guess>, <guess>) }`, and a
//! guess that is wrong is a wrong `Layout` handed to the allocator. An
//! unregistered pointer is a bug, so it is reported in debug builds and leaked
//! rather than freed.
#![warn(missing_docs)]

use std::collections::HashMap;
use std::sync::{LazyLock, Mutex};

use crate::helpers::vec_raw_parts::vec_into_raw_parts;

/// What was recorded about one handed-out allocation.
///
/// `len` and `capacity` are in elements. The element type is not part of the
/// key, so `elem_size`/`elem_align` are kept to catch a pointer registered as
/// one type and reclaimed as another -- a mismatch the address alone cannot
/// reveal, and which would deallocate with the wrong `Layout`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Record {
    len: usize,
    capacity: usize,
    elem_size: usize,
    elem_align: usize,
}

impl Record {
    fn of<T>(len: usize, capacity: usize) -> Self {
        Record {
            len,
            capacity,
            elem_size: size_of::<T>(),
            elem_align: align_of::<T>(),
        }
    }

    fn matches<T>(&self) -> bool {
        self.elem_size == size_of::<T>() && self.elem_align == align_of::<T>()
    }
}

/// Address -> layout, for every allocation this crate has handed out as a raw
/// pointer.
static ALLOCATIONS: LazyLock<Mutex<HashMap<usize, Record>>> = LazyLock::new(Default::default);

/// Hand `v`'s storage out as a raw pointer, recording its layout.
///
/// The returned pointer must eventually reach [`take_registered`] or
/// [`free_registered`]; until then the allocation is owned by nothing and will
/// leak.
pub(crate) fn into_raw_registered<T>(v: Vec<T>) -> *mut T {
    let (p, len, capacity) = vec_into_raw_parts(v);
    ALLOCATIONS
        .lock()
        .unwrap()
        .insert(p as usize, Record::of::<T>(len, capacity));
    p
}

/// Reclaim a pointer handed out by [`into_raw_registered`], as the `Vec` it
/// came from.
///
/// Returns `None` -- leaking the block rather than freeing it with a layout
/// that may be wrong -- for a null pointer, one that was never registered, and
/// one whose element type does not match `T`. The last two are bugs and are
/// reported by a debug assertion.
///
/// # Safety
/// `p` must not be used again after this returns `Some`, and no other live
/// reference may point into the block.
pub(crate) unsafe fn take_registered<T>(p: *mut T) -> Option<Vec<T>> {
    if p.is_null() {
        return None;
    }

    let Some(record) = ALLOCATIONS.lock().unwrap().remove(&(p as usize)) else {
        debug_assert!(false, "{p:p} was never registered with into_raw_registered");
        return None;
    };

    if !record.matches::<T>() {
        debug_assert!(
            false,
            "{p:p} was registered as {}-byte elements but is being reclaimed as {}",
            record.elem_size,
            core::any::type_name::<T>()
        );
        return None;
    }

    // SAFETY: the record is the (len, capacity) the Vec was taken apart with,
    // and the element size and alignment have just been checked to match.
    Some(unsafe { Vec::from_raw_parts(p, record.len, record.capacity) })
}

/// Reclaim and drop a pointer handed out by [`into_raw_registered`].
///
/// The counterpart of C's `free`. A null pointer is ignored, as `free(NULL)` is.
///
/// # Safety
/// As [`take_registered`].
pub(crate) unsafe fn free_registered<T>(p: *mut T) {
    drop(unsafe { take_registered::<T>(p) });
}

/// Move the registration from `old` to `v`'s storage, returning the new
/// pointer.
///
/// This is the `realloc` shape: the caller has already reclaimed `old` (or is
/// replacing it wholesale) and wants the map to describe the replacement. `old`
/// is used only as a key, so passing a stale address is harmless.
pub(crate) fn reregister<T>(old: *mut T, v: Vec<T>) -> *mut T {
    let (p, len, capacity) = vec_into_raw_parts(v);
    let mut map = ALLOCATIONS.lock().unwrap();
    map.remove(&(old as usize));
    map.insert(p as usize, Record::of::<T>(len, capacity));
    p
}

/// The recorded element count of a handed-out allocation, for the callers that
/// need to view it as a slice without giving up ownership.
pub(crate) fn registered_len<T>(p: *const T) -> Option<usize> {
    ALLOCATIONS
        .lock()
        .unwrap()
        .get(&(p as usize))
        .map(|r| r.len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round_trip_preserves_contents_and_capacity() {
        let mut v: Vec<u32> = Vec::with_capacity(16);
        v.extend([1, 2, 3]);

        let p = into_raw_registered(v);
        assert_eq!(registered_len(p.cast_const()), Some(3));

        let back = unsafe { take_registered::<u32>(p) }.unwrap();
        assert_eq!(back, vec![1, 2, 3]);
        assert_eq!(back.capacity(), 16);

        /* the entry is gone, so a second reclaim finds nothing */
        assert_eq!(registered_len(p.cast_const()), None);
    }

    #[test]
    fn test_null_is_ignored() {
        assert!(unsafe { take_registered::<u8>(core::ptr::null_mut()) }.is_none());
        unsafe { free_registered::<u8>(core::ptr::null_mut()) };
    }

    #[test]
    fn test_reregister_replaces_the_entry() {
        let old = into_raw_registered(vec![0u8; 4]);
        /* the caller has taken the old block apart itself */
        let v = unsafe { Vec::from_raw_parts(old, 4, 4) };
        let mut v = v;
        v.resize(9, 7);

        let new = reregister(old, v);
        assert_eq!(registered_len(new.cast_const()), Some(9));

        let back = unsafe { take_registered::<u8>(new) }.unwrap();
        assert_eq!(back.len(), 9);
    }

    /// An unregistered pointer must not be deallocated on a guessed layout.
    /// The debug assertion fires first in a debug build, so this is the
    /// release-build behaviour being pinned.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "the debug assertion fires first")]
    fn test_unregistered_pointer_is_leaked_not_freed() {
        let mut v = vec![1u8, 2, 3];
        let p = v.as_mut_ptr();
        assert!(unsafe { take_registered::<u8>(p) }.is_none());
        /* v still owns its storage and drops normally */
        assert_eq!(v, vec![1, 2, 3]);
    }
}
