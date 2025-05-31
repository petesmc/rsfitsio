use std::alloc::Layout;
use std::alloc::alloc;
use std::error::Error;
use std::mem;
use std::ptr::NonNull;

/// The `AllocError` error indicates an allocation failure
/// that may be due to resource exhaustion or to
/// something wrong when combining the given input arguments with this
/// allocator.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct AllocError;

impl Error for AllocError {}

// (we need this for downstream impl of trait Error)
impl core::fmt::Display for AllocError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("memory allocation failed")
    }
}

pub(crate) fn box_try_new<T>(x: T) -> Result<Box<T>, AllocError> {
    let mut boxed = try_new_uninit()?;
    boxed.write(x);
    unsafe { Ok(boxed.assume_init()) }
}

fn try_new_uninit<T>() -> Result<Box<mem::MaybeUninit<T>>, AllocError> {
    let ptr = if size_of::<T>() == 0 {
        unsafe { mem::transmute::<NonNull<T>, *mut mem::MaybeUninit<T>>(NonNull::dangling()) }
    } else {
        let layout = Layout::new::<mem::MaybeUninit<T>>();
        let ptr = unsafe { alloc(layout) };
        if ptr.is_null() {
            return Err(AllocError);
        }
        (ptr) as *mut mem::MaybeUninit<T>
    };

    unsafe { Ok(Box::from_raw(ptr)) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_new() {
        let boxed: Box<i32> = box_try_new(42).expect("Allocation failed");
        assert_eq!(*boxed, 42);
    }
}
