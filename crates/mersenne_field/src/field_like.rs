use crate::field::Field;

pub trait FieldLikeVectorized: Field {
    type Base: Field;
    const SIZE_FACTOR: usize = core::mem::size_of::<Self>() / core::mem::size_of::<Self::Base>();

    fn constant(value: Self::Base) -> Self;

    fn get_base_element(&self, idx: usize) -> Self::Base;

    fn from_base_elements(input: &[Self::Base]) -> Self;

    fn from_base_array(input: &[Self::Base; Self::SIZE_FACTOR]) -> Self;

    // fn from_base_elements(input: &[Self::Base]) -> Self {
    //     unsafe {
    //         core::slice::from_raw_parts_mut(input.as_ptr() as *mut Self, 1)[0]
    //     }
    // }

    // fn from_base_array(input: &[Self::Base; Self::SIZE_FACTOR]) -> Self {
    //     unsafe {
    //         core::slice::from_raw_parts_mut(input.as_ptr() as *mut Self, 1)[0]
    //     }
    // }

    fn as_base_array(&self) -> [Self::Base; Self::SIZE_FACTOR] {
        let as_slice = core::slice::from_ref(self);
        let as_base_slice =
            unsafe { core::slice::from_raw_parts(as_slice.as_ptr() as *mut Self::Base, 16) };
        as_base_slice.try_into().unwrap()
    }
}
