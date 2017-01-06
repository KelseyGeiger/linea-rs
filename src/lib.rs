mod vector;

#[cfg(test)]
mod tests {

    use vector::vector::Vector;
    
    #[test]
    fn test_vector_new() {
        let mut v: Vector<f32> = Vector::new(3);
        v.v[0] = 1f32;
        v.v[1] = 2f32;
        v.v[2] = 3f32;

        print!("v: ({}, {}, {})", v.v[0], v.v[1], v.v[2]);
    }

    #[test]
    fn test_vector_from_vec() {
        let mut v : Vector<u64> = Vector::from_vec(vec![1u64, 2u64]);

        print!("v: ({}, {})", v.v[0], v.v[1]);
    }

    #[test]
    fn test_vector_from_slice() {
        let slice = &[1.0f64, 0.0f64, 0.0f64, 1.0f64];
        let v = Vector::from_slice(slice);

        print!("v: ({}, {}, {}, {})", v.v[0], v.v[1], v.v[2], v.v[3]);
    }
}
