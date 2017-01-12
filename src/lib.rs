mod vector;
mod matrix;

#[cfg(test)]
mod tests {

    #[allow(unused_imports)]
    use vector::*;
    #[allow(unused_imports)]
    use matrix::*;

    #[test]
    fn test_index() {
        let v2: Vec2f = Vec2f::new(1.0f32, 2.0f32);
        let v3: Vec3f = Vec3f::new(1.0f32, 2.0f32, 3.0f32);
        let v4: Vec4f = Vec4f::new(1.0f32, 2.0f32, 3.0f32, 4.0f32);

        print!("v2: ({}", v2[0]);
        for i in 1..v2.dim() {
            print!(", {}", v2[i]);
        }
        println!(")");

        print!("v3: ({}", v3[0]);
        for i in 1..v3.dim() {
            print!(", {}", v3[i]);
        }
        println!(")");

        print!("v4: ({}", v4[0]);
        for i in 1..v4.dim() {
            print!(", {}", v4[i]);
        }
        println!(")");
    }

}
