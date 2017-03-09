use std::ops::*;
use std::convert::*;
#[allow(unused_imports)]
use std::{f32, f64};
use ::vector::vec2::*;

#[derive(Copy, Clone)]
pub struct Mat2x2<T: Copy + Clone> {
    ///Array of Vector2<T> to facilitate the writing of operations
    pub m: [Vector2<T>; 2],
}

#[allow(dead_code)]
impl<T: Copy + Clone> Mat2x2<T> {

    pub fn new(m00: T, m01: T, m10: T, m11: T) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                Vector2 {x: m00, y: m01},
                Vector2 {x: m10, y: m11}
            ]
        }
    }

    pub fn from_vectors(m0: Vector2<T>, m1: Vector2<T>) -> Mat2x2<T> {
        Mat2x2 {
            m: [m0, m1]
        }
    }

    ///Returns the dimensions of this matrix (for generic usage)
    pub fn dim(&self) -> (usize, usize) {
        (2, 2)
    }

    ///Returns a Mat2x2 which has been transposed, i.e. m[0][1] and m[1][0] swapped for Mat2x2
    pub fn as_transposed(&self) -> Mat2x2<T> {
        let mut m = Mat2x2 {
            m: self.m.clone()
        };

        m.transpose();
        m
    }

    pub fn transpose(&mut self) {
        let tmp_01 = self.m[0][1];
        let tmp_10 = self.m[1][0];

        self.m[0][1] = tmp_10;
        self.m[1][0] = tmp_01;
    }

    pub fn submatrix(&self, i: usize, j: usize) -> T {
        let col = (i+1) % self.dim().0;
        let row = (j+1) % self.dim().1;

        self.m[col][row]
    }

}

impl<T: Copy + Clone> Index<usize> for Mat2x2<T> {
    type Output = Vector2<T>;

    fn index(&self, idx: usize) -> &Vector2<T> {
        &self.m[idx]
    }
}

impl<T: Copy + Clone> Index<(usize, usize)> for Mat2x2<T> {
    type Output = T;

    fn index(&self, idx: (usize, usize)) -> &T {
        return &self.m[idx.0][idx.1];
    }
}

impl<T: Copy + Clone> IndexMut<usize> for Mat2x2<T> {
    fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut Vector2<T> {
        &mut self.m[idx]
    }
}

impl<T: Copy + Clone> IndexMut<(usize, usize)> for Mat2x2<T> {
    fn index_mut<'a>(&'a mut self, idx: (usize, usize)) -> &'a mut T {
        &mut self.m[idx.0][idx.1]
    }
}

impl<T: Copy + Clone + Add<Output = T>> Add for Mat2x2<T> {
    type Output = Mat2x2<T>;

    fn add(self, rhs: Self) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                self[0] + rhs[0],
                self[1] + rhs[1]
            ]
        }
    }
}

impl<T: Copy + Clone + Add<Output = T>> AddAssign for Mat2x2<T> {
    fn add_assign(&mut self, rhs: Self) {
        self[0] = self[0] + rhs[0];
        self[1] = self[1] + rhs[1];
    }
}

impl<T: Copy + Clone + Sub<Output = T>> Sub for Mat2x2<T> {
    type Output = Mat2x2<T>;

    fn sub(self, rhs: Self) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                self[0] - rhs[0],
                self[1] - rhs[1]
            ]
        }
    }
}

impl<T: Copy + Clone + Sub<Output = T>> SubAssign for Mat2x2<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self[0] = self[0] - rhs[0];
        self[1] = self[1] - rhs[1];
    }
}

impl<T: Copy + Clone + Mul<Output = T>> Mul<T> for Mat2x2<T> {
    type Output = Mat2x2<T>;

    fn mul(self, rhs: T) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                self.m[0] * rhs,
                self.m[1] * rhs
            ]
        }
    }
}

impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T>> Mul<Vector2<T>> for Mat2x2<T> {
    type Output = Vector2<T>;

    fn mul(self, rhs: Vector2<T>) -> Vector2<T> {
        let self_transposed = self.as_transposed();

        Vector2 {
            x: self_transposed[0].dot(rhs),
            y: self_transposed[1].dot(rhs)
        }
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> Mul for Mat2x2<T> {
    type Output = Mat2x2<T>;

    fn mul(self, rhs: Self) -> Mat2x2<T> {
        //Get a transposed version of self to do dot(row, column) for multiplication
        let self_transposed = self.as_transposed();

        let mut result: Mat2x2<T> = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                result[(i, j)] = self_transposed.m[i].dot(rhs.m[j]);
            }
        }

        result
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> MulAssign<T> for Mat2x2<T> {
    fn mul_assign(&mut self, rhs: T) {
        self[0] = self[0] * rhs;
        self[1] = self[1] * rhs;
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> MulAssign for Mat2x2<T> {
    fn mul_assign(&mut self, rhs: Self) {
        let self_transposed = self.as_transposed();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                self[(i, j)] = self_transposed.m[i].dot(rhs.m[j]);
            }
        }
    }
}

impl<T: Copy + Clone + Div<Output = T>> Div<T> for Mat2x2<T> {
    type Output = Mat2x2<T>;

    fn div(self, rhs: T) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                self[0] / rhs,
                self[1] / rhs
            ]
        }
    }
}

impl<T: Copy + Clone + Div<Output = T>> DivAssign<T> for Mat2x2<T> {
    fn div_assign(&mut self, rhs: T) {
        self[0] = self[0] / rhs;
        self[1] = self[1] / rhs;
    }
}

impl<T: Copy + Clone + Neg<Output = T>> Neg for Mat2x2<T> {
    type Output = Mat2x2<T>;

    fn neg(self) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                -self[0],
                -self[1]
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Neg<Output = T>> Mat2x2<T> {
    pub fn adjugate(&self) -> Mat2x2<T> {
        Mat2x2 {
            m: [
                Vector2 {x: self[1][1], y: -self[0][1]},
                Vector2 {x: -self[1][0], y: self[0][0]}
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Sub<Output = T>> Mat2x2<T> {
    pub fn from_diagonal(diag: Vector2<T>) -> Mat2x2<T> {
        let zero = diag.x - diag.x;

        Mat2x2 {
            m: [
                Vector2 {x: diag.x, y: zero},
                Vector2 {x: zero, y: diag.y}
            ]
        }
    }

    pub fn scale(s: T) -> Mat2x2<T> {
        let zero: T = s - s;

        Mat2x2 {
            m: [
                Vector2 {x: s, y: zero},
                Vector2 {x: zero, y: s}
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + From<u8>> Mat2x2<T> {
    pub fn identity() -> Mat2x2<T> {
        let one: T = From::from(1u8);
        let zero: T = From::from(0u8);

        Mat2x2 {
            m: [
                Vector2 {x: one, y: zero},
                Vector2 {x: zero, y: one}
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Sub<Output = T>> Mat2x2<T> {
    pub fn det(&self) -> T {
        (self[0][0] * self[1][1]) - (self[1][0] * self[0][1])
    }

    pub fn minor(&self) -> Mat2x2<T> {
        let mut minor = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                let sub = self.submatrix(i, j);

                minor[i][j] = sub;
            }
        }

        minor
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Sub<Output = T> + Div<Output = T> + Neg<Output = T> + PartialEq> Mat2x2<T> {
    pub fn inverse(&self) -> Mat2x2<T> {
        let zero: T = self[0][0] - self[0][0];
        let determinate = self.det();

        if determinate == zero {
            return Mat2x2 {
                m: [
                    Vector2 {x: zero, y: zero},
                    Vector2 {x: zero, y: zero}
                ]
            }
        } else {
            return self.adjugate() / determinate;
        }
    }
}

#[allow(dead_code)]
impl Mat2x2<f32> {
    pub fn rotation(angle_rad: f32) -> Mat2x2<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat2x2 {
            m: [
                Vector2 {x: cos, y: sin},
                Vector2 {x: -sin, y: cos}
            ]
        }
    }


}

#[allow(dead_code)]
impl Mat2x2<f64> {
    pub fn rotation(angle_rad: f64) -> Mat2x2<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat2x2 {
            m: [
                Vector2 {x: cos, y: sin},
                Vector2 {x: -sin, y: cos}
            ]
        }
    }
}

type Mat2x2f = Mat2x2<f32>;
type Mat2x2d = Mat2x2<f64>;

type Mat2x2ll = Mat2x2<i64>;
type Mat2x2i = Mat2x2<i32>;
type Mat2x2s = Mat2x2<i16>;
type Mat2x2b = Mat2x2<i8>;

type Mat2x2ull = Mat2x2<u64>;
type Mat2x2ui = Mat2x2<u32>;
type Mat2x2us = Mat2x2<u16>;
type Mat2x2ub = Mat2x2<u8>;
