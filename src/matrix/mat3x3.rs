use std::ops::*;
use std::convert::*;
#[allow(unused_imports)]
use std::{f32, f64};
use ::matrix::mat2x2::*;
use ::vector::vec2::*;
use ::vector::vec3::*;

#[derive(Copy, Clone)]
pub struct Mat3x3<T: Copy + Clone> {
    ///Array of Vector3<T> to facilitate the writing of operations
    pub m: [Vector3<T>; 3],
}

#[allow(dead_code)]
impl<T: Copy + Clone> Mat3x3<T> {

    pub fn new(m00: T, m01: T, m02: T, m10: T, m11: T, m12: T, m20: T, m21: T, m22: T) -> Mat3x3<T> {
        Mat3x3 {
            m: [
                Vector3 {x: m00, y: m01, z: m02},
                Vector3 {x: m10, y: m11, z: m12},
                Vector3 {x: m20, y: m21, z: m22}
            ]
        }
    }

    pub fn from_vectors(m0: Vector3<T>, m1: Vector3<T>, m2: Vector3<T>) -> Mat3x3<T> {
        Mat3x3 {
            m: [m0, m1, m2]
        }
    }

    ///Returns the dimensions of this matrix (for generic usage)
    pub fn dim(&self) -> (usize, usize) {
        (3, 3)
    }

    ///Returns a Mat3x3 which has been transposed
    pub fn as_transposed(&self) -> Mat3x3<T> {
        let mut m = Mat3x3 {
            m: self.m.clone()
        };

        m.transpose();
        m
    }

    pub fn transpose(&mut self) {
        let copy = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                self.m[i][j] = copy.m[j][i];
            }
        }
    }

    pub fn minor(&self, i: usize, j: usize) -> Mat2x2<T> {
        let col_0 = (i+1) % self.dim().0;
        let col_1 = (i+2) % self.dim().0;
        let row_0 = (j+1) % self.dim().1;
        let row_1 = (j+2) % self.dim().1;

        let (col_0, col_1) = if col_1 < col_0 { (col_1, col_0) } else { (col_0, col_1) };
        let (row_0, row_1) = if row_1 < row_0 { (row_1, row_0) } else { (row_0, row_1) };

        let m0: Vector2<T> = Vector2::new(self[col_0][row_0], self[col_0][row_1]);
        let m1: Vector2<T> = Vector2::new(self[col_1][row_0], self[col_1][row_1]);
        let sub_matrix: Mat2x2<T> = Mat2x2::from_vectors(m0, m1);

        sub_matrix
    }

}

impl<T: Copy + Clone> Index<usize> for Mat3x3<T> {
    type Output = Vector3<T>;

    fn index(&self, idx: usize) -> &Vector3<T> {
        &self.m[idx]
    }
}

impl<T: Copy + Clone> Index<(usize, usize)> for Mat3x3<T> {
    type Output = T;

    fn index(&self, idx: (usize, usize)) -> &T {
        return &self.m[idx.0][idx.1];
    }
}

impl<T: Copy + Clone> IndexMut<usize> for Mat3x3<T> {
    fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut Vector3<T> {
        &mut self.m[idx]
    }
}

impl<T: Copy + Clone> IndexMut<(usize, usize)> for Mat3x3<T> {
    fn index_mut<'a>(&'a mut self, idx: (usize, usize)) -> &'a mut T {
        &mut self.m[idx.0][idx.1]
    }
}

impl<T: Copy + Clone + Add<Output = T>> Add for Mat3x3<T> {
    type Output = Mat3x3<T>;

    fn add(self, rhs: Self) -> Mat3x3<T> {
        Mat3x3 {
            m: [
                self[0] + rhs[0],
                self[1] + rhs[1],
                self[2] + rhs[2]
            ]
        }
    }
}

impl<T: Copy + Clone + Add<Output = T>> AddAssign for Mat3x3<T> {
    fn add_assign(&mut self, rhs: Self) {
        self[0] = self[0] + rhs[0];
        self[1] = self[1] + rhs[1];
        self[2] = self[2] + rhs[2];
    }
}

impl<T: Copy + Clone + Sub<Output = T>> Sub for Mat3x3<T> {
    type Output = Mat3x3<T>;

    fn sub(self, rhs: Self) -> Mat3x3<T> {
        Mat3x3 {
            m: [
                self[0] - rhs[0],
                self[1] - rhs[1],
                self[2] - rhs[2]
            ]
        }
    }
}

impl<T: Copy + Clone + Sub<Output = T>> SubAssign for Mat3x3<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self[0] = self[0] - rhs[0];
        self[1] = self[1] - rhs[1];
        self[2] = self[2] - rhs[2];
    }
}

impl<T: Copy + Clone + Mul<Output = T>> Mul<T> for Mat3x3<T> {
    type Output = Mat3x3<T>;

    fn mul(self, rhs: T) -> Mat3x3<T> {
        Mat3x3 {
            m: [
                self.m[0] * rhs,
                self.m[1] * rhs,
                self.m[2] * rhs
            ]
        }
    }
}

impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T>> Mul<Vector3<T>> for Mat3x3<T> {
    type Output = Vector3<T>;

    fn mul(self, rhs: Vector3<T>) -> Vector3<T> {
        let self_transposed = self.as_transposed();

        Vector3 {
            x: self_transposed[0].dot(rhs),
            y: self_transposed[1].dot(rhs),
            z: self_transposed[2].dot(rhs)
        }
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> Mul for Mat3x3<T> {
    type Output = Mat3x3<T>;

    fn mul(self, rhs: Self) -> Mat3x3<T> {
        //Get a transposed version of self to do dot(row, column) for multiplication
        let self_transposed = self.as_transposed();

        let mut result: Mat3x3<T> = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                result[(i, j)] = self_transposed.m[i].dot(rhs.m[j]);
            }
        }

        result
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> MulAssign<T> for Mat3x3<T> {
    fn mul_assign(&mut self, rhs: T) {
        self[0] = self[0] * rhs;
        self[1] = self[1] * rhs;
        self[2] = self[2] * rhs;
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> MulAssign for Mat3x3<T> {
    fn mul_assign(&mut self, rhs: Self) {
        let self_transposed = self.as_transposed();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                self[(i, j)] = self_transposed.m[i].dot(rhs.m[j]);
            }
        }
    }
}

impl<T: Copy + Clone + Div<Output = T>> Div<T> for Mat3x3<T> {
    type Output = Mat3x3<T>;

    fn div(self, rhs: T) -> Mat3x3<T> {
        Mat3x3 {
            m: [
                self[0] / rhs,
                self[1] / rhs,
                self[2] / rhs
            ]
        }
    }
}

impl<T: Copy + Clone + Div<Output = T>> DivAssign<T> for Mat3x3<T> {
    fn div_assign(&mut self, rhs: T) {
        self[0] = self[0] / rhs;
        self[1] = self[1] / rhs;
        self[2] = self[2] / rhs;
    }
}

impl<T: Copy + Clone + Neg<Output = T>> Neg for Mat3x3<T> {
    type Output = Mat3x3<T>;

    fn neg(self) -> Mat3x3<T> {
        Mat3x3 {
            m: [
                -self[0],
                -self[1],
                -self[2]
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Neg<Output = T> + Mul<Output = T> + Sub<Output = T>> Mat3x3<T> {
    pub fn adjugate(&self) -> Mat3x3<T> {
        let mut adj = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                let negative: bool = (i + j) % 2 == 0;
                let sub_matrix = self.minor(i, j);

                if negative {
                    adj[i][j] = -sub_matrix.det();
                } else {
                    adj[i][j] = sub_matrix.det();
                }
            }
        }

        adj.transpose();

        adj
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Sub<Output = T>> Mat3x3<T> {
    pub fn from_diagonal(diag: Vector3<T>) -> Mat3x3<T> {
        let zero = diag.x - diag.x;

        Mat3x3 {
            m: [
                Vector3 {x: diag.x, y: zero, z: zero},
                Vector3 {x: zero, y: diag.y, z: zero},
                Vector3 {x: zero, y: zero, z: diag.z}
            ]
        }
    }

    pub fn scale(s: T) -> Mat3x3<T> {
        let zero: T = s - s;

        Mat3x3 {
            m: [
                Vector3 {x: s, y: zero, z: zero},
                Vector3 {x: zero, y: s, z: zero},
                Vector3 {x: zero, y: zero, z: s}
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + From<u8>> Mat3x3<T> {
    pub fn identity() -> Mat3x3<T> {
        let one: T = From::from(1u8);
        let zero: T = From::from(0u8);

        Mat3x3 {
            m: [
                Vector3 {x: one, y: zero, z: zero},
                Vector3 {x: zero, y: one, z: zero},
                Vector3 {x: zero, y: zero, z: one}
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Neg<Output = T>> Mat3x3<T> {
    pub fn det(&self) -> T {
        let mut _det: T = self[0][0] - self[0][0];

        for i in 0..self.dim().0 {
            let m = self.minor(i, 0);

            _det = if (i % 2) == 1 {_det + (m.det() * self[i][0])} else {_det - (m.det() * self[i][0])};
        }

        return _det
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Neg<Output = T> + PartialEq> Mat3x3<T> {
    pub fn inverse(&self) -> Mat3x3<T> {
        let zero: T = self[0][0] - self[0][0];
        let determinate = self.det();

        if determinate == zero {
            return Mat3x3 {
                m: [
                    Vector3 {x: zero, y: zero, z: zero},
                    Vector3 {x: zero, y: zero, z: zero},
                    Vector3 {x: zero, y: zero, z: zero}
                ]
            }
        } else {
            return self.adjugate() / determinate;
        }
    }
}

#[allow(dead_code)]
impl Mat3x3<f32> {

    pub fn pitch(angle_rad: f32) -> Mat3x3<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat3x3 {
            m: [
                Vector3 {x: 1.0f32, y: 0.0f32, z: 0.0f32},
                Vector3 {x: 0.0f32, y: cos,    z: sin   },
                Vector3 {x: 0.0f32, y: -sin,   z: cos   }
            ]
        }
    }

    pub fn yaw(angle_rad: f32) -> Mat3x3<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat3x3 {
            m: [
                Vector3 {x: cos,    y: 0.0f32, z: -sin  },
                Vector3 {x: 0.0f32, y: 1.0f32, z: 0.0f32},
                Vector3 {x: sin,    y: 0.0f32, z: cos   }
            ]
        }
    }

    pub fn roll(angle_rad: f32) -> Mat3x3<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat3x3 {
            m: [
                Vector3 {x: cos,    y: sin,    z: 0.0f32},
                Vector3 {x: -sin,   y: cos,    z: 0.0f32},
                Vector3 {x: 0.0f32, y: 0.0f32, z: 1.0f32}
            ]
        }
    }

    pub fn angle_axis(angle_rad: f32, axis: Vector3<f32>) -> Mat3x3<f32> {
        let axis = axis.norm();

        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        let x_sqr = axis.x * axis.x;
        let y_sqr = axis.y * axis.y;
        let z_sqr = axis.z * axis.z;

        let versine = 1.0f32 - cos;

        let x_basis: Vector3<f32> = Vector3 {
            x: cos + (x_sqr * versine),
            y: (axis.y * axis.x * versine) + (axis.z * sin),
            z: (axis.z * axis.x * versine) - (axis.y * sin)
        };

        let y_basis: Vector3<f32> = Vector3 {
            x: (axis.x * axis.y * versine) - (axis.z * sin),
            y: cos + (y_sqr * versine),
            z: (axis.z * axis.y * versine) + (axis.x * sin)
        };

        let z_basis: Vector3<f32> = Vector3 {
            x: (axis.x * axis.z * versine) + (axis.y * sin),
            y: (axis.y * axis.z * versine) - (axis.x * sin),
            z: cos + (z_sqr * versine)
        };

        Mat3x3 {
            m: [
                x_basis,
                y_basis,
                z_basis
            ]
        }
    }

}

#[allow(dead_code)]
impl Mat3x3<f64> {
    pub fn pitch(angle_rad: f64) -> Mat3x3<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat3x3 {
            m: [
                Vector3 {x: 1.0f64, y: 0.0f64, z: 0.0f64},
                Vector3 {x: 0.0f64, y: cos,    z: sin   },
                Vector3 {x: 0.0f64, y: -sin,   z: cos   }
            ]
        }
    }

    pub fn yaw(angle_rad: f64) -> Mat3x3<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat3x3 {
            m: [
                Vector3 {x: cos,    y: 0.0f64, z: -sin  },
                Vector3 {x: 0.0f64, y: 1.0f64, z: 0.0f64},
                Vector3 {x: sin,    y: 0.0f64, z: cos   }
            ]
        }
    }

    pub fn roll(angle_rad: f64) -> Mat3x3<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat3x3 {
            m: [
                Vector3 {x: cos,    y: sin,    z: 0.0f64},
                Vector3 {x: -sin,   y: cos,    z: 0.0f64},
                Vector3 {x: 0.0f64, y: 0.0f64, z: 1.0f64}
            ]
        }
    }

    pub fn angle_axis(angle_rad: f64, axis: Vector3<f64>) -> Mat3x3<f64> {
        let axis = axis.norm();

        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        let x_sqr = axis.x * axis.x;
        let y_sqr = axis.y * axis.y;
        let z_sqr = axis.z * axis.z;

        let versine = 1.0f64 - cos;

        let x_basis: Vector3<f64> = Vector3 {
            x: cos + (x_sqr * versine),
            y: (axis.y * axis.x * versine) + (axis.z * sin),
            z: (axis.z * axis.x * versine) - (axis.y * sin)
        };

        let y_basis: Vector3<f64> = Vector3 {
            x: (axis.x * axis.y * versine) - (axis.z * sin),
            y: cos + (y_sqr * versine),
            z: (axis.z * axis.y * versine) + (axis.x * sin)
        };

        let z_basis: Vector3<f64> = Vector3 {
            x: (axis.x * axis.z * versine) + (axis.y * sin),
            y: (axis.y * axis.z * versine) - (axis.x * sin),
            z: cos + (z_sqr * versine)
        };

        Mat3x3 {
            m: [
                x_basis,
                y_basis,
                z_basis
            ]
        }
    }
}

type Mat3x3f = Mat3x3<f32>;
type Mat3x3d = Mat3x3<f64>;

type Mat3x3ll = Mat3x3<i64>;
type Mat3x3i = Mat3x3<i32>;
type Mat3x3s = Mat3x3<i16>;
type Mat3x3b = Mat3x3<i8>;

type Mat3x3ull = Mat3x3<u64>;
type Mat3x3ui = Mat3x3<u32>;
type Mat3x3us = Mat3x3<u16>;
type Mat3x3ub = Mat3x3<u8>;
