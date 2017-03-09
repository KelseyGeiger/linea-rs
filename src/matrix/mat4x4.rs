use std::ops::*;
use std::convert::*;
#[allow(unused_imports)]
use std::{f32, f64};
use ::matrix::mat3x3::*;
use ::vector::vec3::*;
use ::vector::vec4::*;

#[derive(Copy, Clone)]
pub struct Mat4x4<T: Copy + Clone> {
    ///Array of Vector4<T> to facilitate the writing of operations
    pub m: [Vector4<T>; 4],
}

#[allow(dead_code)]
impl<T: Copy + Clone> Mat4x4<T> {

    pub fn new(m00: T, m01: T, m02: T, m03: T, m10: T, m11: T, m12: T, m13: T, m20: T, m21: T, m22: T, m23: T, m30: T, m31: T, m32: T, m33: T) -> Mat4x4<T> {
        Mat4x4 {
            m: [
                Vector4 {x: m00, y: m01, z: m02, w: m03},
                Vector4 {x: m10, y: m11, z: m12, w: m13},
                Vector4 {x: m20, y: m21, z: m22, w: m23},
                Vector4 {x: m30, y: m31, z: m32, w: m33}
            ]
        }
    }

    pub fn from_vectors(m0: Vector4<T>, m1: Vector4<T>, m2: Vector4<T>, m3: Vector4<T>) -> Mat4x4<T> {
        Mat4x4 {
            m: [m0, m1, m2, m3]
        }
    }

    ///Returns the dimensions of this matrix (for generic usage)
    pub fn dim(&self) -> (usize, usize) {
        (4, 4)
    }

    ///Returns a Mat4x4 which has been transposed, i.e. m[0][1] and m[1][0] swapped for Mat4x4
    pub fn as_transposed(&self) -> Mat4x4<T> {
        let mut m = Mat4x4 {
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

    pub fn submatrix(&self, i: usize, j: usize) -> Mat3x3<T> {
        let col_0 = (i+1) % self.dim().0;
        let col_1 = (i+2) % self.dim().0;
        let col_2 = (i+3) % self.dim().0;

        let row_0 = (j+1) % self.dim().1;
        let row_1 = (j+2) % self.dim().1;
        let row_2 = (j+3) % self.dim().1;

        let (col_0, col_1, col_2) = if col_1 < col_0 {
            if col_2 < col_1 {
                (col_2, col_1, col_0)
            } else {
                if col_0 < col_2 {
                    (col_1, col_0, col_2)
                } else {
                    (col_1, col_2, col_0)
                }
            }
        } else {
            if col_2 < col_0 {
                (col_2, col_0, col_1)
            } else {
                if col_2 < col_1 {
                    (col_0, col_2, col_1)
                } else {
                    (col_0, col_1, col_2)
                }
            }
        };

        let (row_0, row_1, row_2) = if row_1 < row_0 {
            if row_2 < row_1 {
                (row_2, row_1, row_0)
            } else {
                if row_0 < row_2 {
                    (row_1, row_0, row_2)
                } else {
                    (row_1, row_2, row_0)
                }
            }
        } else {
            if row_2 < row_0 {
                (row_2, row_0, row_1)
            } else {
                if row_2 < row_1 {
                    (row_0, row_2, row_1)
                } else {
                    (row_0, row_1, row_2)
                }
            }
        };

        let m0: Vector3<T> = Vector3::new(self[col_0][row_0], self[col_0][row_1], self[col_0][row_2]);
        let m1: Vector3<T> = Vector3::new(self[col_1][row_0], self[col_1][row_1], self[col_1][row_2]);
        let m2: Vector3<T> = Vector3::new(self[col_2][row_0], self[col_2][row_1], self[col_2][row_2]);
        let sub_matrix: Mat3x3<T> = Mat3x3::from_vectors(m0, m1, m2);

        sub_matrix
    }

}

impl<T: Copy + Clone> Index<usize> for Mat4x4<T> {
    type Output = Vector4<T>;

    fn index(&self, idx: usize) -> &Vector4<T> {
        &self.m[idx]
    }
}

impl<T: Copy + Clone> Index<(usize, usize)> for Mat4x4<T> {
    type Output = T;

    fn index(&self, idx: (usize, usize)) -> &T {
        return &self.m[idx.0][idx.1];
    }
}

impl<T: Copy + Clone> IndexMut<usize> for Mat4x4<T> {
    fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut Vector4<T> {
        &mut self.m[idx]
    }
}

impl<T: Copy + Clone> IndexMut<(usize, usize)> for Mat4x4<T> {
    fn index_mut<'a>(&'a mut self, idx: (usize, usize)) -> &'a mut T {
        &mut self.m[idx.0][idx.1]
    }
}

impl<T: Copy + Clone + Add<Output = T>> Add for Mat4x4<T> {
    type Output = Mat4x4<T>;

    fn add(self, rhs: Self) -> Mat4x4<T> {
        Mat4x4 {
            m: [
                self[0] + rhs[0],
                self[1] + rhs[1],
                self[2] + rhs[2],
                self[3] + rhs[3]
            ]
        }
    }
}

impl<T: Copy + Clone + Add<Output = T>> AddAssign for Mat4x4<T> {
    fn add_assign(&mut self, rhs: Self) {
        self[0] = self[0] + rhs[0];
        self[1] = self[1] + rhs[1];
        self[2] = self[2] + rhs[2];
        self[3] = self[3] + rhs[3];
    }
}

impl<T: Copy + Clone + Sub<Output = T>> Sub for Mat4x4<T> {
    type Output = Mat4x4<T>;

    fn sub(self, rhs: Self) -> Mat4x4<T> {
        Mat4x4 {
            m: [
                self[0] - rhs[0],
                self[1] - rhs[1],
                self[2] - rhs[2],
                self[3] - rhs[3]
            ]
        }
    }
}

impl<T: Copy + Clone + Sub<Output = T>> SubAssign for Mat4x4<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self[0] = self[0] - rhs[0];
        self[1] = self[1] - rhs[1];
        self[2] = self[2] - rhs[2];
        self[3] = self[3] - rhs[3];
    }
}

impl<T: Copy + Clone + Mul<Output = T>> Mul<T> for Mat4x4<T> {
    type Output = Mat4x4<T>;

    fn mul(self, rhs: T) -> Mat4x4<T> {
        Mat4x4 {
            m: [
                self.m[0] * rhs,
                self.m[1] * rhs,
                self.m[2] * rhs,
                self.m[3] * rhs
            ]
        }
    }
}

impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T>> Mul<Vector4<T>> for Mat4x4<T> {
    type Output = Vector4<T>;

    fn mul(self, rhs: Vector4<T>) -> Vector4<T> {
        let self_transposed = self.as_transposed();

        Vector4 {
            x: self_transposed[0].dot(rhs),
            y: self_transposed[1].dot(rhs),
            z: self_transposed[2].dot(rhs),
            w: self_transposed[3].dot(rhs)
        }
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> Mul for Mat4x4<T> {
    type Output = Mat4x4<T>;

    fn mul(self, rhs: Self) -> Mat4x4<T> {
        //Get a transposed version of self to do dot(row, column) for multiplication
        let self_transposed = self.as_transposed();

        let mut result: Mat4x4<T> = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                result[(i, j)] = self_transposed.m[i].dot(rhs.m[j]);
            }
        }

        result
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> MulAssign<T> for Mat4x4<T> {
    fn mul_assign(&mut self, rhs: T) {
        self[0] = self[0] * rhs;
        self[1] = self[1] * rhs;
        self[2] = self[2] * rhs;
        self[3] = self[3] * rhs;
    }
}

impl<T: Copy + Clone + Add<Output = T> + Mul<Output = T>> MulAssign for Mat4x4<T> {
    fn mul_assign(&mut self, rhs: Self) {
        let self_transposed = self.as_transposed();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                self[(i, j)] = self_transposed.m[i].dot(rhs.m[j]);
            }
        }
    }
}

impl<T: Copy + Clone + Div<Output = T>> Div<T> for Mat4x4<T> {
    type Output = Mat4x4<T>;

    fn div(self, rhs: T) -> Mat4x4<T> {
        Mat4x4 {
            m: [
                self[0] / rhs,
                self[1] / rhs,
                self[2] / rhs,
                self[3] / rhs
            ]
        }
    }
}

impl<T: Copy + Clone + Div<Output = T>> DivAssign<T> for Mat4x4<T> {
    fn div_assign(&mut self, rhs: T) {
        self[0] = self[0] / rhs;
        self[1] = self[1] / rhs;
        self[2] = self[2] / rhs;
        self[3] = self[3] / rhs;
    }
}

impl<T: Copy + Clone + Neg<Output = T>> Neg for Mat4x4<T> {
    type Output = Mat4x4<T>;

    fn neg(self) -> Mat4x4<T> {
        Mat4x4 {
            m: [
                -self[0],
                -self[1],
                -self[2],
                -self[3]
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Sub<Output = T>> Mat4x4<T> {
    pub fn from_diagonal(diag: Vector4<T>) -> Mat4x4<T> {
        let zero = diag.x - diag.x;

        Mat4x4 {
            m: [
                Vector4 {x: diag.x, y: zero,   z: zero,   w: zero  },
                Vector4 {x: zero,   y: diag.y, z: zero,   w: zero  },
                Vector4 {x: zero,   y: zero,   z: diag.z, w: zero  },
                Vector4 {x: zero,   y: zero,   z: zero,   w: diag.w}
            ]
        }
    }

    pub fn scale(s: T) -> Mat4x4<T> {
        let zero: T = s - s;

        Mat4x4 {
            m: [
                Vector4 {x: s, y: zero, z: zero, w: zero},
                Vector4 {x: zero, y: s, z: zero, w: zero},
                Vector4 {x: zero, y: zero, z: s, w: zero},
                Vector4 {x: zero, y: zero, z: zero, w: s}
            ]
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + From<u8>> Mat4x4<T> {
    pub fn identity() -> Mat4x4<T> {
        let one: T = From::from(1u8);
        let zero: T = From::from(0u8);

        Mat4x4 {
            m: [
                Vector4 {x: one, y: zero, z: zero, w: zero},
                Vector4 {x: zero, y: one, z: zero, w: zero},
                Vector4 {x: zero, y: zero, z: one, w: zero},
                Vector4 {x: zero, y: zero, z: zero, w: one}
            ]
        }
    }

    pub fn translate(pos: Vector3<T>) -> Mat4x4<T> {
        let mut t = Self::identity();
        t[3] = Vector4 {x: pos[0], y: pos[1], z: pos[2], w: From::from(1u8)};

        t
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Neg<Output = T>> Mat4x4<T> {
    pub fn det(&self) -> T {
        let mut _det: T = self[0][0] - self[0][0];

        for i in 0..self.dim().0 {
            let m = self.submatrix(i, 0);

            _det = if (i % 2) == 1 {_det + (m.det() * self[i][0])} else {_det - (m.det() * self[i][0])};
        }

        return _det
    }

    pub fn minor(&self) -> Mat4x4<T> {
        let mut minor = self.clone();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                let sub = self.submatrix(i, j);

                minor[i][j] = sub.det();
            }
        }

        minor
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Sub<Output = T>> Mat4x4<T> {
    pub fn adjugate(&self) -> Mat4x4<T> {
        let mut adj = self.minor();

        for i in 0..self.dim().0 {
            for j in 0..self.dim().1 {
                let negative: bool = (i + j) % 2 == 0;

                if negative {
                    adj[i][j] = -adj[i][j];
                }
            }
        }

        adj.transpose();

        adj
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Neg<Output = T> + PartialEq> Mat4x4<T> {
    pub fn inverse(&self) -> Mat4x4<T> {
        let zero: T = self[0][0] - self[0][0];
        let determinate = self.det();

        if determinate == zero {
            return Mat4x4 {
                m: [
                    Vector4 {x: zero, y: zero, z: zero, w: zero},
                    Vector4 {x: zero, y: zero, z: zero, w: zero},
                    Vector4 {x: zero, y: zero, z: zero, w: zero},
                    Vector4 {x: zero, y: zero, z: zero, w: zero}
                ]
            }
        } else {
            return self.adjugate() / determinate;
        }
    }
}

#[allow(dead_code)]
impl Mat4x4<f32> {

    pub fn pitch(angle_rad: f32) -> Mat4x4<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat4x4 {
            m: [
                Vector4 {x: 1.0f32, y: 0.0f32, z: 0.0f32, w: 0.0f32},
                Vector4 {x: 0.0f32, y: cos,    z: sin,    w: 0.0f32},
                Vector4 {x: 0.0f32, y: -sin,   z: cos,    w: 0.0f32},
                Vector4 {x: 0.0f32, y: 0.0f32, z: 0.0f32, w: 1.0f32}
            ]
        }
    }

    pub fn yaw(angle_rad: f32) -> Mat4x4<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat4x4 {
            m: [
                Vector4 {x: cos,    y: 0.0f32, z: -sin,   w: 0.0f32},
                Vector4 {x: 0.0f32, y: 1.0f32, z: 0.0f32, w: 0.0f32},
                Vector4 {x: sin,    y: 0.0f32, z: cos,    w: 0.0f32},
                Vector4 {x: 0.0f32, y: 0.0f32, z: 0.0f32, w: 1.0f32}
            ]
        }
    }

    pub fn roll(angle_rad: f32) -> Mat4x4<f32> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat4x4 {
            m: [
                Vector4 {x: cos,    y: sin,    z: 0.0f32, w: 0.0f32},
                Vector4 {x: -sin,   y: cos,    z: 0.0f32, w: 0.0f32},
                Vector4 {x: 0.0f32, y: 0.0f32, z: 1.0f32, w: 0.0f32},
                Vector4 {x: 0.0f32, y: 0.0f32, z: 0.0f32, w: 1.0f32}
            ]
        }
    }

    pub fn angle_axis(angle_rad: f32, axis: Vector4<f32>) -> Mat4x4<f32> {
        let axis = axis.norm();

        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        let x_sqr = axis.x * axis.x;
        let y_sqr = axis.y * axis.y;
        let z_sqr = axis.z * axis.z;

        let versine = 1.0f32 - cos;

        let x_basis: Vector4<f32> = Vector4 {
            x: cos + (x_sqr * versine),
            y: (axis.y * axis.x * versine) + (axis.z * sin),
            z: (axis.z * axis.x * versine) - (axis.y * sin),
            w: 0.0f32
        };

        let y_basis: Vector4<f32> = Vector4 {
            x: (axis.x * axis.y * versine) - (axis.z * sin),
            y: cos + (y_sqr * versine),
            z: (axis.z * axis.y * versine) + (axis.x * sin),
            w: 0.0f32
        };

        let z_basis: Vector4<f32> = Vector4 {
            x: (axis.x * axis.z * versine) + (axis.y * sin),
            y: (axis.y * axis.z * versine) - (axis.x * sin),
            z: cos + (z_sqr * versine),
            w: 0.0f32
        };

        let w_basis: Vector4<f32> = Vector4 {
            x: 0.0f32,
            y: 0.0f32,
            z: 0.0f32,
            w: 1.0f32
        };

        Mat4x4 {
            m: [
                x_basis,
                y_basis,
                z_basis,
                w_basis
            ]
        }
    }

    fn perspective_fov(field_of_view_degree: f32, aspect_ratio: f32, near: f32, far: f32) -> Mat4x4<f32> {
        let cot_half_fov = 1.0f32 / (field_of_view_degree / 2.0f32).to_radians().tan();
        let depth = far - near;
        let n_f_sum = near + far;
        let n_f_2 = 2.0f32 * near * far;

        let x_basis = Vector4 { x: cot_half_fov / aspect_ratio, y: 0.0f32,       z: 0.0f32,        w:  0.0f32};
        let y_basis = Vector4 { x: 0.0f32,                      y: cot_half_fov, z: 0.0f32,        w:  0.0f32};
        let z_basis = Vector4 { x: 0.0f32,                      y: 0.0f32,       z: n_f_sum/depth, w: -1.0f32};
        let w_basis = Vector4 { x: 0.0f32,                      y: 0.0f32,       z: n_f_2 / depth, w:  0.0f32};

        Mat4x4 {
            m: [
                x_basis,
                y_basis,
                z_basis,
                w_basis
            ]
        }
    }

    fn orthographic(left: f32, right: f32, bottom: f32, top: f32, near: f32, far: f32) -> Mat4x4<f32> {
        let width = right - left;
        let height = top - bottom;
        let depth = far - near;
        let tx = -(right + left)/width;
        let ty = -(top + bottom)/height;
        let tz = -(far + near)/depth;

        let x_basis = Vector4 { x: 2.0f32 / width, y: 0.0f32,          z: 0.0f32,          w: 0.0f32};
        let y_basis = Vector4 { x: 0.0f32,         y: 2.0f32 / height, z: 0.0f32,          w: 0.0f32};
        let z_basis = Vector4 { x: 0.0f32,         y: 0.0f32,          z: -2.0f32 / depth, w: 0.0f32};
        let w_basis = Vector4 { x: tx,             y: ty,              z: tz,              w: 1.0f32};

        Mat4x4 {
            m: [
                x_basis,
                y_basis,
                z_basis,
                w_basis
            ]
        }
    }

    fn viewport(left: f32, right: f32, bottom: f32, top: f32) -> Mat4x4<f32> {
        let width = right - left;
        let height = top - bottom;

        Mat4x4 {
            m: [
                Vector4::new(width/2.0f32,        0.0f32,              0.0f32, 0.0f32),
                Vector4::new(0.0f32,              height/2.0f32,       0.0f32, 0.0f32),
                Vector4::new(0.0f32,              0.0f32,              0.5f32, 0.0f32),
                Vector4::new((left+right)/2.0f32, (bottom+top)/2.0f32, 0.5f32, 1.0f32)
            ]
        }
    }
}

#[allow(dead_code)]
impl Mat4x4<f64> {
    pub fn pitch(angle_rad: f64) -> Mat4x4<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat4x4 {
            m: [
                Vector4 {x: 1.0f64, y: 0.0f64, z: 0.0f64, w: 0.0f64},
                Vector4 {x: 0.0f64, y: cos,    z: sin,    w: 0.0f64},
                Vector4 {x: 0.0f64, y: -sin,   z: cos,    w: 0.0f64},
                Vector4 {x: 0.0f64, y: 0.0f64, z: 0.0f64, w: 1.0f64}
            ]
        }
    }

    pub fn yaw(angle_rad: f64) -> Mat4x4<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat4x4 {
            m: [
                Vector4 {x: cos,    y: 0.0f64, z: -sin,    w: 0.0f64},
                Vector4 {x: 0.0f64, y: 1.0f64, z: 0.0f64,  w: 0.0f64},
                Vector4 {x: sin,    y: 0.0f64, z: cos,     w: 0.0f64},
                Vector4 {x: 0.0f64, y: 0.0f64, z: 0.0f64,  w: 1.0f64}
            ]
        }
    }

    pub fn roll(angle_rad: f64) -> Mat4x4<f64> {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        Mat4x4 {
            m: [
                Vector4 {x: cos,    y: sin,    z: 0.0f64, w: 0.0f64},
                Vector4 {x: -sin,   y: cos,    z: 0.0f64, w: 0.0f64},
                Vector4 {x: 0.0f64, y: 0.0f64, z: 1.0f64, w: 0.0f64},
                Vector4 {x: 0.0f64, y: 0.0f64, z: 0.0f64, w: 1.0f64}
            ]
        }
    }

    pub fn angle_axis(angle_rad: f64, axis: Vector4<f64>) -> Mat4x4<f64> {
        let axis = axis.norm();

        let cos = angle_rad.cos();
        let sin = angle_rad.sin();

        let x_sqr = axis.x * axis.x;
        let y_sqr = axis.y * axis.y;
        let z_sqr = axis.z * axis.z;

        let versine = 1.0f64 - cos;

        let x_basis: Vector4<f64> = Vector4 {
            x: cos + (x_sqr * versine),
            y: (axis.y * axis.x * versine) + (axis.z * sin),
            z: (axis.z * axis.x * versine) - (axis.y * sin),
            w: 0.0f64
        };

        let y_basis: Vector4<f64> = Vector4 {
            x: (axis.x * axis.y * versine) - (axis.z * sin),
            y: cos + (y_sqr * versine),
            z: (axis.z * axis.y * versine) + (axis.x * sin),
            w: 0.0f64
        };

        let z_basis: Vector4<f64> = Vector4 {
            x: (axis.x * axis.z * versine) + (axis.y * sin),
            y: (axis.y * axis.z * versine) - (axis.x * sin),
            z: cos + (z_sqr * versine),
            w: 0.0f64
        };

        let w_basis: Vector4<f64> = Vector4 {
            x: 0.0f64,
            y: 0.0f64,
            z: 0.0f64,
            w: 1.0f64
        };

        Mat4x4 {
            m: [
                x_basis,
                y_basis,
                z_basis,
                w_basis
            ]
        }
    }

    fn perspective_fov(field_of_view_degree: f64, aspect_ratio: f64, near: f64, far: f64) -> Mat4x4<f64> {
        let cot_half_fov = 1.0f64 / (field_of_view_degree / 2.0f64).to_radians().tan();
        let depth = far - near;
        let n_f_sum = near + far;
        let n_f_2 = 2.0f64 * near * far;

        let x_basis = Vector4 { x: cot_half_fov / aspect_ratio, y: 0.0f64,       z: 0.0f64,        w:  0.0f64};
        let y_basis = Vector4 { x: 0.0f64,                      y: cot_half_fov, z: 0.0f64,        w:  0.0f64};
        let z_basis = Vector4 { x: 0.0f64,                      y: 0.0f64,       z: n_f_sum/depth, w: -1.0f64};
        let w_basis = Vector4 { x: 0.0f64,                      y: 0.0f64,       z: n_f_2 / depth, w:  0.0f64};

        Mat4x4 {
            m: [
                x_basis,
                y_basis,
                z_basis,
                w_basis
            ]
        }
    }

    fn orthographic(left: f64, right: f64, bottom: f64, top: f64, near: f64, far: f64) -> Mat4x4<f64> {
        let width = right - left;
        let height = top - bottom;
        let depth = far - near;
        let tx = -(right + left)/width;
        let ty = -(top + bottom)/height;
        let tz = -(far + near)/depth;

        let x_basis = Vector4 { x: 2.0f64 / width, y: 0.0f64,          z: 0.0f64,          w: 0.0f64};
        let y_basis = Vector4 { x: 0.0f64,         y: 2.0f64 / height, z: 0.0f64,          w: 0.0f64};
        let z_basis = Vector4 { x: 0.0f64,         y: 0.0f64,          z: -2.0f64 / depth, w: 0.0f64};
        let w_basis = Vector4 { x: tx,             y: ty,              z: tz,              w: 1.0f64};

        Mat4x4 {
            m: [
                x_basis,
                y_basis,
                z_basis,
                w_basis
            ]
        }
    }

    fn viewport(left: f64, right: f64, bottom: f64, top: f64) -> Mat4x4<f64> {
        let width = right - left;
        let height = top - bottom;

        Mat4x4 {
            m: [
                Vector4::new(width/2.0f64,        0.0f64,              0.0f64, 0.0f64),
                Vector4::new(0.0f64,              height/2.0f64,       0.0f64, 0.0f64),
                Vector4::new(0.0f64,              0.0f64,              0.5f64, 0.0f64),
                Vector4::new((left+right)/2.0f64, (bottom+top)/2.0f64, 0.5f64, 1.0f64)
            ]
        }
    }
}

type Mat4x4f = Mat4x4<f32>;
type Mat4x4d = Mat4x4<f64>;

type Mat4x4ll = Mat4x4<i64>;
type Mat4x4i = Mat4x4<i32>;
type Mat4x4s = Mat4x4<i16>;
type Mat4x4b = Mat4x4<i8>;

type Mat4x4ull = Mat4x4<u64>;
type Mat4x4ui = Mat4x4<u32>;
type Mat4x4us = Mat4x4<u16>;
type Mat4x4ub = Mat4x4<u8>;
