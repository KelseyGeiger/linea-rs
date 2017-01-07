use std::ops::*;
use std::fmt::*;
use std::{f32, f64};
use vector::vec2::*;
use vector::vec3::*;

#[derive(Copy, Clone, Debug)]
pub struct Vector4<T: Copy + Clone> {
	pub x: T,
	pub y: T,
    pub z: T,
    pub w: T,
}

#[allow(dead_code)]
impl<T: Copy + Clone> Vector4<T> {

	pub fn new(x: T, y: T, z: T, w: T) -> Vector4<T> {
		Vector4 {
			x: x,
			y: y,
            z: z,
            w: w
		}
	}

    pub fn from_vec3(v3: Vector3<T>, w: T) -> Vector4<T> {
        Vector4 {
            x: v3.x,
            y: v3.y,
            z: v3.z,
            w: w
        }
    }

    pub fn from_vec2(v2: Vector2<T>, z: T, w: T) -> Vector4<T> {
        Vector4 {
            x: v2.x,
            y: v2.y,
            z: z,
            w: w
        }
    }

    pub fn from_xy_zw(xy: Vector2<T>, zw: Vector2<T>) -> Vector4<T> {
        Vector4 {
            x: xy.x,
            y: xy.y,
            z: zw.x,
            w: zw.y
        }
    }


    /*
     *  NOTE: Swizzle functions below- very repetitive and prone to error.
     *  These were not generated.
    */


    pub fn xy(self) -> Vector2<T> {
        Vector2 {
            x: self.x,
            y: self.y
        }
    }

    pub fn xz(self) -> Vector2<T> {
        Vector2 {
            x: self.x,
            y: self.z
        }
    }

    pub fn xw(self) -> Vector2<T> {
        Vector2 {
            x: self.x,
            y: self.w
        }
    }

    pub fn yx(self) -> Vector2<T> {
        self.xy().yx()
    }

    pub fn yz(self) -> Vector2<T> {
        Vector2 {
            x: self.y,
            y: self.z
        }
    }

    pub fn yw(self) -> Vector2<T> {
        Vector2 {
            x: self.y,
            y: self.w
        }
    }

    pub fn zx(self) -> Vector2<T> {
        self.xz().yx()
    }

    pub fn zy(self) -> Vector2<T> {
        self.yz().yx()
    }

    pub fn zw(self) -> Vector2<T> {
        Vector2 {
            x: self.z,
            y: self.w
        }
    }

    pub fn wx(self) -> Vector2<T> {
        self.xw().yx()
    }

    pub fn wy(self) -> Vector2<T> {
        self.yw().yx()
    }

    pub fn wz(self) -> Vector2<T> {
        self.zw().yx()
    }

    pub fn xyz(self) -> Vector3<T> {
        Vector3 {
            x: self.x,
            y: self.y,
            z: self.z
        }
    }

    pub fn xyw(self) -> Vector3<T> {
        Vector3 {
            x: self.x,
            y: self.y,
            z: self.w
        }
    }

    pub fn xzy(self) -> Vector3<T> {
        Vector3 {
            x: self.x,
            y: self.z,
            z: self.y
        }
    }

    pub fn xzw(self) -> Vector3<T> {
        Vector3 {
            x: self.x,
            y: self.z,
            z: self.w
        }
    }

    pub fn xwy(self) -> Vector3<T> {
        self.xyw().xzy()
    }

    pub fn xwz(self) -> Vector3<T> {
        self.xzw().xzy()
    }

    pub fn yxz(self) -> Vector3<T> {
        Vector3 {
            x: self.y,
            y: self.x,
            z: self.z
        }
    }

    pub fn yxw(self) -> Vector3<T> {
        Vector3 {
            x: self.y,
            y: self.x,
            z: self.w
        }
    }

    pub fn yzx(self) -> Vector3<T> {
        Vector3 {
            x: self.y,
            y: self.z,
            z: self.x
        }
    }

    pub fn yzw(self) -> Vector3<T> {
        Vector3 {
            x: self.y,
            y: self.z,
            z: self.w
        }
    }

    pub fn ywx(self) -> Vector3<T> {
        self.xyw().yzx()
    }

    pub fn ywz(self) -> Vector3<T> {
        self.yzw().xzy()
    }

    pub fn zxy(self) -> Vector3<T> {
        Vector3 {
            x: self.z,
            y: self.x,
            z: self.y
        }
    }

    pub fn zxw(self) -> Vector3<T> {
        Vector3 {
            x: self.z,
            y: self.x,
            z: self.w
        }
    }

    pub fn zyx(self) -> Vector3<T> {
        self.zxy().xzy()
    }

    pub fn zyw(self) -> Vector3<T> {
        Vector3 {
            x: self.z,
            y: self.y,
            z: self.w
        }
    }

    pub fn zwx(self) -> Vector3<T> {
        self.zxw().xzy()
    }

    pub fn zwy(self) -> Vector3<T> {
        self.zyw().xzy()
    }

    pub fn wxy(self) -> Vector3<T> {
        self.xyw().zxy()
    }

    pub fn wxz(self) -> Vector3<T> {
        self.xzw().zxy()
    }

    pub fn wyx(self) -> Vector3<T> {
        self.xyw().zyx()
    }

    pub fn wyz(self) -> Vector3<T> {
        self.yzw().zxy()
    }

    pub fn wzx(self) -> Vector3<T> {
        self.xzw().zyx()
    }

    pub fn wzy(self) -> Vector3<T> {
        self.yzw().zyx()
    }

    pub fn xywz(self) -> Vector4<T> {
        Vector4 {
            x: self.x,
            y: self.y,
            z: self.w,
            w: self.z
        }
    }

    pub fn xzyw(self) -> Vector4<T> {
        Vector4 {
            x: self.x,
            y: self.z,
            z: self.y,
            w: self.w
        }
    }

    pub fn xzwy(self) -> Vector4<T> {
        Vector4 {
            x: self.x,
            y: self.z,
            z: self.w,
            w: self.y
        }
    }

    pub fn xwyz(self) -> Vector4<T> {
        Vector4 {
            x: self.x,
            y: self.w,
            z: self.y,
            w: self.z
        }
    }

    pub fn xwzy(self) -> Vector4<T> {
        Vector4 {
            x: self.x,
            y: self.w,
            z: self.z,
            w: self.y
        }
    }

    pub fn yxzw(self) -> Vector4<T> {
        Vector4 {
            x: self.y,
            y: self.x,
            z: self.z,
            w: self.w
        }
    }

    pub fn yxwz(self) -> Vector4<T> {
        Vector4 {
            x: self.y,
            y: self.x,
            z: self.w,
            w: self.z
        }
    }

    pub fn yzxw(self) -> Vector4<T> {
        Vector4 {
            x: self.y,
            y: self.z,
            z: self.x,
            w: self.w
        }
    }

    pub fn yzwx(self) -> Vector4<T> {
        Vector4 {
            x: self.y,
            y: self.z,
            z: self.w,
            w: self.x
        }
    }

    pub fn ywxz(self) -> Vector4<T> {
        Vector4 {
            x: self.y,
            y: self.w,
            z: self.x,
            w: self.z
        }
    }

    pub fn ywzx(self) -> Vector4<T> {
        Vector4 {
            x: self.y,
            y: self.w,
            z: self.z,
            w: self.x
        }
    }

    pub fn zxyw(self) -> Vector4<T> {
        Vector4 {
            x: self.z,
            y: self.y,
            z: self.y,
            w: self.w
        }
    }

    pub fn zxwy(self) -> Vector4<T> {
        Vector4 {
            x: self.z,
            y: self.x,
            z: self.w,
            w: self.y
        }
    }

    pub fn zyxw(self) -> Vector4<T> {
        Vector4 {
            x: self.z,
            y: self.y,
            z: self.x,
            w: self.w
        }
    }

    pub fn zywx(self) -> Vector4<T> {
        Vector4 {
            x: self.z,
            y: self.y,
            z: self.w,
            w: self.x
        }
    }

    pub fn zwxy(self) -> Vector4<T> {
        Vector4 {
            x: self.z,
            y: self.w,
            z: self.x,
            w: self.y
        }
    }

    pub fn zwyx(self) -> Vector4<T> {
        Vector4 {
            x: self.z,
            y: self.w,
            z: self.y,
            w: self.x
        }
    }

    pub fn wxyz(self) -> Vector4<T> {
        Vector4 {
            x: self.w,
            y: self.x,
            z: self.y,
            w: self.z
        }
    }

    pub fn wxzy(self) -> Vector4<T> {
        Vector4 {
            x: self.w,
            y: self.x,
            z: self.z,
            w: self.y
        }
    }

    pub fn wyxz(self) -> Vector4<T> {
        Vector4 {
            x: self.w,
            y: self.y,
            z: self.x,
            w: self.z
        }
    }

    pub fn wyzx(self) -> Vector4<T> {
        Vector4 {
            x: self.w,
            y: self.y,
            z: self.z,
            w: self.x
        }
    }

    pub fn wzxy(self) -> Vector4<T> {
        Vector4 {
            x: self.w,
            y: self.z,
            z: self.x,
            w: self.y
        }
    }

    pub fn wzyx(self) -> Vector4<T> {
        Vector4 {
            x: self.w,
            y: self.z,
            z: self.y,
            w: self.x
        }
    }
}

impl<T: Copy + Clone + Debug> Index<usize> for Vector4<T> {
	type Output = T;

	fn index(&self, idx: usize) -> &T {
		if idx >= 4 {
			panic!("Index {} out of bounds for {:?}", idx, *self);
		} else if idx == 0 {
			return &self.x;
		} else if idx == 1{
			return &self.y;
		} else if idx == 2 {
            return &self.z;
        } else {
            return &self.w;
        }
	}

}

impl<T: Copy + Clone + Debug> IndexMut<usize> for Vector4<T> {

	fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut T {
		if idx >= 4 {
			panic!("Index {} out of bounds for {:?}", idx, *self);
		} else if idx == 0 {
			return &mut self.x;
		} else if idx == 1 {
			return &mut self.y;
		} else if idx == 2 {
            return &mut self.z;
        } else {
            return &mut self.w;
        }
	}

}

impl<T: Copy + Clone + Add<Output = T>> Add for Vector4<T> {
	type Output = Vector4<T>;

	fn add(self, other: Self) -> Vector4<T> {
		Vector4 {
			x: self.x + other.x,
			y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w
		}
	}

}

impl<T: Copy + Clone + Add<Output = T>> AddAssign for Vector4<T> {

	fn add_assign(&mut self, other: Self) {
		self.x = self.x + other.x;
		self.y = self.y + other.y;
        self.z = self.z + other.z;
        self.w = self.w + other.w;
	}

}

impl<T: Copy + Clone + Sub<Output = T>> Sub for Vector4<T> {
	type Output = Vector4<T>;

	fn sub(self, other: Self) -> Vector4<T> {
		Vector4 {
			x: self.x - other.x,
			y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w
		}
	}

}

impl<T: Copy + Clone + Sub<Output = T>> SubAssign for Vector4<T> {

	fn sub_assign(&mut self, other: Self) {
		self.x = self.x - other.x;
		self.y = self.y - other.y;
        self.z = self.z - other.z;
        self.w = self.w - other.w;
	}

}

impl<T: Copy + Clone + Mul<Output = T>> Mul<T> for Vector4<T> {
	type Output = Vector4<T>;

	fn mul(self, rhs: T) -> Vector4<T> {
		Vector4 {
			x: self.x * rhs,
			y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs
		}
	}

}

impl<T: Copy + Clone + Mul<Output = T>> Mul<Vector4<T>> for Vector4<T> {
	type Output = Vector4<T>;

	fn mul(self, rhs: Vector4<T>) -> Vector4<T> {
		Vector4 {
			x: self.x * rhs.x,
			y: self.y * rhs.y,
            z: self.z * rhs.z,
            w: self.w * rhs.w
		}
	}

}

impl<T: Copy + Clone + Mul<Output = T>> MulAssign<T> for Vector4<T> {

	fn mul_assign(&mut self, rhs: T) {
		self.x = self.x * rhs;
		self.y = self.y * rhs;
        self.z = self.z * rhs;
        self.w = self.w * rhs;
	}

}

impl<T: Copy + Clone + Mul<Output = T>> MulAssign for Vector4<T> {

	fn mul_assign(&mut self, rhs: Vector4<T>) {
		self.x = self.x * rhs.x;
		self.y = self.y * rhs.y;
        self.z = self.z * rhs.z;
        self.w = self.w * rhs.w;
	}

}

impl<T: Copy + Clone + Div<Output = T>> Div<T> for Vector4<T> {
	type Output = Vector4<T>;

	fn div(self, rhs: T) -> Vector4<T> {
		Vector4 {
			x: self.x / rhs,
			y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs
		}
	}

}

impl<T: Copy + Clone + Div<Output = T>> DivAssign<T> for Vector4<T> {

	fn div_assign(&mut self, rhs: T) {
		self.x = self.x / rhs;
		self.y = self.y / rhs;
        self.z = self.z / rhs;
        self.w = self.w / rhs;
	}
}

impl<T: Copy + Clone + Neg<Output = T>> Neg for Vector4<T> {
	type Output = Vector4<T>;

	fn neg(self) -> Vector4<T> {
		Vector4 {
			x: -self.x,
			y: -self.y,
            z: -self.z,
            w: -self.w
		}
	}

}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T>> Vector4<T> {

	pub fn dot(self,  r: Vector4<T>) -> T {
		let product: Vector4<T> = self * r;
		let dot: T = product.x + product.y + product.z + product.w;

		dot
	}

}

#[allow(dead_code)]
impl<T: Copy + Clone + Into<f64> + From<f64>> Vector4<T> {

	pub fn mag(self) -> f64 {
		let x: f64 = self.x.into();
		let y: f64 = self.y.into();
        let z: f64 = self.z.into();
        let w: f64 = self.w.into();

		((x * x) + (y * y) + (z * z) + (w * w)).sqrt()
	}

	pub fn proj(self, rhs: Vector4<T>) -> Vector4<T> {
		let l_converted: Vector4<f64> = Vector4::new(self.x.into(), self.y.into(), self.z.into(), self.w.into());
		let r_converted: Vector4<f64> = Vector4::new(rhs.x.into(), rhs.y.into(), rhs.z.into(), rhs.w.into());

		let r_normal = r_converted / r_converted.mag();
		let projected = r_normal * (l_converted.dot(r_converted));

		Vector4 {
			x: projected.x.into(),
			y: projected.y.into(),
            z: projected.z.into(),
            w: projected.w.into()
		}
	}
}

#[allow(dead_code)]
impl Vector4<f64> {

	pub fn norm(self) -> Vector4<f64> {

		let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z) + (self.w * self.w);
		let mag = mag_sqr.sqrt();

		let vec_scaled = self / mag;

		vec_scaled
	}

	pub fn normalize(&mut self) {
		let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z) + (self.w * self.w);
		let mag = mag_sqr.sqrt();

		self.x = self.x / mag;
		self.y = self.y / mag;
        self.z = self.z / mag;
        self.w = self.w / mag;
	}

}

#[allow(dead_code)]
impl Vector4<f32> {

	pub fn norm(self) -> Vector4<f32> {

        let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z) + (self.w * self.w);
		let mag = mag_sqr.sqrt();

		let vec_scaled = self / mag;

		vec_scaled
	}

	pub fn normalize(&mut self) {
		let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z) + (self.w * self.w);
		let mag = mag_sqr.sqrt();

		self.x = self.x / mag;
		self.y = self.y / mag;
        self.z = self.z / mag;
        self.w = self.w / mag;
	}

}

pub type Vec4d = Vector4<f64>;
pub type Vec4 = Vector4<f32>;
pub type Vec4f = Vector4<f32>;

pub type Vec4ll = Vector4<i64>;
pub type Vec4i = Vector4<i32>;
pub type Vec4s = Vector4<i16>;
pub type Vec4b = Vector4<i8>;

pub type Vec4ull = Vector4<u64>;
pub type Vec4ui = Vector4<u32>;
pub type Vec4us = Vector4<u16>;
pub type Vec4ub = Vector4<u8>;
