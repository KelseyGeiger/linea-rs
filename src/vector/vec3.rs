use std::ops::*;
use std::{f32, f64};
use vector::vec2::*;

#[derive(Copy, Clone, Debug)]
pub struct Vector3<T: Copy + Clone> {
	pub x: T,
	pub y: T,
    pub z: T,
}

#[allow(dead_code)]
impl<T: Copy + Clone> Vector3<T> {

	pub fn new(x: T, y: T, z: T) -> Vector3<T> {
		Vector3 {
			x: x,
			y: y,
            z: z
		}
	}

    pub fn from_vec2(v2: Vector2<T>, z: T) -> Vector3<T> {
        Vector3 {
            x: v2.x,
            y: v2.y,
            z: z
        }
    }

	pub fn dim(self) -> usize {
		3
	}

    /*
        NOTE: Swizzle functions below- very repetitive and prone to error.
        These were not generated.
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

    pub fn yx(self) -> Vector2<T> {
        self.xy().yx()
    }

    pub fn yz(self) -> Vector2<T> {
        Vector2 {
            x: self.y,
            y: self.z
        }
    }

    pub fn zx(self) -> Vector2<T> {
        self.xz().yx()
    }

    pub fn zy(self) -> Vector2<T> {
        self.yz().yx()
    }

    pub fn xzy(self) -> Vector3<T> {
        Vector3 {
            x: self.x,
            y: self.z,
            z: self.y
        }
    }

    pub fn yxz(self) -> Vector3<T> {
        Vector3 {
            x: self.y,
            y: self.x,
            z: self.z
        }
    }

    pub fn yzx(self) -> Vector3<T> {
        Vector3 {
            x: self.y,
            y: self.z,
            z: self.x
        }
    }

    pub fn zxy(self) -> Vector3<T> {
        Vector3 {
            x: self.z,
            y: self.x,
            z: self.y
        }
    }

    pub fn zyx(self) -> Vector3<T> {
        Vector3 {
            x: self.z,
            y: self.y,
            z: self.x
        }
    }
}

impl<T: Copy + Clone> Index<usize> for Vector3<T> {
	type Output = T;

	fn index(&self, idx: usize) -> &T {
		if idx >= 3 {
			panic!("index out of bounds: the len is 3 but the index is {}", idx);
		} else if idx == 0 {
			return &self.x;
		} else if idx == 1{
			return &self.y;
		} else {
            return &self.z;
        }
	}

}

impl<T: Copy + Clone> IndexMut<usize> for Vector3<T> {

	fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut T {
		if idx >= 3 {
			panic!("index out of bounds: the len is 3 but the index is {}", idx);
		} else if idx == 0 {
			return &mut self.x;
		} else if idx == 1 {
			return &mut self.y;
		} else {
            return &mut self.z;
        }
	}

}

impl<T: Copy + Clone + Add<Output = T>> Add for Vector3<T> {
	type Output = Vector3<T>;

	fn add(self, other: Self) -> Vector3<T> {
		Vector3 {
			x: self.x + other.x,
			y: self.y + other.y,
            z: self.z + other.z
		}
	}

}

impl<T: Copy + Clone + Add<Output = T>> AddAssign for Vector3<T> {

	fn add_assign(&mut self, other: Self) {
		self.x = self.x + other.x;
		self.y = self.y + other.y;
        self.z = self.z + other.z;
	}

}

impl<T: Copy + Clone + Sub<Output = T>> Sub for Vector3<T> {
	type Output = Vector3<T>;

	fn sub(self, other: Self) -> Vector3<T> {
		Vector3 {
			x: self.x - other.x,
			y: self.y - other.y,
            z: self.z - other.z
		}
	}

}

impl<T: Copy + Clone + Sub<Output = T>> SubAssign for Vector3<T> {

	fn sub_assign(&mut self, other: Self) {
		self.x = self.x - other.x;
		self.y = self.y - other.y;
        self.z = self.z - other.z;
	}

}

impl<T: Copy + Clone + Mul<Output = T>> Mul<T> for Vector3<T> {
	type Output = Vector3<T>;

	fn mul(self, rhs: T) -> Vector3<T> {
		Vector3 {
			x: self.x * rhs,
			y: self.y * rhs,
            z: self.z * rhs
		}
	}

}

impl<T: Copy + Clone + Mul<Output = T>> Mul<Vector3<T>> for Vector3<T> {
	type Output = Vector3<T>;

	fn mul(self, rhs: Vector3<T>) -> Vector3<T> {
		Vector3 {
			x: self.x * rhs.x,
			y: self.y * rhs.y,
            z: self.z * rhs.z
		}
	}

}

impl<T: Copy + Clone + Mul<Output = T>> MulAssign<T> for Vector3<T> {

	fn mul_assign(&mut self, rhs: T) {
		self.x = self.x * rhs;
		self.y = self.y * rhs;
        self.z = self.z * rhs;
	}

}

impl<T: Copy + Clone + Mul<Output = T>> MulAssign for Vector3<T> {

	fn mul_assign(&mut self, rhs: Vector3<T>) {
		self.x = self.x * rhs.x;
		self.y = self.y * rhs.y;
        self.z = self.z * rhs.z;
	}

}

impl<T: Copy + Clone + Div<Output = T>> Div<T> for Vector3<T> {
	type Output = Vector3<T>;

	fn div(self, rhs: T) -> Vector3<T> {
		Vector3 {
			x: self.x / rhs,
			y: self.y / rhs,
            z: self.z / rhs
		}
	}

}

impl<T: Copy + Clone + Div<Output = T>> DivAssign<T> for Vector3<T> {

	fn div_assign(&mut self, rhs: T) {
		self.x = self.x / rhs;
		self.y = self.y / rhs;
        self.z = self.z / rhs;
	}
}

impl<T: Copy + Clone + Neg<Output = T>> Neg for Vector3<T> {
	type Output = Vector3<T>;

	fn neg(self) -> Vector3<T> {
		Vector3 {
			x: -self.x,
			y: -self.y,
            z: -self.z
		}
	}

}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T>> Vector3<T> {

	pub fn dot(self,  r: Vector3<T>) -> T {
		let product: Vector3<T> = self * r;
		let dot: T = product.x + product.y + product.z;

		dot
	}

}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Sub<Output = T>> Vector3<T> {

    pub fn cross(self, r: Vector3<T>) -> Vector3<T> {
        let crss_x = (self.y * r.z) - (self.z * r.y);
        let crss_y = (self.z * r.x) - (self.x * r.z);
        let crss_z = (self.x * r.y) - (self.y * r.x);

        Vector3 {
            x: crss_x,
            y: crss_y,
            z: crss_z
        }
    }
}

#[allow(dead_code)]
impl<T: Copy + Clone + Into<f64> + From<f64>> Vector3<T> {

	pub fn mag(self) -> f64 {
		let x: f64 = self.x.into();
		let y: f64 = self.y.into();
        let z: f64 = self.z.into();

		((x * x) + (y * y) + (z * z)).sqrt()
	}

	pub fn proj(self, rhs: Vector3<T>) -> Vector3<T> {
		let l_converted: Vector3<f64> = Vector3::new(self.x.into(), self.y.into(), self.z.into());
		let r_converted: Vector3<f64> = Vector3::new(rhs.x.into(), rhs.y.into(), rhs.z.into());

		let r_normal = r_converted / r_converted.mag();
		let projected = r_normal * (l_converted.dot(r_converted));

		Vector3 {
			x: projected.x.into(),
			y: projected.y.into(),
            z: projected.z.into()
		}
	}
}

#[allow(dead_code)]
impl Vector3<f64> {

	pub fn norm(self) -> Vector3<f64> {

		let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z);
		let mag = mag_sqr.sqrt();

		let vec_scaled = self / mag;

		vec_scaled
	}

	pub fn normalize(&mut self) {
		let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z);
		let mag = mag_sqr.sqrt();

		self.x = self.x / mag;
		self.y = self.y / mag;
        self.z = self.z / mag;
	}

}

#[allow(dead_code)]
impl Vector3<f32> {

	pub fn norm(self) -> Vector3<f32> {

        let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z);
		let mag = mag_sqr.sqrt();

		let vec_scaled = self / mag;

		vec_scaled
	}

	pub fn normalize(&mut self) {
		let mag_sqr = (self.x * self.x) + (self.y * self.y) + (self.z * self.z);
		let mag = mag_sqr.sqrt();

		self.x = self.x / mag;
		self.y = self.y / mag;
        self.z = self.z / mag;
	}

}

pub type Vec3d = Vector3<f64>;
pub type Vec3 = Vector3<f32>;
pub type Vec3f = Vector3<f32>;

pub type Vec3ll = Vector3<i64>;
pub type Vec3i = Vector3<i32>;
pub type Vec3s = Vector3<i16>;
pub type Vec3b = Vector3<i8>;

pub type Vec3ull = Vector3<u64>;
pub type Vec3ui = Vector3<u32>;
pub type Vec3us = Vector3<u16>;
pub type Vec3ub = Vector3<u8>;
