use std::ops::*;
use std::fmt::*;
use std::{f32, f64};

#[derive(Copy, Clone, Debug)]
pub struct Vector2<T: Copy + Clone> {
	pub x: T,
	pub y: T,
}

#[allow(dead_code)]
impl<T: Copy + Clone> Vector2<T> {

	pub fn new(x: T, y: T) -> Vector2<T> {
		Vector2 {
			x: x,
			y: y
		}
	}

}

impl<T: Copy + Clone + Debug> Index<usize> for Vector2<T> {
	type Output = T;

	fn index(&self, idx: usize) -> &T {
		if idx >= 2 {
			panic!("Index {} out of bounds for {:?}", idx, *self);
		} else if idx == 0 {
			return &self.x;
		} else {
			return &self.y;
		}
	}

}

impl<T: Copy + Clone + Debug> IndexMut<usize> for Vector2<T> {

	fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut T {
		if idx >= 2 {
			panic!("Index {} out of bounds for {:?}", idx, *self);
		} else if idx == 0 {
			return &mut self.x;
		} else {
			return &mut self.y;
		}
	}

}

impl<T: Copy + Clone + Add<Output = T>> Add for Vector2<T> {
	type Output = Vector2<T>;

	fn add(self, other: Self) -> Vector2<T> {
		Vector2 {
			x: self.x + other.x,
			y: self.y + other.y
		}
	}

}

impl<T: Copy + Clone + Add<Output = T>> AddAssign for Vector2<T> {

	fn add_assign(&mut self, other: Self) {
		self.x = self.x + other.x;
		self.y = self.y + other.y;
	}

}

impl<T: Copy + Clone + Sub<Output = T>> Sub for Vector2<T> {
	type Output = Vector2<T>;

	fn sub(self, other: Self) -> Vector2<T> {
		Vector2 {
			x: self.x - other.x,
			y: self.y - other.y
		}
	}

}

impl<T: Copy + Clone + Sub<Output = T>> SubAssign for Vector2<T> {

	fn sub_assign(&mut self, other: Self) {
		self.x = self.x - other.x;
		self.y = self.y - other.y;
	}

}

impl<T: Copy + Clone + Mul<Output = T>> Mul<T> for Vector2<T> {
	type Output = Vector2<T>;

	fn mul(self, rhs: T) -> Vector2<T> {
		Vector2 {
			x: self.x * rhs,
			y: self.y * rhs
		}
	}

}

impl<T: Copy + Clone + Mul<Output = T>> Mul for Vector2<T> {
	type Output = Vector2<T>;

	fn mul(self, rhs: Vector2<T>) -> Vector2<T> {
		Vector2 {
			x: self.x * rhs.x,
			y: self.y * rhs.y
		}
	}

}

impl<T: Copy + Clone + Mul<Output = T>> MulAssign<T> for Vector2<T> {

	fn mul_assign(&mut self, rhs: T) {
		self.x = self.x * rhs;
		self.y = self.y * rhs;
	}

}

impl<T: Copy + Clone + Mul<Output = T>> MulAssign for Vector2<T> {

	fn mul_assign(&mut self, rhs: Vector2<T>) {
		self.x = self.x * rhs.x;
		self.y = self.y * rhs.y;
	}

}

impl<T: Copy + Clone + Div<Output = T>> Div<T> for Vector2<T> {
	type Output = Vector2<T>;

	fn div(self, rhs: T) -> Vector2<T> {
		Vector2 {
			x: self.x / rhs,
			y: self.y / rhs,
		}
	}

}

impl<T: Copy + Clone + Div<Output = T>> DivAssign<T> for Vector2<T> {

	fn div_assign(&mut self, rhs: T) {
		self.x = self.x / rhs;
		self.y = self.y / rhs;
	}
}

impl<T: Copy + Clone + Neg<Output = T>> Neg for Vector2<T> {
	type Output = Vector2<T>;

	fn neg(self) -> Vector2<T> {
		Vector2 {
			x: -self.x,
			y: -self.y
		}
	}

}

#[allow(dead_code)]
impl<T: Copy + Clone + Mul<Output = T> + Add<Output = T>> Vector2<T> {

	pub fn dot(l: Vector2<T>,  r: Vector2<T>) -> T {
		let product: Vector2<T> = l * r;
		let dot: T = product.x + product.y;

		dot
	}

}

#[allow(dead_code)]
impl Vector2<f64> {

	pub fn norm(self) -> Vector2<f64> {

		let mag_sqr = (self.x * self.x) + (self.y * self.y);
		let mag = mag_sqr.sqrt();

		let vec_scaled = self / mag;

		vec_scaled
	}

	pub fn normalize(&mut self) {
		let mag_sqr = (self.x * self.x) + (self.y * self.y);
		let mag = mag_sqr.sqrt();

		self.x = self.x / mag;
		self.y = self.y / mag;
	}

	pub fn angle_rad(self) -> f64 {
		self.y.atan2(self.x)
	}

	pub fn angle_deg(self) -> f64 {
		self.angle_rad().to_degrees()
	}

}

#[allow(dead_code)]
impl Vector2<f32> {
	pub fn norm(self) -> Vector2<f32> {

		let mag_sqr = (self.x * self.x) + (self.y * self.y);
		let mag = mag_sqr.sqrt();

		let vec_scaled = self / mag;

		vec_scaled
	}

	pub fn normalize(&mut self) {
		let mag_sqr = (self.x * self.x) + (self.y * self.y);
		let mag = mag_sqr.sqrt();

		self.x = self.x / mag;
		self.y = self.y / mag;
	}

	pub fn angle_rad(self) -> f32 {
		self.y.atan2(self.x)
	}

	pub fn angle_deg(self) -> f32 {
		self.angle_rad().to_degrees()
	}

}

pub type Vec2d = Vector2<f64>;
pub type Vec2 = Vector2<f32>;

pub type Vec2l = Vector2<i64>;
pub type Vec2i = Vector2<i32>;
pub type Vec2s = Vector2<i16>;
pub type Vec2b = Vector2<i8>;

pub type Vec2ul = Vector2<u64>;
pub type Vec2ui = Vector2<u32>;
pub type Vec2us = Vector2<u16>;
pub type Vec2ub = Vector2<u8>;
