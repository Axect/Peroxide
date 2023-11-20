# Ver 0.1.3 (2023-11-09)

- Add `Group` and `Ring`
  - `Group`: require `Add<Self, Output=Self>` and has `zero()`
  - `Ring`: require `Group + Mul<Self, Output=Self>` and has `one()`
  - `f32` and `f64` are `Ring`

# Ver 0.1.2 (2023-11-09)

- Remove unnecessary constraints for `Numeric<T>`
  - `T: Add<Self, Output=Self> + Sub<Self, Output=Self> + Mul<Self, Output=Self> + Div<Self, Output=Self>`
- More specify `Neg` for `Numeric` & `Float`: `Neg` -> `Neg<Output=Self>` 

# Ver 0.1.1 (2023-11-09)

- Add constraints for `Numeric<T>`
  - `T: Add<Self, Output=Self> + Sub<Self, Output=Self> + Mul<Self, Output=Self> + Div<Self, Output=Self>`
- Add example : `Vec3D`
