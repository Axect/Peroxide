# Peroxide

<!-- ![travis](https://api.travis-ci.org/Axect/Peroxide.svg?branch=master) -->

Rust numeric library with R Syntax.

## Usage

### Matrix Declaration

```R
# R
a = matrix(1:4, 2, 2, True)
```

```rust
// Peroxide
let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
```

### Print

```R
# R
a = matrix(1:4, 2, 2, True)
print(a)
#      [,1] [,2]
# [1,]    1    2
# [2,]    3    4
```

```rust
// Peroxide
let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
println!("{}", a);
//       c[0] c[1]
// r[0]     1    2
// r[1]     3    4
```

### Matrix operation

* If you want to do multiple operations on same matrix, then you should use `clone` because Rust `std::ops` consume value. 

```R
# R
a = matrix(1:4, 2, 2, T)
b = matrix(1:4, 2, 2, F)
print(a + b)
print(a - b)
print(a %*% b)
```

```rust
// Peroxide
let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
let b = Matrix::new(vec![1,2,3,4], 2, 2, Col);
println!("{}", a.clone() + b.clone());
println!("{}", a.clone() - b.clone());
println!("{}", a * b); // Consume -> You can't use a,b anymore.
```

### Functional Programming

```rust
// Peroxide
let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
println!("{}", a.fmap(|x| x + 1.0));
println!("{}", a.fmap(|x| x - 1.0));
println!("{}", a.fmap(|x| x * 2.0));

// Results
//
//       c[0] c[1]
// r[0]     2    3
// r[1]     4    5
//
//       c[0] c[1]
// r[0]     0    1
// r[1]     2    3
//
//       c[0] c[1]
// r[0]     2    4
// r[1]     6    8
```

## TODO

* Extract row & col operator