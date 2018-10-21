# Peroxide

![travis](https://api.travis-ci.org/Axect/Peroxide.svg?branch=master)

Rust numeric library with R Syntax.

## Latest README version

Corresponds with `0.3.1`.

## Usage

### Initial Import

```rust
extern crate peroxide;
use peroxide::*;
```

### Vec\<f64\> Declaration

```R
# R
a = c(1,2,3,4)
b = seq(1,5,2) # (=c(1,3,5))
```

```rust
// Peroxide
let a = c!(1,2,3,4);
let b = seq!(1,5,2) // (=c!(1,3,5))
```

### Matrix Declaration

```R
# R
a = matrix(1:4, 2, 2, True)
```

```rust
// Peroxide (All belows are same)
// matrix function
let a = matrix(vec![1,2,3,4], 2, 2, Row);
let b = matrix(c!(1,2,3,4), 2, 2, Row);
let c = matrix(seq!(1,4,1), 2, 2, Row);

// matrix macro (More convenient)
let c = matrix!(1;4;1, 2, 2, Row);
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
let a = matrix!(1;4;1,  2, 2, Row);
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
print(a * b)
print(a %*% b)
```

```rust
// Peroxide
let a = matrix!(1;4;1, 2, 2, Row);
let b = matrix!(1;4;1, 2, 2, Col);
println!("{}", a.clone() + b.clone());
println!("{}", a.clone() - b.clone());
println!("{}", a.clone() * b.clone()); // Element-wise multiplication
println!("{}", a % b);  // Matrix multiplication
// Consume -> You can't use a,b anymore.
```

### Extract Column or Row

```R
a = matrix(1:4, 2, 2, T)
print(a[,1])
print(a[,2])
print(a[1,])
print(a[2,])
```

```rust
//Peroxide
let a = matrix!(1;4;1, 2, 2, Row);
println!("{}", a.col(0));
println!("{}", a.col(1));
println!("{}", a.row(0));
println!("{}", a.row(1));
```

### Functional Programming

```rust
// Peroxide
let a = matrix!(1;4;1, 2, 2, Row);
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

* ~~Extract row & col operator~~