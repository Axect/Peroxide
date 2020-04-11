extern crate rand;
use rand::prelude::*;

/// Extract no-duplicate sample from Vector
pub trait SampleRNG {
    type Item;
    fn sample(&self, n: usize) -> Vec<Self::Item>;
}

impl SampleRNG for Vec<usize> {
    type Item = usize;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<u32> {
    type Item = u32;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<u64> {
    type Item = u64;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<isize> {
    type Item = isize;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<i32> {
    type Item = i32;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<i64> {
    type Item = i64;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<f64> {
    type Item = f64;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<f32> {
    type Item = f32;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for Vec<char> {
    type Item = char;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl<'a> SampleRNG for Vec<&'a str> {
    type Item = &'a str;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.iter().map(|x| *x).choose_multiple(&mut rng, n)
    }
}

impl SampleRNG for String {
    type Item = char;
    fn sample(&self, n: usize) -> Vec<Self::Item> {
        let mut rng = thread_rng();
        self.chars().choose_multiple(&mut rng, n)
    }
}