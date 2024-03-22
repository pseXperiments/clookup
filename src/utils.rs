pub fn log(n: usize) -> u32 {
    let mut k = 0;
    let mut m = n;
    while m > 1 {
        m >>= 1;
        k += 1;
    }
    if n & (n - 1) == 0 {
        k
    } else {
        k + 1
    }
}