#[macro_export]
macro_rules! hstack {
    // Two
    ( $x:expr, $y:expr ) => {
        {
            let l = $x.len();
            let mut temp = $x;
            temp.extend(&$y);
            matrix(temp, l, 2, Col)
        }
    };

    // Multi
    ( $x0:expr, $( $x: expr ),* ) => {
        {
            let r = $x0.len();
            let mut temp0 = $x0;
            let mut c = 1usize;

            $(
                let mut temp = $x;
                // Must equal row
                assert_eq!(r, temp.len());
                // Add column
                c += 1;
                temp0.extend_from_slice(&temp[..]);
            )*
            matrix(temp0, r, c, Col)
        }
    };
}

#[macro_export]
macro_rules! vstack {
    // Two
    ( $x:expr, $y:expr ) => {
        {
            let l = $x.len();
            let mut temp = $x;
            temp.extend(&$y);
            matrix(temp, 2, l, Row)
        }
    };

    // Multi
    ( $x0:expr, $( $x: expr ),* ) => {
        {
            let c = $x0.len();
            let mut temp0 = $x0;
            let mut r = 1usize;

            $(
                let mut temp = $x;
                // Must equal row
                assert_eq!(c, temp.len());
                // Add column
                r += 1;
                temp0.extend_from_slice(&temp[..]);
            )*
            matrix(temp0, r, c, Row)
        }
    };
}
