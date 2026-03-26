use extendr_api::prelude::*;
use ndarray::ArrayView2; // ndarray for R matrices
use anyhow::{anyhow, Result};

fn print_2x2(matrix: ArrayView2<f64>) {
    let nrows = matrix.nrows().min(2);
    let ncols = matrix.ncols().min(2);

    rprintln!("Top-left {}x{} block:", nrows, ncols);

    for i in 0..nrows {
        let mut row = Vec::with_capacity(ncols);
        for j in 0..ncols {
            row.push(matrix[[i, j]]);
        }
        rprintln!("{:?}", row);
    }
}
