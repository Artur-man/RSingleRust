use extendr_api::prelude::*;
use ndarray::ArrayView2; // ndarray for R matrices
use anndata_memory::{IMAnnData, IMArrayElement}; // in memory anndata for feature observation matrices
use anndata::ArrayData; // array data

pub mod io;

// /// Reads an h5ad file into memory
// /// @export
// #[extendr]
// fn read_h5ad_memory(file_path: String) {
//
//     // Load data and run complete analysis pipeline
//     let adata = io::read_h5ad_memory(&file_path);
//
// }

/// read dense matrix as an in-memory AnnData
/// @export
#[extendr]
fn read_matrix(matrix: ArrayView2<f64>, cells: Vec<String>, genes: Vec<String>) {

    // convert to anndata::ArrayData
    let array_data: ArrayData = matrix.to_owned().into();

    // Create the AnnData object
    let adata = IMAnnData::new_basic(
        array_data,
        cells,
        genes
    ).unwrap();

    // get the shape of x
    let adata_x: IMArrayElement = adata.x();
    let shape = adata_x.get_shape().unwrap();
    println!("Working with matrix of shape {:?}", shape);
}


/// temp function
/// @export
#[extendr]
fn read_matrix_hdf5(file_path: String) {

}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod RSingleRust;
    fn read_matrix;
    fn read_matrix_hdf5;
}
