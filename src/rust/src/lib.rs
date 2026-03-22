use extendr_api::prelude::*;
use ndarray::ArrayView2; // ndarray for R matrices

use anndata::{
    ArrayData,
    Readable,
    backend::{Backend, DataContainer, GroupOp},
    data::DynCsrMatrix,
    data::DynCscMatrix
};
use anndata_memory::{IMAnnData, IMArrayElement}; // in memory anndata for feature observation matrices
use anndata_hdf5::H5;

use nalgebra_sparse::CsrMatrix;
use nalgebra_sparse::CscMatrix;

use anyhow::Result;

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
fn read_matrix_hdf5(file_path: &str, group_path: &str) -> Result<()> {

    // read from h5
    let file = H5::open(file_path)?;
    let group = file.open_group(group_path)?;
    let container = DataContainer::<H5>::Group(group);
    let csr = CsrMatrix::<f64>::read(&container)?;
    let matrix = DynCsrMatrix::from(csr);

    // convert to anndata::ArrayData
    let array_data = ArrayData::CsrMatrix(matrix);

    Ok(())
}

/// temp function
/// @export
#[extendr]
fn read_matrix_SCE(obj: S4, slot: String, cells: Vec<String>, genes: Vec<String>) {

    // read singlecellexperiment data
    let assays: S4 = obj.get_slot("assays").unwrap().try_into().unwrap();
    let data_slot: S4 = assays.get_slot(slot).unwrap().try_into().unwrap();
    let list_data: List = data_slot.get_slot("listData").unwrap().try_into().unwrap();

    // visualize type
    println!("{}", std::any::type_name_of_val(&assays));

    let counts: ArrayView2<i32> = list_data.dollar("counts").unwrap().try_into().unwrap();

    // visualize type
    println!("{}", std::any::type_name_of_val(&counts));

    // get shape
    // let shape = counts.shape();
    // println!("Working with matrix of shape {:?}", shape);

    // convert to anndata::ArrayData
    let array_data: ArrayData = counts.to_owned().into();

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

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod RSingleRust;
    fn read_matrix;
    fn read_matrix_hdf5;
    fn read_matrix_SCE;
}
