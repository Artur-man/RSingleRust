use extendr_api::prelude::*;
use ndarray::ArrayView2; // ndarray for R matrices
use anyhow::{anyhow, Result};

// anndata
use anndata::{
    ArrayData,
    Readable,
    backend::{Backend, DataContainer, GroupOp},
    data::DynCsrMatrix,
    data::DynCscMatrix
};
use anndata_memory::{IMAnnData, IMArrayElement}; // in memory anndata for feature observation matrices
use anndata_hdf5::H5;
use polars_core::frame::DataFrame;

// arrow extendr for polars integration
// use arrow_extendr::IntoArrowRobj;

// sparse matrix libraries
use nalgebra_sparse::CsrMatrix;
use nalgebra_sparse::CscMatrix;

fn convert_to_CscMatrix(obj: S4) -> Result<CscMatrix<f64>> {
    let row_indices_r: Vec<i32> = obj.get_slot("i").unwrap().try_into().unwrap();
    let col_offsets_r: Vec<i32> = obj.get_slot("p").unwrap().try_into().unwrap();
    let dim: Vec<i32> = obj.get_slot("Dim").unwrap().try_into().unwrap();
    let values: Vec<f64> = obj.get_slot("x").unwrap().try_into().unwrap();

    let nrows = dim[0] as usize;
    let ncols = dim[1] as usize;

    let row_indices: Vec<usize> = row_indices_r.into_iter().map(|v| v as usize).collect();
    let col_offsets: Vec<usize> = col_offsets_r.into_iter().map(|v| v as usize).collect();

    let csc = CscMatrix::try_from_csc_data(nrows, ncols, col_offsets, row_indices, values).unwrap();

    rprintln!("constructed CSC: {} x {}, nnz = {}", csc.nrows(), csc.ncols(), csc.nnz());

    Ok(csc)
}

pub fn read_matrix_SCE(obj: S4, cells: Vec<String>, genes: Vec<String>) -> IMAnnData {

    // convert to cscmatrix first
    let csr = convert_to_CscMatrix(obj).unwrap();

    // define as IMAnnData::ArrayData
    let matrix = DynCscMatrix::from(csr);
    let array_data = ArrayData::CscMatrix(matrix);

    // Create the AnnData object
    let adata = IMAnnData::new_basic(
        array_data,
        cells,
        genes
    ).unwrap();

    // get the shape of x
    let adata_x: IMArrayElement = adata.x();
    let shape = adata_x.get_shape().unwrap();
    println!("IMAnnData object with shape {:?}", shape);

    return adata;
}
