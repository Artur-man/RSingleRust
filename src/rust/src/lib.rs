use extendr_api::prelude::*;
use ndarray::ArrayView2; // ndarray for R matrices
// use ndarray::Array2;

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
// use arrow_extendr::IntoArrowRobj;

use nalgebra_sparse::CsrMatrix;
use nalgebra_sparse::CscMatrix;

use anyhow::Result;

pub mod io;
pub mod statistics;

/// read dense matrix as an in-memory AnnData
/// @export
#[extendr]
fn read_matrix_dense(matrix: ArrayView2<f64>, cells: Vec<String>, genes: Vec<String>) {

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

fn read_matrix_SCE(obj: S4, cells: Vec<String>, genes: Vec<String>) -> IMAnnData {

    // read singlecellexperiment data
    // let assays: S4 = obj.get_slot("assays").unwrap().try_into().unwrap();
    // let data_slot: S4 = assays.get_slot(slot).unwrap().try_into().unwrap();
    // let list_data: List = data_slot.get_slot("listData").unwrap().try_into().unwrap();

    // visualize type
    // let counts: S4 = list_data.dollar("counts").unwrap().try_into().unwrap();
    // let csr = convert_to_CscMatrix(counts).unwrap();
    let csr = convert_to_CscMatrix(obj).unwrap();
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

/// temp function
/// @export
#[extendr]
fn get_qc_metrics(obj: S4, slot: String, cells: Vec<String>, genes: Vec<String>) -> anyhow::Result<()> {

  let adata: IMAnnData = read_matrix_SCE(obj, cells, genes);

  // qc, qc_metrics gives error because defaults top values are too low
  statistics::custom_qc_metrics(&adata)?;

  // get obs matrix
  let obs_df: DataFrame = adata.obs().get_data();
  println!("Working with column names {:?}", obs_df.get_column_names());
  //let df = obs_df.into_arrow_robj().map_err(|e| anyhow!("{e:?}"));

  // let n_genes_by_counts = obs_df.column("n_genes_by_counts").unwrap();
  // let n_genes_by_counts = obs_df.column("n_genes_by_counts").unwrap().i32().unwrap().into_iter()
  //                   .map(|v| v.unwrap_or(i32::MIN))
  //                  .collect();

  // let ndarray = obs_df.to_ndarray::<Float64Type>(IndexOrder::Fortran).unwrap();
  // println!("{}", std::any::type_name_of_val(&n_genes_by_counts));

  //Ok(data_frame!(n_genes_by_counts = n_genes_by_counts, total_counts = total_counts))
  // Ok(df.unwrap())
  Ok(())
}

/// temp function
/// @export
#[extendr]
fn update_S4(obj: S4, slot: String) -> S4 {

    // read singlecellexperiment data
    let assays: S4 = obj.get_slot("assays").unwrap().try_into().unwrap();
    let data_slot: S4 = assays.get_slot(slot).unwrap().try_into().unwrap();
    let list_data: List = data_slot.get_slot("listData").unwrap().try_into().unwrap();

    // visualize type
    println!("{}", std::any::type_name_of_val(&assays));

    // get counts
    // let counts: ArrayView2<i32> = list_data.dollar("counts").unwrap().try_into().unwrap();

    // update counts:
    // Make a new matrix with all values doubled
    //let doubled: Array2<i32> = &counts * 2;
    //let new_counts: Robj = doubled.try_into().unwrap();
    // find the index of "counts"
    //let idx = list_data
    //  .iter()
    //  .position(|(name, _)| name == "counts")
    //  .expect("counts not found");
    //list_data.set_elt(idx, new_counts).unwrap();

    return assays;
}

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

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod RSingleRust;
    fn read_matrix_dense;
    fn read_matrix_hdf5;
    fn update_S4;
    fn get_qc_metrics;
}
