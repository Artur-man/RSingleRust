use extendr_api::prelude::*;

pub mod io;

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world2!"
}

/// @export
#[extendr]
fn read_h5ad_memory(file_path: String) {

     // Load data and run complete analysis pipeline
     let adata = io::read_h5ad_memory(&file_path);

}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rusttest;
    fn hello_world;
    fn read_h5ad_memory;
}
