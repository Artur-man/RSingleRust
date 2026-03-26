use anndata_memory::IMAnnData;

use single_rust::memory::statistics::calculate_qc_metrics;
use polars::prelude::Column;

pub fn custom_qc_metrics(adata: &IMAnnData) -> anyhow::Result<()> {
    let var_names = adata.var_names();
    let mito_mask: Vec<bool> = var_names
        .iter()
        .map(|name| name.starts_with("MT-") || name.starts_with("mt-"))
        .collect();

    let mut var_df = adata.var().get_data();
    var_df.with_column(Column::new("mito".into(), mito_mask))?;
    adata.var().set_data(var_df)?;

    calculate_qc_metrics(
        adata,
        Some("counts"),
        Some("genes"),
        Some(vec!["mito"]),
        Some(vec![]),
        None,
        false,
        true,
        true,
    )?;

    Ok(())
}
