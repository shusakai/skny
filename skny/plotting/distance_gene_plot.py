import plotly.graph_objects as go
import plotly
import pandas as pd
import seaborn as sns


def distance_gene_plot(
    grid, gene_symbol, 
    x_range=None, y_range=None, std=False, show=True, #min_grids=100
):
    '''
    grid: AnnData
    gene: str
    x_range, y_range: list
    std: bool
    min_grids: int
    '''
    ## get col and row length data from grid object --------------------------
    N_ROW = len(grid.uns["grid_yedges"]) - 1
    N_COL = len(grid.uns["grid_xedges"]) - 1

    # define variable which is name of table of gene expression by region--------------
    # if roi is selected
    if (x_range != None) | (y_range != None):
        name = "gene_expresion_by_region" \
        + "_x_["+", ".join([str(i) for i in x_range]) + "]" \
        + "_y_["+", ".join([str(i) for i in y_range]) + "]"
    else:
        name = "gene_expresion_by_region"

    # If it has already been calculated -----------------------------------------------
    if name in grid.uns.keys():
        df_grid_region_gene = grid.uns[name]
    else:
    
        ## Load shortest path and gene expression -------------------------------------
        df_shotest = grid.uns["shortest"]
        df_grid = grid.to_df()
        
        # complement nan grid
        df_grid = pd.merge(
            pd.DataFrame(index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)]), 
            df_grid, right_index=True, left_index=True, how="left"
        ).fillna(np.nan)
        
        df_region = pd.DataFrame(
            np.array(df_shotest["region"]).reshape(N_ROW, N_COL).T.reshape(N_ROW * N_COL), 
            index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)], columns=["region"]
        )
    
        ## ROI selection -------------------------------------------------
        if (x_range != None) | (y_range != None):
            df_roi = pd.DataFrame(
                np.array(df_shotest.index).reshape(N_ROW, N_COL).T.reshape(N_ROW * N_COL), 
                index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)], columns=["region"]
            )
            roi_grids = df_roi[df_roi["region"].isin(
                [i for i in df_roi.region if (x_range[0] <= i[0] <= x_range[1]) & (y_range[0] <= i[1] <= y_range[1])]
            )].index
            # extract ROI
            df_region = df_region.loc[roi_grids, :]
            df_grid = df_grid.loc[roi_grids, :]
    
        ## Merge shortest table and expression table --------------------------
        df_grid_region = pd.merge(
            df_grid, df_region, 
            right_index=True, left_index=True, how="left"
        )
        # Shaping
        df_grid_region = df_grid_region.dropna()
        df_grid_region_gene = df_grid_region.groupby("region").mean().T
        region_ls = df_region.region.value_counts()[df_region.region.value_counts() > 100].index # min grids
        df_grid_region_gene = df_grid_region_gene.loc[:,region_ls]
        df_grid_region_gene.columns = [i * 10 for i in df_grid_region_gene.columns]
        df_grid_region_gene = df_grid_region_gene.T.sort_index().T

        ## Save in uns----------------------------------------------------------------
        grid.uns[name] = df_grid_region_gene
    
    ## Histgram ------------------------------------------------------------
    # When standardizing for each gene
    if std:
        df_grid_region_gene = df_grid_region_gene.apply(lambda x: scipy.stats.zscore(x), axis=1)
        
    df_dist_gene = pd.DataFrame(df_grid_region_gene.loc[gene_symbol, :])
    df_dist_gene.index = [str(i) for i in df_dist_gene.index]
    fig = px.bar(
        df_dist_gene,
        y=gene_symbol,
    )

    # When standardizing for each gene
    if std:
        fig.update_layout(xaxis=dict(title='Region'),
                         yaxis=dict(title=f'{gene_symbol} (z-score)'))
        fig.update_layout(title=f'Z-score of {gene_symbol} expession of each region')
        
    else:
        fig.update_layout(xaxis=dict(title='Region'),
                         yaxis=dict(title=f'{gene_symbol}/μm<sup>2</sup>'))
        fig.update_layout(title=f'Density of {gene_symbol} expession of each region')
        
    fig.update_layout(plot_bgcolor="white")
    
    #plotly.io.write_html(fig, "QC_cxcl9_barplot.html")
    if show:
        fig.show()