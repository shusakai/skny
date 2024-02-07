import plotly.graph_objects as go
import plotly
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats
import numpy as np


def clustering_heatmap(
    grid, x_range=None, y_range=None, std=True, # min_grids=100
):
    '''
    grid : AnnData
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
        
        ## Merge shortest path and count of gene -------------------------------------
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
        
        df_grid_region = pd.merge(
            df_grid, df_region, 
            right_index=True, left_index=True, how="left"
        )
        
        df_grid_region = df_grid_region.dropna()
        df_grid_region_gene = df_grid_region.groupby("region").mean().T
        region_ls = df_grid_region.dropna().region.value_counts()[df_grid_region.dropna().region.value_counts() > 100].index
        df_grid_region_gene = df_grid_region_gene.loc[:,region_ls]
        df_grid_region_gene.columns = [i * 10 for i in df_grid_region_gene.columns]
        df_grid_region_gene = df_grid_region_gene.T.sort_index().T
    
        ## Save in uns----------------------------------------------------------------
        grid.uns[name] = df_grid_region_gene
    
    
    # Clustering ---------------------------------------------------------------
    # delete zero expression gene
    df_grid_region_gene = df_grid_region_gene[df_grid_region_gene.sum(axis=1) != 0]
    
    # if standardizing
    if std:
        df_grid_region_gene = df_grid_region_gene.apply(lambda x: scipy.stats.zscore(x), axis=1)
    
    g = sns.clustermap(
        df_grid_region_gene, 
        figsize=(8, 12), cmap="viridis", col_cluster=False, method="ward"
    )
    ax = g.ax_heatmap
    ax.set_ylabel("")
    plt.close() # not show
    # extract order of clustering
    row_order = g.dendrogram_row.reordered_ind
    df_grid_region_gene = df_grid_region_gene.iloc[row_order, :] # reorder by clustering

    ## Plotly code ----------------------------------------------------------------
    hover_text = []
    for e1, i in enumerate(df_grid_region_gene.columns):
        temp_ls = []
        for e2, gene in enumerate(df_grid_region_gene.index):
            temp_ls += [
                f"""{gene} {i}<br> 
    Mean density={(df_grid_region_gene.iloc[e2, e1] / (10 * 10)).round(3)}/μm<sup>2</sup><br> 
    Z-score={df_grid_region_gene.iloc[e2, e1].round(3)}"""
            ]
        hover_text += [temp_ls]
    hover_text = np.array(hover_text).T
    
    # restriction of horizontal zoom
    layout = go.Layout(
        xaxis=dict(
            fixedrange=True
        )
    )
    
    x_col = [int(float(i.split(", ")[0][1:])) for i in df_grid_region_gene.columns.astype(str)] +\
    [int(float(df_grid_region_gene.columns.astype(str)[-1].split(", ")[1][:-1]))]
    
    fig = go.Figure(
        data=go.Heatmap(
            z=df_grid_region_gene, 
            x=x_col,
            y=df_grid_region_gene.index,
            colorscale = 'Viridis', 
            hoverinfo="text",  # ホバー時にテキストを表示
            text=hover_text,  # ホバーテキストを指定, 
        ), layout=layout, 
    )
    fig.update_xaxes(
        tickvals=x_col,  # 表示する位置のインデックス
        ticktext=x_col  # ラベルのリスト
    )
    fig.update_layout(title='Clustering based on transcript density of each region')
    fig.update_layout(xaxis=dict(title='Region'), height=800,)
    fig.show()

    #plotly.io.write_html(fig, "grid_region_gene_std.html")