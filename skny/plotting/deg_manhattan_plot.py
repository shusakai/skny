from scipy import stats
from numpy import std, mean, sqrt
import plotly.graph_objects as go
import plotly
import pandas as pd
import seaborn as sns

#correct if the population S.D. is expected to be equal for the two groups.
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (mean(x) - mean(y)) / sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)

def calculate_q(p_seq):
  p_arr   = np.asarray(p_seq)
  N       = len(p_arr)
  idx_arr = np.argsort(p_arr)
  q_arr   = p_arr[idx_arr] * N / (np.arange(N) + 1)
  q_arr   = np.minimum.accumulate(q_arr[::-1])[::-1]
  q_arr[idx_arr] = q_arr.copy()
  return q_arr

def deg_manhattan_plot(
    grid, 
    x_range=None, y_range=None, show=True, #min_grids=100
):
    '''
    grid: AnnData
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
        name = "deg_by_region" \
        + "_x_["+", ".join([str(i) for i in x_range]) + "]" \
        + "_y_["+", ".join([str(i) for i in y_range]) + "]"
    else:
        name = "deg_by_region"

    # If it has already been calculated -----------------------------------------------
    if name in grid.uns.keys():
       df_q_res = grid.uns[name]
        
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
        # region selection
        df_grid_region = df_grid_region.dropna()
        region_ls = df_region.region.value_counts()[df_region.region.value_counts() > 100].index # min grids
        
        # control is all
        df_control = df_grid_region.drop("region", axis=1)
        df_control = df_control.applymap(lambda x: np.log(x+1))
        gene_ls = df_grid.columns.tolist()

        # place holder
        p_res_ls = []
        # target is expression which is belong to each region 
        for i in region_ls:
            df_target = df_grid_region.groupby("region").get_group(i).drop("region", axis=1)
            df_target = df_target.applymap(lambda x: np.log(x+1))
            
            temp_p_res_ls = []        
            for s in gene_ls:
                A = df_target[s]
                B = df_control[s]
                # one sided test
                temp_p_res_ls += [stats.ttest_ind(A, B, equal_var = False, alternative="greater").pvalue]
        
            p_res_ls += [temp_p_res_ls]

        # dataframe
        df_p_res = pd.DataFrame(
            p_res_ls, 
            index=region_ls, 
            columns=gene_ls
        ).T

        # adjust
        df_q_res = pd.DataFrame(
            calculate_q(df_p_res.values.reshape(len(region_ls) * df_grid.shape[1])).reshape(df_grid.shape[1], len(region_ls)).T,
            index=region_ls, 
            columns=gene_ls
        ).T.fillna(1)

        
        df_q_res = df_q_res.apply(lambda x: np.log10(x) * -1)
        df_q_res = df_q_res.T.sort_index().T

        ## Save in uns----------------------------------------------------------------
        grid.uns[name] = df_q_res

    # Plot -----------------------------------------------------------------
    df_q_res_stack = pd.DataFrame(
        df_q_res.stack()
    ).reset_index()
    df_q_res_stack.columns = ["Gene", "Region", "-log(qvalue)"]
    df_q_res_stack.Region = df_q_res_stack.Region.astype(str)
    
    # ジッタープロットを作成
    fig = px.strip(df_q_res_stack, x='Region', y="-log(qvalue)", color="Gene", 
                  )
    
    # edit layouts
    fig.update_layout(
        title='Differential abundance gene of each region',
        yaxis_title='-log(qvalue)',
    )
    
    for i in range(-1, len(df_q_res_stack.Region.unique()) -1 ):
        fig.add_vrect(
            x0=0.5+i, x1=0.5+i, row="all", col=1, line_width=0.2
        )
    
    
    fig.update_xaxes(matches=None)
    fig.update_layout(plot_bgcolor="white")
    fig.update_layout(xaxis=dict(title=None))
    fig.for_each_xaxis(lambda xax: xax.update(title=None))

    fig.show()