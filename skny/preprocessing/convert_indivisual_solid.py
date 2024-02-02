import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import cv2


def convert_indivisual_solid(grid):
    '''
    grid: AnnData after `distance_calculator`
    '''
    
    ## get col and row length data from grid object --------------------------
    N_ROW = len(grid.uns["grid_yedges"]) - 1
    N_COL = len(grid.uns["grid_xedges"]) - 1
    
    # show separation of indivisual tumor solid using rectangle -------------------------------------------
    # load contour image from skny object
    img_color = grid.uns["marker_delineation"].copy()
    img_gray = cv2.cvtColor(img_color, cv2.COLOR_BGR2GRAY)
    
    # extraction of indivisual tumor-solid from contoured image
    retval, labels, stats, centroids = cv2.connectedComponentsWithStats(img_gray)
    
    # delineation of indivisual tumor-solid using rectangle
    for i in range(1, retval):
        x, y, width, height, area = stats[i] # extract regtangle information
        if area > 10: # exclude area<100μm^2
            cv2.rectangle(img_color, (x - 3, y - 3), (x+width+3, y+height+3), (0, 0, 255), thickness=1)
            cv2.putText(img_color, f"{i}", (x + 3, y + 3), cv2.FONT_HERSHEY_PLAIN, 0.8, (0, 0, 255), 1, cv2.LINE_8)
    
    # Load region annotation dataframe
    df_shotest = grid.uns["shortest"]
    df_shotest["right"] = df_shotest["region"].apply(lambda x: x.right)
    
    # convert to dataframe
    df_solid = pd.DataFrame(
        labels.reshape(N_COL, N_ROW).reshape(N_ROW * N_COL), 
        index=df_shotest.index, columns=["solid"]
    )
    
    # merge to df_shortest
    df_shotest = pd.merge(
        df_shotest, df_solid, right_index=True, left_index=True, how="left"
    )
    
    # extract (−∞, 0] section from indivisual tumor solid ----------------------------------------
    solid_ls = []
    for i, s in zip(df_shotest["right"], df_shotest["solid"]):
        if (i > 0 ) & (s != 0 ): # select solid regeion and (−∞, 0] section 
            solid_ls += [0]
        else:
            solid_ls += [s]
    # add to dataframe
    df_shotest["solid"] = solid_ls
    
    
    # merge region annotations and gene expression data -----------------------------------
        # load gene expression of each grid 
    df_grid = grid.to_df()
    
    # complement nan grid
    df_grid = pd.merge(
        pd.DataFrame(index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)]), 
        df_grid, right_index=True, left_index=True, how="left"
    ).fillna(np.nan)
    # complement nan grid
    df_solid = pd.DataFrame(
        np.array(df_shotest["solid"]).reshape(N_ROW, N_COL).T.reshape(N_ROW * N_COL), 
        index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)], columns=["solid"]
    )
    # merge
    df_grid_solid = pd.merge(
        df_grid, df_solid, 
        right_index=True, left_index=True, how="left"
    )
    
    # mean based on (−∞, 0] section of indivisual tumor solid
    df_solid = df_grid_solid.dropna().groupby("solid").mean()
    # exclude untargeted coodinates
    df_solid = df_solid[df_solid.index != 0]
    
    # generate annotation matrix of each tumor solid 
    df_solid_centroid = pd.DataFrame(
        centroids[df_solid.index], 
        index=df_solid.index,
        columns=["imagecol", "imagerow"]
    )
    df_solid_centroid["solid"] = df_solid_centroid.index
    
    
    # create AnnData object of (−∞, 0] section of indivisual tumor solid ----------------------------
    counts = df_solid.values
    solid = ad.AnnData(counts, uns={"spatial": {"data": {
        'scalefactors': {'tissue_hires_scalef': 1, 'spot_diameter_fullres': 10}, 
        'use_quality': 'hires',
        'images': {'hires':np.array([[255, 255, 255, 255]])}
    }}})
    
    # add cell and gene name
    solid.obs_names = df_solid.index.tolist()
    solid.var_names = df_solid.columns.tolist()
    
    # add centroid to metadata
    solid.obs = df_solid_centroid
    
    # annotation processing image to AnnData object ------------------------------------
    solid.uns["indivisual_tumor_solid"] = img_color
    # annotation of dataframe
    solid.uns["shortest"] = df_shotest

    return solid