import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import cv2

def convert_indivisual_peri(grid):
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

            
    # take all-region-XOR based on all rectangle ---------------------------------
    xor_result = None # place holder
    for i in range(1, retval): # each rectangle
        x, y, width, height, area = stats[i] # extract regtangle information
        if area > 10: # exclude area<100μm^2
            image = np.zeros((N_ROW, N_COL), dtype=np.uint8) * 255 # generate zero matrix of same size
            cv2.rectangle(image, (x - 3, y - 3), (x+width+3, y+height+3), 255, thickness=-1) # delineate rectangle
            
            # rectangle is added to zero matrix
            if xor_result is None:
                xor_result = image
            else:
                xor_result = xor_result + image
            

    # select coordinates where rectangle is unique
    xor_result = np.where(xor_result != 255, 0, 255)
    # convert datatype
    xor_result = np.array(xor_result, dtype="uint8")

    # to attribute XOR region to indivisual tumor solid we take AND between all-region-XOR image and indivisual rectangle----------------
    solid_peri = None # place holder
    for i in range(1, retval): # each rectangle
        x, y, width, height, area = stats[i] # extract regtangle information
        if area > 10: # exclude area<100μm^2
            image = np.zeros((N_ROW, N_COL), dtype=np.uint8) * 255 # generate zero matrix of same size
            cv2.rectangle(image, (x - 3, y - 3), (x+width+3, y+height+3), 255, thickness=-1) # delineate rectangle
            and_result = cv2.bitwise_and(xor_result, image) # AND

            # white value (255) is converted to rectangle label, and these labels on coodinate added to zero matrix
            if solid_peri is None:
                solid_peri = (and_result / 255) * i
            else:
                solid_peri = solid_peri + ((and_result / 255) * i)

    # convert to matrix to vector
    solid_peri = solid_peri.reshape(N_COL, N_ROW).reshape(N_ROW * N_COL)

    # dataframe
    df_temp_solid = pd.DataFrame(
        solid_peri, 
        index=df_shotest.index, columns=["solid_peri"]
    )

    # Load region annotation dataframe
    df_shotest = grid.uns["shortest"]
    df_shotest["right"] = df_shotest["region"].apply(lambda x: x.right)
    # merge to dataframe
    df_shotest = pd.merge(
        df_shotest, df_temp_solid, right_index=True, left_index=True, how="left"
    )


    # extract (0, +30] section from XOR region of indivisual tumor solid ----------------------------------------
    solid_peri_ls = [] # place holder
    for i, s in zip(df_shotest["right"], df_shotest["solid_peri"]):
        if (i != 3 ) & (s != 0 ): # select XOR regeion and (0, +30] section 
            solid_peri_ls += [0]
        else:
            solid_peri_ls += [s]
    # add to dataframe
    df_shotest["solid_peri"] = solid_peri_ls


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
        np.array(df_shotest["solid_peri"]).reshape(N_ROW, N_COL).T.reshape(N_ROW * N_COL), 
        index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)], columns=["solid_peri"]
    )
    # merge
    df_grid_solid = pd.merge(
        df_grid, df_solid, 
        right_index=True, left_index=True, how="left"
    )

    # mean based on (0, +30] section of indivisual tumor solid
    df_solid = df_grid_solid.dropna().groupby("solid_peri").mean()
    # exclude untargeted coodinates
    df_solid = df_solid[df_solid.index != 0]

    # generate annotation matrix of each tumor solid 
    df_solid_centroid = pd.DataFrame(
        centroids[[int(i) for i in df_solid.index]], 
        index=[int(i) for i in df_solid.index],
        columns=["imagecol", "imagerow"])
    df_solid_centroid["solid_peri"] = df_solid_centroid.index

    
    # create AnnData object of (0, +30] section of indivisual tumor solid ----------------------------
    counts = df_solid.values
    solid_peri = ad.AnnData(counts, uns={"spatial": {"data": {
        'scalefactors': {'tissue_hires_scalef': 1, 'spot_diameter_fullres': 10}, 
        'use_quality': 'hires',
        'images': {'hires':np.array([[255, 255, 255, 255]])}
    }}})

    # add cell and gene name
    solid_peri.obs_names = df_solid.index.tolist()
    solid_peri.var_names = df_solid.columns.tolist()

    # add centroid to metadata
    solid_peri.obs = df_solid_centroid

    # annotation processing image to AnnData object ------------------------------------
    solid_peri.uns["indivisual_tumor_solid"] = img_color
    solid_peri.uns["xor_result"] = xor_result
    # annotation of dataframe
    solid_peri.uns["shortest"] = df_shotest

    return solid_peri
