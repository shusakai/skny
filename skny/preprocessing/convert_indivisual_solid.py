import numpy as np
import pandas as pd
import anndata as ad
import cv2


def convert_indivisual_solid(grid, section=None):
    """Pseudo-bulk gene expression of each individual tumor nest ("single-TME" analysis).

    Connected components of the tumor contour image are treated as individual tumor nests
    (solids). For every nest the mean expression of all genes over its bins is returned as
    one observation, so the result can be analysed like a bulk / pseudo-bulk dataset
    (clustering of nests, trajectory, correlation with the peri-tumoral band from
    :func:`skny.pp.convert_indivisual_peri`).

    Parameters
    ----------
    grid : AnnData
        Output of :func:`skny.pp.calculate_distance`.
    section : int or None
        Distance band (um, inside the tumor) to aggregate: one of ``0, -30, -60, -90, -120``
        (e.g. ``-30`` = the band 30 to 0 um inside the contour). ``None`` aggregates the whole
        interior of each nest, i.e. the interval (-inf, 0].

    Returns
    -------
    AnnData
        One observation per tumor nest (``obs_names`` ``"solid_<id>"``), ``var`` = genes,
        ``obs["imagecol"]`` / ``obs["imagerow"]`` = nest centroid, ``obs["solid"]`` = nest id,
        ``uns["indivisual_tumor_solid"]`` = image with the nests boxed and numbered.
        Nests smaller than 10 bins are ignored.

    Examples
    --------
    >>> solid = sk.pp.convert_indivisual_solid(grid)             # whole interior
    >>> solid_30 = sk.pp.convert_indivisual_solid(grid, section=-30)
    """
    
    ## get col and row length data from grid object --------------------------
    N_ROW = len(grid.uns["grid_yedges"]) - 1
    N_COL = len(grid.uns["grid_xedges"]) - 1

    if section in [0, -30, -60, -90, -120]:
        section = int(section / 10)

    elif section is None:
        print("Calculatting gene expression of (-∞, 0]")
    
    else:
        return "Section must be 0, -30, -60, -90, or -120"
    
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
            cv2.rectangle(img_color, (x, y), (x+width, y+height), (0, 0, 255), thickness=1)
            cv2.putText(img_color, f"{i}", (x, y), cv2.FONT_HERSHEY_PLAIN, 0.8, (0, 0, 255), 1, cv2.LINE_8)
    
    # Load region annotation dataframe
    #df_shotest = grid.uns["shortest"]
    df_shotest = getattr(grid, "shortest")
    df_shotest["right"] = [x.right if isinstance(x, pd.Interval) else np.nan for x in df_shotest["region"]]  # NaN-safe (pandas 3 maps NaN too)
    
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
    if section is None:
        for i, s in zip(df_shotest["right"], df_shotest["solid"]):
            if (i > 0 ) & (s != 0 ): # select solid regeion and (−∞, 0] section 
                solid_ls += [0]
            else:
                solid_ls += [s]

    else:
        for i, s in zip(df_shotest["right"], df_shotest["solid"]):
            if (i != section ) & (s != 0 ): # selected section 
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
        index=[f"solid_{str(int(i))}" for i in df_solid.index],
        columns=["imagecol", "imagerow"]
    )
    df_solid_centroid["solid"] = df_solid.index
    
    
    # create AnnData object of (−∞, 0] section of indivisual tumor solid ----------------------------
    counts = df_solid.values
    solid = ad.AnnData(counts, uns={"spatial": {"data": {
        'scalefactors': {'tissue_hires_scalef': 1, 'spot_diameter_fullres': 10}, 
        'use_quality': 'hires',
        'images': {'hires':np.array([[255, 255, 255, 255]])}
    }}})
    
    # add cell and gene name
    solid.obs_names = [f"solid_{str(int(i))}" for i in df_solid.index]
    solid.var_names = df_solid.columns.tolist()
    
    # add centroid to metadata
    solid.obs = df_solid_centroid
    
    # annotation processing image to AnnData object ------------------------------------
    solid.uns["indivisual_tumor_solid"] = img_color
    # annotation of dataframe
    #solid.uns["shortest"] = df_shotest
    setattr(solid, 'shortest', df_shotest)

    return solid