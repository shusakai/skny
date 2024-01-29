# load library
import scanpy as sc
import stlearn as st
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import cv2
from PIL import Image
import networkx as nx
import math


# function
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
      

def cv2pil(image):
    ''' OpenCV type -> PIL type '''
    new_image = image.copy()
    if new_image.ndim == 2: 
        pass
    elif new_image.shape[2] == 3: 
        new_image = cv2.cvtColor(new_image, cv2.COLOR_BGR2RGB)
    elif new_image.shape[2] == 4: 
        new_image = cv2.cvtColor(new_image, cv2.COLOR_BGRA2RGBA)
    new_image = Image.fromarray(new_image)
    return new_image

def calculate_distance(grid, pos_marker_ls, neg_marker_ls=None):
    '''
    grid : AnnData
    pos_marker_ls, neg_marker_ls: list
    '''

    ## get col and row length data from grid object --------------------------
    N_ROW = len(grid.uns["grid_yedges"]) - 1
    N_COL = len(grid.uns["grid_xedges"]) - 1
    
    ## Difine positive grid ---------------------------
    pos_series = sum([ grid.to_df()[i] for i in pos_marker_ls])
    pos_series = pd.Series([1 if i>=1 else 0 for i in pos_series], index=pos_series.index)
    
    # distract negative grid
    if neg_marker_ls != None:
        neg_series = sum([ grid.to_df()[i] for i in neg_marker_ls])
        neg_series = pd.Series([1 if i>=1 else 0 for i in neg_series], index=neg_series.index)
        pos_series = pos_series - neg_series
        pos_series = pd.Series([-1 if i==0 else i for i in pos_series], index=pos_series.index)
    
    else:
        pos_series = pd.Series([-1 if i==0 else i for i in pos_series], index=pos_series.index)
    
    # expand information of all grid
    pos_series.name = "tumor_grid"
    df_grid_tumor = pd.merge(
        pd.DataFrame(index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)]), 
        pos_series, right_index=True, left_index=True, how="left"
    ).fillna(-1)


    ## Generate data of positive grid---------------------------
    fig, ax = plt.subplots(figsize=(N_COL, N_ROW), dpi=1, tight_layout=True)
    cmap = matplotlib.cm.viridis
    cmap.set_bad('black',1.)
    cmap.set_under(color='black') 
    ax.imshow(
        np.array(df_grid_tumor["tumor_grid"]).reshape(N_COL, N_ROW).T, cmap=cmap,
        vmin=0,
        vmax=1,
    )
    ax.axis('off')
    # convert from fig to numpy
    fig.canvas.draw()
    img = np.array(fig.canvas.renderer.buffer_rgba())
    img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)
    plt.close()
    grid.uns["marker"] = img.copy() # add to grid object

    
    ## Median filter---------------------------
    #img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)
    img_med = cv2.medianBlur(img, ksize=3)
    grid.uns["marker_median"] = img_med.copy() # add to grid object

    
    ## Delineation --------------------------------
    threshold = 90 # threshold
    # convert glay scale and blur
    img_gray = cv2.cvtColor(img_med, cv2.COLOR_BGR2GRAY)
    img_blur = cv2.blur(img_gray, (3,3))
    # delineation
    ret, img_binary = cv2.threshold(img_blur, threshold, 255, cv2.THRESH_BINARY)
    contours, hierarchy = cv2.findContours(img_binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    
    contours_ = []
    # check each delineation
    for i, s in zip(contours, hierarchy[0]):
        # select top of hierarchy polygon
        #if s[-1] == -1:
        contours_ += [i]
    
    img_color_with_contours = cv2.drawContours(img_med, contours_, -1, (0,255,0), 1)
    grid.uns["marker_delineation"] = img_color_with_contours # add to grid object


    ## Converting image data to graph stracture ---------------------------------
    # OpenCV -> PIL
    image = cv2pil(img_color_with_contours)

    # extract size
    width, height = image.size
    # generate graph
    graph = nx.Graph()
    
    tumor_contour_pixel_ls = []
    tumor_pixel_ls = []
    edge_pixel_ls = []
    # ピクセルをノードに変換
    for y in range(height):
        for x in range(width):
            pixel_value = image.getpixel((x, y))
            graph.add_node((x, y), color=pixel_value)
    
            # extract tumor coutour coordination
            if pixel_value == (0, 255, 0):
                tumor_contour_pixel_ls += [(x, y)]
    
            # extract tumor coordination
            elif pixel_value == (253, 231, 36):
            #elif pixel_value == (255, 255, 0):

                tumor_pixel_ls += [(x, y)]
    
            if (y == 0) | (y == height-1) | (x == 0) | (x == width-1):
                edge_pixel_ls += [(x, y)]
    
    # get nodes
    nodes_ls = graph.nodes
    nodes_ls = list(nodes_ls)
    
    # check node inside contour
    inside_contour_ls = []
    for node in nodes_ls:
        for i in contours_:
            if cv2.pointPolygonTest(i, node, False) == 1:
                inside_contour_ls += [node]
                break
    
    # extract tumor inside contour
    tumor_pixel_ls = [i for i in tumor_pixel_ls if i in inside_contour_ls]
    
    # extract non-tumor inside contour
    #nontumor_pixel_inside_contour_ls = [i for i in inside_contour_ls if i not in tumor_pixel_ls]
    
    # 隣接ピクセル間にエッジを作成
    for y in range(height):
        for x in range(width):
            current_node = (x, y)
            neighbors = [
                (x, y - 1),  # 上
                (x, y + 1),  # 下
                (x - 1, y),  # 左
                (x + 1, y)   # 右
            ]
            # 上下左右のエッジ
            for neighbor in neighbors:
                if neighbor in graph.nodes:
                    graph.add_edge(current_node, neighbor)
    
            # 斜めのエッジ
            if (y != height) & (x != width):
                node1 = (x, y)
                node2 = (x + 1, y + 1)
    
                # ノードの座標からユークリッド距離を計算
                distance = math.sqrt((node1[0] - node2[0])**2 + (node1[1] - node2[1])**2)
                graph.add_edge(node1, node2, weight=distance)
    
                node1 = (x + 1, y)
                node2 = (x, y + 1)
    
                # ノードの座標からユークリッド距離を計算
                distance = math.sqrt((node1[0] - node2[0])**2 + (node1[1] - node2[1])**2)
                graph.add_edge(node1, node2, weight=distance)
    
    
    #remove non-tumor inside contour from graph
    #graph.remove_nodes_from(nontumor_pixel_inside_contour_ls)
    
    # delete contour at the edge of image 
    tumor_edge_contour_pixel_ls = list(set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
    
    tumor_pixel_ls = list(
        set(tumor_pixel_ls) | (set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
    )
    
    tumor_contour_pixel_ls = list(
        set(tumor_contour_pixel_ls) - (set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
    )

    #tumor_contour_pixel_ls = [i for i in tumor_contour_pixel_ls if i in list(graph.nodes)]
    
    ## Calculate distance of each nodes from tumor surface ---------------------------------------------------
    shortest_paths = nx.multi_source_dijkstra(graph, tumor_contour_pixel_ls, cutoff=16, weight="weight")
    
    # shotest length from contour
    df_shotest = pd.DataFrame.from_dict(shortest_paths[0], orient="index", columns=["euclidean"])
    
    # negative value for inside contour
    df_shotest.loc[df_shotest.index[df_shotest.index.isin(tumor_pixel_ls)]] = df_shotest.loc[df_shotest.index[df_shotest.index.isin(tumor_pixel_ls)]] * -1
    df_shotest["euclidean_round"] = df_shotest["euclidean"].round(0)
    df_nodes = pd.DataFrame(index=nodes_ls)
    
    df_shotest = pd.merge(
        df_nodes, df_shotest, right_index=True, left_index=True, how="left"
    ).fillna(np.nan)

    # debag
    #return grid
    
    ## Color scale of distance from surface ------------------------
    fig, ax = plt.subplots(figsize=(N_COL, N_ROW), dpi=1, tight_layout=True)
    cmap = matplotlib.cm.Spectral
    cmap.set_bad('black',1.)
    cmap.set_under(color='black') 
    elev_min=df_shotest["euclidean"].min()
    elev_max=df_shotest["euclidean"].max()
    mid_val=0
    
    ax.imshow(
        np.array(df_shotest["euclidean"]).reshape(N_ROW, N_COL), cmap=cmap, 
        clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=mid_val,vmin=elev_min, vmax=elev_max)
    )
    ax.axis('off')
    
    # convert from fig to numpy
    fig.canvas.draw()
    img = np.array(fig.canvas.renderer.buffer_rgba())
    plt.close()
    grid.uns["shotest"] = img # add to grid object

    ## Delineation ----------------------------------
    col_ls = []
    for i in df_shotest["euclidean_round"]:
        if i == 0:
            col_ls += [(0, 255, 0)]
        elif i < 0:
            col_ls += [(0, 255, 255)]
        else:
            col_ls += [(0, 0, 0)]
    
    col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
    grid.uns["marker_median_delineation"] = col_arr.copy() # add to grid object

    ## Delineation with 30μm interval -------------
    col_ls = []
    for i in df_shotest["euclidean_round"]:
        if i == 0:
            col_ls += [(0, 255, 0)]
        elif (i % 3 == 0) & (i > 0):
            col_ls += [(0, 0, 255)]
        elif (i % 3 == 0) & (i < 0):
            col_ls += [(0, 150, 255)]
        elif i < 0:
            col_ls += [(0, 255, 255)]
        else:
            col_ls += [(0, 0, 0)]
    
    col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
    grid.uns["shotest_30_delineation"] = col_arr.copy() # add to grid object

    ## Discretization ---------------------------------
    df_shotest["region"] = pd.cut(
        df_shotest.dropna()["euclidean"], 
        bins=list(range(-15, 16, 3))
    )
    
    for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index:
            
        col_ls = []
        for s in df_shotest["region"]:
            if s == i:
                col_ls += [(255, 255, 255)]
            else:
                col_ls += [(0, 0, 0)]
        
        col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
        grid.uns[f"shotest_region_{i}_delineation"] = col_arr # add to grid object
    
    # annotation of dataframe of shortest distance from surface
    grid.uns["shortest"] = df_shotest

    return grid
