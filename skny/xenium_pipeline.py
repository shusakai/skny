## load library ---------------------------------------------------------
import stlearn as st
import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import cv2
import pandas as pd
import scipy.stats
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots
from PIL import Image
import networkx as nx
import matplotlib.pyplot as plt
import math
import sys
import os
import matplotlib.colors as colors


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


### main ###
# arg
args = sys.argv
cell_feature_matrix_h5 = args[1]
cells_csv = args[2]
pos_marker = args[3] # marker
pos_marker_ls = pos_marker.split(";")
neg_marker = args[4] # marker
neg_marker_ls = neg_marker.split(";")
input_dir = args[5]
output_dir = args[6]


## Load data and gridding ------------------------------------------------------------
# load the data
adata = st.ReadXenium(feature_cell_matrix_file=f"{os.path.join(input_dir, cell_feature_matrix_h5)}",
                     cell_summary_file=f"{os.path.join(input_dir, cells_csv)}",
                     library_id="aaa",
                     image_path=None,
                     scale=1,
                     spot_diameter_fullres=10
                     )

# Gridding at 10μm interval
N_COL = int((adata.obs.imagecol.max() - adata.obs.imagecol.min()) / 10)
N_ROW = int((adata.obs.imagerow.max() - adata.obs.imagerow.min()) / 10)
# Gridding
grid = st.tl.cci.grid(adata, n_row=N_ROW, n_col=N_COL)


## Extract tumor grids ------------------------------------------------------------

# sum positive grid
pos_series = sum([ grid.to_df()[i] for i in pos_marker_ls])
pos_series = pd.Series([1 if i>=1 else 0 for i in pos_series], index=pos_series.index)


# distract negative grid
if neg_marker != "":
    neg_series = sum([ grid.to_df()[i] for i in neg_marker_ls])
    neg_series = pd.Series([1 if i>=1 else 0 for i in neg_series], index=neg_series.index)
    pos_series = pos_series - neg_series
    pos_series = pd.Series([-1 if i==0 else i for i in pos_series], index=pos_series.index)

else:
    pos_series = pd.Series([-1 if i==0 else i for i in pos_series], index=pos_series.index)

# all grid
pos_series.name = "tumor_grid"
df_grid_tumor = pd.merge(
    pd.DataFrame(index=["grid_" + str(i+1) for i in range(N_ROW * N_COL)]), 
    pos_series, right_index=True, left_index=True, how="left"
).fillna(-1)

# save figure
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
plt.savefig(
    os.path.join(output_dir, "marker.png")
)


## Delineation of tumor surface-------------------------------------------------------------------------
# median filter
img_path = os.path.join(output_dir, "marker.png")
img = cv2.imread(img_path)
img_med = cv2.medianBlur(img, ksize=3)

# save
cv2.imwrite(
    os.path.join(output_dir, "marker_median.png"), img_med
)

# delineation
img_path = os.path.join(output_dir, "marker_median.png")
threshold = 90 #二値化に用いる閾値

# glay scale
img_color = cv2.imread(img_path)
img_gray = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
img_blur = cv2.blur(img_gray,(3,3))

# delineation
ret, img_binary = cv2.threshold(img_blur, threshold, 255, cv2.THRESH_BINARY)
contours, hierarchy = cv2.findContours(img_binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

contours_ = []
# check each delineation
for i, s in zip(contours, hierarchy[0]):
    # select top of hierarchy polygon
    if s[-1] == -1:
        contours_ += [i]

img_color_with_contours = cv2.drawContours(img_color, contours_, -1, (0,255,0), 1)
# save
cv2.imwrite(
    os.path.join(output_dir, "marker_delineation.png"), img_color_with_contours
)


## Conversion from grids to graph -------------------------------------------
# 画像を読み込む
image = Image.open(
    os.path.join(output_dir, "marker_delineation.png")
)

# 画像サイズを取得
width, height = image.size

# 新しい有向グラフを作成
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
nontumor_pixel_inside_contour_ls = [i for i in inside_contour_ls if i not in tumor_pixel_ls]

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
graph.remove_nodes_from(nontumor_pixel_inside_contour_ls)

# delete contour at the edge of image 
tumor_edge_contour_pixel_ls = list(set(tumor_contour_pixel_ls) & set(edge_pixel_ls))

tumor_pixel_ls = list(
    set(tumor_pixel_ls) | (set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
)

tumor_contour_pixel_ls = list(
    set(tumor_contour_pixel_ls) - (set(tumor_contour_pixel_ls) & set(edge_pixel_ls))
)

## Calculate distance of each nodes from tumor surface ----------------------------------------------
# Calculate distance of each nodes from tumor surface
# ある一点からの最短距離を計算（cutoff=16）
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
     
# save figure
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
plt.savefig(
    os.path.join(output_dir, "shotest.png")
)


# delineation of contours
img_path = os.path.join(output_dir, "shotest.png")
threshold = 90 #二値化に用いる閾値

# glay scale
img_color = cv2.imread(img_path)

img_color_with_contours = cv2.drawContours(img_color, contours_, -1, (0,0,0), 1)
# save
cv2.imwrite(
    os.path.join(output_dir, "shotest_delineation.png"), img_color_with_contours
)

# plot each region
col_ls = []
for i in df_shotest["euclidean_round"]:
    if i == 0:
        col_ls += [(0, 255, 0)]
    elif i < 0:
        col_ls += [(0, 255, 255)]
    else:
        col_ls += [(0, 0, 0)]

col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
#img_color_with_contours = cv2.drawContours(col_arr, contours_, -1, (0,255,0), 1)

# save
cv2.imwrite(
    os.path.join(output_dir, "marker_median_delineation.png"), col_arr
)

# plot each region
col_ls = []
for i in df_shotest["euclidean_round"]:
    if i == 0:
        col_ls += [(0, 255, 0)]
    elif i % 3 == 0:
        col_ls += [(0, 0, 255)]
    else:
        col_ls += [(0, 0, 0)]

col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)

# save
cv2.imwrite(
    os.path.join(output_dir, "shotest_30_delineation.png"), col_arr
)


## Region separation--------------------------------------------------------------
# 離散化
df_shotest["region"] = pd.cut(
    df_shotest.dropna()["euclidean"], 
    bins=list(range(-15, 16, 3))
)

# plot each region
for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index:
        
    col_ls = []
    for s in df_shotest["region"]:
        if s == i:
            col_ls += [(255, 255, 255)]
        else:
            col_ls += [(0, 0, 0)]
    
    col_arr = np.array(col_ls, dtype=np.uint8).reshape(N_ROW, N_COL, 3)
    
    # save
    cv2.imwrite(
        os.path.join(output_dir, f"shotest_region_{i}_delineation.png"), col_arr
    )

## Region annotation to expression data of Xenium --------------------------------------------
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

df_grid_region = pd.merge(
    df_grid, df_region, 
    right_index=True, left_index=True, how="left"
)

## Clustering -------------------------------------------------------------------------
df_grid_region = df_grid_region.dropna()
df_grid_region_gene = df_grid_region.groupby("region").mean().T
region_ls = df_shotest.region.value_counts()[df_shotest.region.value_counts() > 100].index
df_grid_region_gene = df_grid_region_gene.loc[:,region_ls]
df_grid_region_gene.columns = [i * 10 for i in df_grid_region_gene.columns]
df_grid_region_gene = df_grid_region_gene.T.sort_index().T



# draw clustering map
g = sns.clustermap(
    df_grid_region_gene.apply(lambda x: scipy.stats.zscore(x), axis=1).dropna(), 
    figsize=(8, 12), cmap="viridis", col_cluster=False, method="ward"
)
ax = g.ax_heatmap
ax.set_ylabel("")
#ax.tick_params(axis="x", labelrotation=45)

# 行と列のクラスタリングの順番を抽出
row_order = g.dendrogram_row.reordered_ind
df_grid_region_gene_std = df_grid_region_gene.apply(lambda x: scipy.stats.zscore(x), axis=1).iloc[row_order, :]
df_grid_region_gene = df_grid_region_gene.iloc[row_order, :]

hover_text = []
for e1, i in enumerate(df_grid_region_gene_std.columns):
    temp_ls = []
    for e2, gene in enumerate(df_grid_region_gene_std.index):
        temp_ls += [
            f"""{gene} {i}<br> 
Mean density={(df_grid_region_gene.iloc[e2, e1] / (10 * 10)).round(3)}/μm<sup>2</sup><br> 
Z-score={df_grid_region_gene_std.iloc[e2, e1].round(3)}"""
        ]
    hover_text += [temp_ls]
hover_text = np.array(hover_text).T

# zoomを制限
layout = go.Layout(
    xaxis=dict(
        fixedrange=True  # y軸のズームを固定
    )
)

x_col = [int(float(i.split(", ")[0][1:])) for i in df_grid_region_gene_std.columns.astype(str)] +\
[int(float(df_grid_region_gene_std.columns.astype(str)[-1].split(", ")[1][:-1]))]

fig = go.Figure(
    data=go.Heatmap(
        z=df_grid_region_gene_std, 
        x=x_col,
        #x=[-60, -30, 0, 30, 60, 90, 120, 150],
        y=df_grid_region_gene_std.index,
        colorscale = 'Viridis', 
        hoverinfo="text",  # ホバー時にテキストを表示
        text=hover_text  # ホバーテキストを指定
    ), layout=layout
)

fig.update_xaxes(
    tickvals=x_col,  # 表示する位置のインデックス
    ticktext=x_col  # ラベルのリスト
)
fig.update_layout(title='Clustering based on transcript density of each region')

fig.update_layout(xaxis=dict(title='Region'))

plotly.io.write_html(
    fig, os.path.join(output_dir, "grid_region_gene_std.html")
)

# merge images----------------------------------------------------------------------------------
# list 
image_ls = [
    "marker.png", 
    "marker_median.png", 
    "marker_median_delineation.png", 
    "shotest.png", 
    "shotest_30_delineation.png"
]

image_ls += [f"shotest_region_{i}_delineation.png" for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index]

image_name_tp = ("Tumor", "Median filter", "Auto-delineation", "Distance", "Distance /30μm"
) + tuple([str(i * 10) for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index])

fig = make_subplots(
    rows=math.ceil(len(image_name_tp) / 4), cols=4, vertical_spacing=0.05, horizontal_spacing=0.02, shared_yaxes='all',shared_xaxes='all',
    row_heights=[0.5] * math.ceil(len(image_name_tp) / 4), 
    column_widths=[0.5, 0.5, 0.5, 0.5], #vertical_spacing=0.02
    subplot_titles=image_name_tp, 
    x_title='X-axis',
    y_title='Y-axis',

)

for e, i in enumerate(image_ls):
    result_image = cv2.imread(
        os.path.join(output_dir, i)
    )
    result_image = cv2.cvtColor(result_image, cv2.COLOR_BGR2RGB)

    px.imshow(result_image)
    fig.add_trace(
        px.imshow(
            result_image,
        ).data[0],
        row=int(e/4)+1, col=e%4+1
    )
    
fig.update_layout(height=1200, width=1000,
                  dragmode='pan', 
                  title_text="Delineation, ellipse fitting, and separation of region")

fig.write_html(
    os.path.join(output_dir,"merge_image.html"), 
    config = {
    'scrollZoom': True,
})

# save dataframe
df_grid_region_gene_std.to_csv(
    os.path.join(output_dir, "df_grid_region_gene_std.txt"), sep="\t"
)

print("Done: plot figures")
