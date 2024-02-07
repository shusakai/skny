# subplot
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go
import cv2
import math


def preprocessed_images(
    grid, figsize=(None, None)
):
    '''
    grid: AnnData
    figsize: tupple
    '''

    # Load data -----------------------------
    df_shotest = grid.uns["shortest"]
    
    # list 
    image_ls = [
        "marker", 
        "marker_median", 
        "marker_median_delineation", 
        "shotest_30_delineation"
    ]
    
    image_ls += [f"shotest_region_{i}_delineation" for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index]

    image_name_tp = ("Tumor", "Median filter", "Delineation", "Equidistant curves"
    ) + tuple([str(i * 10) for i in df_shotest.dropna().region.value_counts()[df_shotest.dropna().region.value_counts() > 100].sort_index().index])

    # Plot --------------------------------------
    fig = make_subplots(
        rows=math.ceil(len(image_name_tp) / 4), cols=4, vertical_spacing=0.05, horizontal_spacing=0.02, shared_yaxes='all',shared_xaxes='all',
        row_heights=[0.5] * math.ceil(len(image_name_tp) / 4), 
        column_widths=[0.5, 0.5, 0.5, 0.5], #vertical_spacing=0.02
        subplot_titles=image_name_tp, 
        x_title='X-axis',
        y_title='Y-axis',
    
    )
    
    for e, i in enumerate(image_ls):
        result_image = grid.uns[i]
        result_image = cv2.cvtColor(result_image, cv2.COLOR_BGR2RGB)
    
        px.imshow(result_image)
        fig.add_trace(
            px.imshow(
                result_image,
            ).data[0],
            row=int(e/4)+1, col=e%4+1
        )
        
    fig.update_layout(height=figsize[0], width=figsize[1],
                      dragmode='pan', 
                      title_text="Preprocessed images")
    
    fig.show(config = {
        'scrollZoom': True,
        'toImageButtonOptions': {
            'format': 'png', # one of png, svg, jpeg, webp
            'filename': 'preprocessed_images',
            'height': figsize[0],
            'width': figsize[1],
            'scale': 3 } # Multiply title/legend/axis/canvas sizes by this factor
    })
    
    #fig.write_html("merge_image.html", config = {
    #    'scrollZoom': True,
    #})



