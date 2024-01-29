# import
import plotly.express as px
import cv2


def equidistant_curves(
    grid, figsize=(1000, 1000)
):

    data = grid.uns["shotest_30_delineation"]
    
    fig = make_subplots()
    fig.add_trace(
        px.imshow(
            cv2.cvtColor(data, cv2.COLOR_BGR2RGB)
        ).data[0]
    )
    fig.update_layout(height=figsize[0], width=figsize[1],
                      dragmode='pan', 
                      title_text="Equidistant curves")
    
    fig.show(config = {
        'scrollZoom': True,
        'toImageButtonOptions': {
            'format': 'png', # one of png, svg, jpeg, webp
            'filename': 'preprocessed_images',
            'height': figsize[0],
            'width': figsize[1],
            'scale': 3 } # Multiply title/legend/axis/canvas sizes by this factor
    })
    