# quick visualize functions
from . import compartment
import pandas as pd
import numpy as np

def plotPsCurve(mcoolsPath:list,celltypeNames:list,chroms:list,resolution=100000,title="P(s) curve",plotType="interaction"):
    """
    plotPsCurve function take bin and value
    """
    import plotly.express as px
    from IPython.display import Image

    #Calculate P(s) data, get a 3 column pd.DataFrame with (bin,resolution,celltype)
    psDataAll = []
    for i in range(len(mcoolsPath)):
        psDataAll.append(compartment.getPsData(mcoolsPath[i],["chr"+str(i+1) for i in range(len(chroms))],resolution=resolution,celltype=celltypeNames[i])) 
    merged = pd.concat(psDataAll)

    data =  pd.merge(merged,merged.groupby("celltype").sum(),how="left",on="celltype").assign(prob= lambda df: df.aveCount_x/df.aveCount_y)

    fig = px.line(x=data["bin_x"]*resolution,y=data["prob"],color=data["celltype"],title=title,log_x=True,log_y=True).update_layout(template='simple_white')
    fig.update_layout(width=800,height=600)
    fig.update_layout(xaxis_title="Genomic Distance(bp)",
                    yaxis_title="Contact Probability")
    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))

def plotMatrix(matrix:np.ndarray,if_log=False,title="Matrix",plotType="static"):
    """
    plotMatrix function for plot hic contact matrix
    """
    import plotly.express as px 
    from IPython.display import Image

    if(if_log == True): 
        matrix = np.log10(matrix)

    fig = px.imshow(matrix,color_continuous_scale=px.colors.sequential.Bluered)
    fig = fig.update_layout(template='simple_white').update_layout(width=800,height=600)

    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))

def genomecoord2human(n):
    symbols = ('bp','Kb','Mb','Gb')
    exp = int(np.log10(n)/3)
    return str(round(n/(1000**exp),2))+symbols[exp]

def plotRegionFromMCOOLS(filepath:str,resolution:int,genome_coord1:str,genome_coord2=None,if_log=False,balance=False,title="Matrix",plotType="static"):
    """
    plotMatrix function for plot hic contact matrix
    """
    import plotly.express as px 
    from IPython.display import Image
    import re
    import cooler
    import numpy as np

    cool = filepath+"::/resolutions/"+str(resolution)

    if (genome_coord2 == None):
        genome_coord2 = genome_coord1
    c = cooler.Cooler(cool)
    matrix = c.matrix(balance=balance).fetch(genome_coord1,genome_coord2).astype("double")

    if(if_log == True): 
        matrix = np.log10(matrix+1)

    #fig = px.imshow(matrix,color_continuous_scale=px.colors.sequential.Bluered)
    fig = px.imshow(matrix)
    fig = fig.update_layout(title=title)
    fig = fig.update_layout(template='simple_white').update_layout(width=600,height=500)
    fig = fig.update_layout(xaxis_title=genome_coord2,yaxis_title=genome_coord1)

    #manually change axis
    posx = re.split("[:-]",genome_coord2)
    xvals = np.percentile([np.round(i) for i in range(0,matrix.shape[1])],(0,25,50,75,100),interpolation='midpoint')
    xtexts = xvals*resolution + int(posx[1].replace(",","")) + resolution/2
    xtexts = [genomecoord2human(i) for i in xtexts]

    posy = re.split("[:-]",genome_coord1)
    yvals = np.percentile([np.round(i) for i in range(0,matrix.shape[0])],(0,25,50,75,100),interpolation='midpoint')
    ytexts = yvals*resolution + int(posy[1].replace(",","")) + resolution/2
    ytexts = [genomecoord2human(i) for i in ytexts]

    fig = fig.update_xaxes(ticktext = xtexts,tickvals = xvals).update_yaxes(ticktext = ytexts,tickvals = yvals)

    # static plot have better performance in jupyter
    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))