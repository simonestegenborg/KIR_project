#!/usr/bin/env python
from gettext import install
import pandas as pd
import numpy as np


#install plotly
import plotly


#load VDJ data and create two arrays


#Source – TCRB only clonotypes
#Target – TCRA_TCRB clones (clone_id)
 
#Just two vectors with the same order of gems

def item_series(item, indexed=None):
    """
    Creates a series filled with item with indexes from indexed (if Series-like) or numerical indexes (size=indexed)
    :param item: value for filling
    :param indexed:
    :return:
    """
    if indexed is not None:
        if hasattr(indexed, 'index'):
           return pd.Series([item] * len(indexed), index=indexed.index)
        elif type(indexed) is int and indexed > 0:
            return pd.Series([item] * indexed, index=np.arange(indexed))
    return pd.Series()
 

def pivot_vectors(vec1, vec2, na_label_1=None, na_label_2=None):
    """
    Aggregates 2 vectors into a table with amount of pairs (vec1.x, vec2.y) in a cell
    Both series must have same index.
    Else different indexes values will be counted in a_label_1/na_label_2 columns if specified or ignored
    :param vec1: pd.Series
    :param vec2: pd.Series
    :param na_label_1: How to name NA column
    :param na_label_2: How to name NA row
    :return: pivot table
    """
 
    name1 = str(vec1.name)
    if vec1.name is None:
        name1 = 'V1'
 
    name2 = str(vec2.name)
    if vec2.name is None:
        name2 = 'V2'
 
    if name1 == name2:
        name1 += '_1'
        name2 += '_2'
 
    sub_df = pd.DataFrame({name1: vec1, name2: vec2})
    # FillNAs
    fill_dict = {}
    if na_label_1 is not None:
        fill_dict[name1] = na_label_1
    if na_label_2 is not None:
        fill_dict[name2] = na_label_2
    sub_df.fillna(value=fill_dict, inplace=True)
 
    sub_df = sub_df.assign(N=item_series(1, sub_df))
 
    return (
        pd.pivot_table(data=sub_df, columns=name1, index=name2, values='N', aggfunc=sum)
        .fillna(0)
        .astype(int)
    )
 


def ordered_series(indexed=None):
    """
    Creates a series 0-N; N=len(indexed) with indexes copied from indexed (if Series-like) or indexed as index
    :param indexed:
    :return:
    """
    if indexed is not None:
        if hasattr(indexed, 'index'):
            return pd.Series(np.arange(len(indexed)), index=indexed.index)
        else:
            return pd.Series(np.arange(len(indexed)), index=indexed)
 
    return pd.Series()




def sankey_plot(source, target, title='', font_size=15, width=1000, height=1000):
    """
    Sankey plot for 2 series. Names MUST be unique. Utilizes plotly instead of matplotlib
    :param source:
    :param target:
    :param title:
    :param font_size:
    :param width:
    :param height:
    :return:
    """

    """
    with warnings.catch_warnings():
        warnings.simplefilter('always', DeprecationWarning)
        warnings.warn(
            'Function `sankey_plot` will be removed in future versions.'
            ' Use `sankey_multilevel_plot` instead',
            DeprecationWarning,
        ) 
    """


    import plotly.graph_objects as go
 
    #from bioreactor.utils import ordered_series




 
    melted = pd.melt(
        pivot_vectors(target.rename('source'), source.rename('target')).reset_index(),
        id_vars='target',
   )
    melted.columns = ['source', 'target', 'value']
    uindexes = ordered_series(
        np.sort(pd.concat([melted.source, melted.target]).unique())
    )
    plot_data = pd.concat(
        [
            melted.source.map(uindexes.to_dict()),
            melted.target.map(uindexes.to_dict()),
            melted.value,
        ],
        axis=1,
    )
 
    layout = go.Layout(
        autosize=False,
        width=width,
        height=height,
        title_text=title,
        font_size=font_size,
    )
 
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=uindexes.index,
                    color="salmon",
                ),
                link=plot_data.to_dict(orient='list'),
            )
        ],
        layout=layout,
    )
 
    return fig


#vdj_data = pd.read_csv(r'/Users/simone/Desktop/MASTER/VDJ_table.csv')
vdj_data = pd.read_csv(r'/Volumes/T-cells-and-cancer/SRH group/Group members/Simone/data/metadata/vdj_table.csv',
                      index_col=0)


clone_size = vdj_data['clone_id'].value_counts()
big_clones = clone_size.index[clone_size > 10]


vdj_data_bc = vdj_data.query('clone_id in @big_clones')
TRCB_clones = vdj_data_bc["TRB_clone"].astype(str).str.replace('.0', '')
clone_id = vdj_data_bc["clone_id"].astype(str).str.replace('.0', '')


sankey_plot(TRCB_clones, clone_id ,title='Sankey Plot', font_size=15, width=1000, height=700)