from flask import Flask, render_template, send_from_directory, request, Response

import matplotlib
import numpy as np

from datastore import Datastore
matplotlib.use('agg')
import matplotlib.pyplot as plt

import seaborn as sns

import plotly
import plotly.graph_objects as go

import os

import argparse
parser = argparse.ArgumentParser(description='Flask-app dashboard for pande-mic.')
parser.add_argument('--summary_path', help='full path to summary directory', required=True)
parser.add_argument('--config_path', help='full path to config directory', required=True)
parser.add_argument('--reference', help='path to reference, if it is not in config directory named ref.fasta')
parser.add_argument('--amplicons', help='path to amplicon scheme, if it is not in config directory named amplicons.json')
parser.add_argument('--run_name', help='run name or ID', default="unknown", required=True)
parser.add_argument('--smoothing_window_size', help='window size used for smoothing of coverage graphs', default="10")
parser.add_argument('--smoothing_window_step', help='step between smoothing windows for coverage graphs (calculating madians for each position might be slow)', default="10")
parser.add_argument('--low_cutoff', help='coverage value that is considered too low', default = "20")
parser.add_argument('--high_cutoff', help='coverage value that is considered too high', default = "3000")
parser.add_argument('--coverage_cutoff', help='cut values higher than high_cutoff from the graph', dest='coverage_cutoff', action='store_true')
parser.add_argument('--port', default = 5000, required=False)
parser.add_argument('--debug', help='run app in debug mode', action='store_true', required=False)
args = parser.parse_args()


summary_full_path = args.summary_path
config_full_path = args.config_path
run_name = args.run_name
smooth = int(args.smoothing_window_size)
smooth_step = int(args.smoothing_window_step)
cutoff_low = int(args.low_cutoff)
cutoff_high = int(args.high_cutoff)
if args.reference:
    reference = args.reference
else:
    reference = config_full_path+"ref.fasta"
if args.amplicons:
    amplicons = args.amplicons
else:
    amplicons = config_full_path+"amplicons.json"
coverage_cutoff_isset = False
if args.coverage_cutoff:
    coverage_cutoff_isset = True
debug = False
if args.debug:
    debug = True
specified_port = args.port

app = Flask(__name__)

datastore = Datastore(summary_full_path, config_full_path, run_name, reference, amplicons, smooth, smooth_step)

@app.route('/')
def index():
    """render the index page"""
    global datastore
    if(os.path.isfile(datastore.summary_path+"basecalled/count.csv")): # if the summary file from base calling exists
        datastore.read_stats()
        return render_template('index.html',run_name=datastore.get_run_name())
    else: # no data yet
        return render_template('index_preload.html', run_name=datastore.get_run_name())

@app.route('/stats')
def stats():
    """read number of reads and bases and render template"""
    global datastore
    datastore.read_stats()
    return render_template('stats.html', reads=datastore.get_num_of_reads() , bases=datastore.get_num_of_bases())

@app.route('/variants_table')
def variants_table():
    """get variant summary and render template with variants table"""
    return render_template('variants_table.html', variants=datastore.get_variant_summary())

@app.route('/summary_graphs')
def summary_graphs():
    """read the summary statistics for all barcodes and render the graphs in the summary part of the page"""
    global datastore
    return render_template('summary.html', barcodes=datastore.get_barcodes(), colors=datastore.get_colors(), processed_reads=datastore.get_processed_reads())

@app.route('/force_reload_summary_data')
def force_reload_summary_data():
    """force reload of summary data"""
    global datastore
    datastore.read_stats()
    datastore.read_base_counts()
    datastore.read_barcode_summary()
    datastore.read_variants()
    datastore.variant_counts()
    return Response("OK", status=200)

@app.route('/barcodes')
def barcodes():
    """get data for individual barcodes and render the template containing list of them"""
    global datastore
    return render_template('barcodes.html', barcodes=datastore.get_barcodes(), colors=datastore.get_colors(), processed_reads=datastore.get_processed_reads(), mapped_reads=datastore.get_mapped_reads(), avg_read_lengths = datastore.get_avg_read_lengths())


@app.route('/summary_file/<path:fname>')
def static_file(fname):
    """load static file from the server"""
    global datastore
    return send_from_directory(datastore.summary_path, fname)

@app.route('/read_depth_plotly')
def read_depth_plotly():
    """draw a coverage graph using plotly library"""
    global datastore

    def draw_amplicons():
        amplicon_pools = datastore.divide_amplicons_into_pools()

        row = 0
        for pool in amplicon_pools:
            row+=1
            for amplicon in pool:
                top_y = (min_value-row*padding) 
                bottom_y = (min_value-(row+1)*padding)
                left_x = amplicon[0]
                right_x = amplicon[1]
                fig.add_shape(type="rect",
                    x0=left_x, y0=top_y, x1=right_x, y1=bottom_y,
                    line=dict(
                        color="RoyalBlue",
                        width=1,
                    ),
                    fillcolor="RoyalBlue",
                )

    def coverage_cutoff():
        nonlocal padding
        if coverage_cutoff_isset:
            lowest = min_value
            highest = min(cutoff_high, max_value+padding)
            if highest > lowest:
                padding = (highest - lowest)/50
                fig.update_yaxes(range=[lowest-4*padding, highest])
            else: #low coverage (max_value < low coverage threshold)
                padding = (max_value-min_value)/50
                fig.update_yaxes(range=[min_value-4*padding, max_value+padding])
        else:
            padding = (max_value-min_value)/50
            fig.update_yaxes(range=[min_value-4*padding, max_value+padding])

    def add_variant_vlines():
        """add a vertical line to the graph for each SNP found for a variant"""
        SNPs = datastore.get_variant_SNPs(barcode)
        for snp in SNPs:
            pos = snp["position"]
            fig.add_vline(x=pos, line_color='red', line_width = 0.5)


    barcode = request.args.get('barcode', default = None, type = str)
    if barcode!=None:
        datastore.reload_barcode_data(barcode)
    data = datastore.get_coverage_sliding_window_median()
    barcodes = datastore.get_barcodes()
    data_colors = datastore.get_colors()


    fig = go.Figure()
    fig.update_layout(
        title="Coverage per position ("+str(smooth)+" bases median with step "+str(smooth_step)+")",
        xaxis_title="position",
        yaxis_title="coverage",
        legend_title="barcodes",
        font=dict(
            family="Arial",
            size=12,
            color="Black"
        ),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=0.05, gridcolor='#dbdbdb')
    fig.update_yaxes(showgrid=True, gridwidth=0.05, gridcolor='#dbdbdb')
    fig.add_hline(y=cutoff_high, line_color="red", line_width = 0.5, line_dash="dash")
    fig.add_hline(y=cutoff_low, line_color="red", line_width = 0.5, line_dash="dash")
    
    min_value = 99999999999
    max_value = 0
    if barcode==None: # draw graph containing all barcodes
        for barcode in barcodes:
            if barcode!="unclassified" or len(barcodes)==1:
                d = np.array(data[barcode])
                x_data,y_data = d.T
                mi = min(y_data)
                if mi<min_value:
                    min_value = mi
                ma = max(y_data)
                if ma>max_value:
                    max_value = ma
                fig.add_trace(go.Scatter(x=x_data, y=y_data,
                        mode='lines',
                        name=barcode,
                        line=dict(color=data_colors[barcode])))
    else: # graph only for one barcode
        d = np.array(data[barcode])
        x_data,y_data = d.T
        min_value = min(y_data)
        max_value = max(y_data)
        fig.add_trace(go.Scatter(x=x_data, y=y_data,
                mode='lines',
                name=barcode,
                line=dict(color=data_colors[barcode])))
        add_variant_vlines()


    padding = (max_value-min_value)/50
#    padding = 0
    coverage_cutoff()

    draw_amplicons()

#    fig.update_yaxes(type="log")
#    fig.update_yaxes(range=[1, 50000])

    return plotly.io.to_html(fig, full_html=False)


@app.route('/read_stats')
def read_stats():
    """read summary stats"""
    global datastore
    datastore.read_stats()

@app.route('/mapped_and_processed.png')    
def mapped_and_processed_reads():
    """create figure for mapped and processed reads graph and save it to a png file"""
    global datastore
    sorted_barcodes = datastore.get_barcodes()
    if len(sorted_barcodes)==1:
        fig_y_size = 2
    elif len(sorted_barcodes)<5:
        fig_y_size = 3
    elif len(sorted_barcodes)<10:
        fig_y_size = 5
    elif len(sorted_barcodes)<96:
        fig_y_size = 0.2*len(sorted_barcodes)
    else:
        fig_y_size = 10
    fig, ax = plt.subplots(figsize=(5,fig_y_size))
    draw_mapped_and_processed_reads(ax)
    fig.tight_layout()
    fig.savefig("./static/img/mapped_and_processed_reads.png")
    plt.close()
    return "/static/img/mapped_and_processed_reads.png"

def draw_mapped_and_processed_reads(ax):
    """draw a bar plot with num of processed reads and mapped for each sample, read from barcode_summary/summary.csv, aligned/mapped.csv"""
    global datastore
    data_colors = datastore.get_colors()
    sorted_barcodes = datastore.get_barcodes()
    processed_reads = datastore.get_processed_reads()
    mapped_reads = datastore.get_mapped_reads()
    barcodes=[]
    counts=[]
    colors=[]
    counts_mapped=[]
    colors_mapped=[]
    for barcode in sorted_barcodes:
        if barcode!="unclassified" or len(sorted_barcodes)==1:
            barcodes.append(barcode)
            counts.append(processed_reads[barcode])
            counts_mapped.append(mapped_reads[barcode])
            colors.append('#808080')
            colors_mapped.append(data_colors[barcode])
    ax.margins(y=0)
    ax.barh(barcodes, counts,  color=colors, align='center')
    ax.barh(barcodes, counts_mapped,  color=colors_mapped, align='center')
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_yticks(np.arange(len(barcodes)), labels=barcodes, fontsize = 7)
    ax.set_ylabel('Barcode')
    ax.set_xlabel('# of reads')
    ax.set_title('Processed (grey) and mapped reads (color) \n per sample')    

@app.route('/barcode_amplicon.png')
def barcode_amplicon_graph():
    """create a figure for barcode - amplicon heatmap and save it to png"""
    global datastore
    sorted_barcodes = datastore.get_barcodes()
    if len(sorted_barcodes)==1:
        fig_y_size = 3
    elif len(sorted_barcodes)<5:
        fig_y_size = 3
    elif len(sorted_barcodes)<10:
        fig_y_size = 5
    elif len(sorted_barcodes)<96:
        fig_y_size = 0.3*len(sorted_barcodes)
    else:
        fig_y_size = 10
    fig, ax = plt.subplots(figsize=(5,fig_y_size))
    draw_barcode_amplicon_graph(fig,ax,250)
    fig.savefig("./static/img/barcode_amplicon.png")
    plt.close()
    return "/static/img/barcode_amplicon.png"

def draw_barcode_amplicon_graph(fig,ax,max_cov):
    """draw a baecode-amplicon heatmap"""
    global datastore
    amplicon_captions, barcode_captions, barcode_amplicon = datastore.create_amplicon_data()
    colors = ["#000000","#830029","#e40017","#ff8b02","#74b643"]
    ax = sns.heatmap( barcode_amplicon , linewidths = .5 , cmap = colors, vmin = 0, vmax = max_cov) # .. square = True
    # Show all ticks and label them with the respective list entries
    ax.set_yticks(np.arange(len(barcode_captions))+0.5, labels=barcode_captions, fontsize = 10)
    ax.set_xticks(np.arange(len(amplicon_captions))+0.5, labels=amplicon_captions, fontsize = 10)
    # Rotate the tick labels and set their alignment.
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    for tick in ax.get_yticklabels():
        tick.set_rotation(0)
    # Loop over data dimensions and create text annotations.
    #for i in range(len(barcodes)):
    #    for j in range(len(amplicons)):
    #        text = ax.text(j, i, barcode_amplicon[i, j],
    #                       ha="center", va="center", color="w")
    ax.set_title("Avg amplicon coverage for barcodes")
    fig.tight_layout()

if __name__ == "__main__":
    if debug:
        app.run(use_reloader=True, debug=True)  
    else:
        from waitress import serve
        serve(app, host="0.0.0.0", port=specified_port)
