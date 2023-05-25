import csv
import json
import random

import numpy as np
from statistics import median
import time


class Datastore:
    def __init__(self, summary_path, config_path, run_name, reference, ampliconsjson,  smoothing_window_size, smoothing_window_step):
        self.summary_path = summary_path
        self.config_path = config_path
        self.sorted_barcodes = []
        self.barcode_colors = {}
        self.reads = 0
        self.bases = 0
        self.processed_reads = {}
        self.barcode_bases = {}
        self.ref_length = 0
        self.ref = ""
        self.ref_name = ""
        self.pos_coverage = {}
        self.pos_coverage_sliding_window_median = {}
        self.run_name = run_name
        self.smoothing_window_size = smoothing_window_size
        self.smoothing_window_step = smoothing_window_step
        self.variants = {}
        self.variant_summary = {}
        self.avg_read_lengths = {}
        self.SNPs = {}
        self.variant_SNPs = {}
        self.amplicons = []
        self.amplicon_pools = []
        self.mapped_reads = {}
        self.color_picker_pos = 0
        self.reference = reference
        self.ampliconsjson = ampliconsjson
        self.colors = ['#e82c17', '#e8b417', '#2d851d', '#1d857e', '#162b73', '#411a6e', '#6e1a66', '#4f0814', '#3D81B8', '#33995F', '#A8384F', '#C7C557', '#C69453', '#4E2267', '#302A7E', '#1D4258', '#194D1A', '#8244C1', '#5C391F', '#602038', '#67C95E', '#3E55BB', '#267329', '#D98C8C', '#4FC4C0', '#B83DB6', '#090A1B', '#627CCB', '#2A687E']

####################################################################################
##### BASIC STATS ###########
####################################################################################
    
    def read_stats(self):
        self.read_ref_from_file()
        self.read_total_reads_and_bases()
        self.read_barcode_summary()
        self.read_mapped_reads()

    def get_run_name(self):
        return self.run_name

####################################################################################

    def get_ref_length(self, force_reload=False):
        if self.ref_length == 0 or force_reload:
            self.read_ref_from_file()
        return self.ref_length

    def read_ref_from_file(self):
        with open(self.reference, 'r') as f:
            self.ref_name = f.readline() #reference genome
            self.ref = f.readline().strip()[1:]
            self.ref_length = len(self.ref)

####################################################################################

    def get_num_of_reads(self, force_reload=False):
        if self.reads == 0 or force_reload:
            self.read_total_reads_and_bases()
        return self.reads

    def get_num_of_bases(self, force_reload=False):
        if self.bases == 0 or force_reload:
            self.read_total_reads_and_bases()
        return self.bases

    def read_total_reads_and_bases(self):
        with open(self.summary_path+"basecalled/count.csv", 'r') as f:
            line = f.readline()
            l = line.split(",") #num of reads, num of bases (total)
            self.reads = int(l[0])
            self.bases = int(l[1])

####################################################################################
########### COVERAGE  ###########
####################################################################################

    def get_pos_coverage(self, force_reload=False):
        if not self.pos_coverage or force_reload:
            self.read_base_counts()
        return self.pos_coverage
    
    def get_barcode_pos_coverage(self, barcode, force_reload=False):
        if not self.pos_coverage or not barcode in self.pos_coverage or force_reload:
            self.read_base_counts()
        return self.pos_coverage[barcode]

    def get_coverage_sliding_window_median(self, force_reload=False):
        if not self.pos_coverage_sliding_window_median or force_reload:
            self.read_base_counts()
        return self.pos_coverage_sliding_window_median


    def read_base_counts(self):
        start = time.time()
        if not self.sorted_barcodes: #dict. empty
            self.read_barcode_summary()
        for barcode in self.sorted_barcodes:
            self.read_base_counts_for_barcode(barcode)
        end = time.time()
        print("csv processing finished in:",end-start)

    def read_base_counts_for_barcode(self,barcode):
        cov, cov_sliding_w= self.read_base_counts_file(barcode+".csv")
        self.pos_coverage[barcode] = cov
        self.pos_coverage_sliding_window_median[barcode] = cov_sliding_w

    def read_base_counts_file(self,file):
        print("processing file "+file+"...")
        cov = [] #[pos, coverage]
        cov_sliding_w = []
        with open(self.summary_path+"coverage_per_base/"+file, 'r') as f:
            csvreaded = csv.reader(f)
            header = next(csvreaded) #position, A, C, G, T
            for line in csvreaded:
                c = int(line[1])+int(line[2])+int(line[3])+int(line[4])
                cov.append([int(line[0]), c])
        cov_sliding_w = self.median_window(cov,self.smoothing_window_size, step = self.smoothing_window_step)
        return cov, cov_sliding_w

    def median_window(self,data, w_size, step = 1):
        cov_sliding_w_xb = [] #[pos, soverage for sliding window of length 10]
        data_length = len(data)
        for i in range(0, data_length-w_size, step):
            cov_sliding_w_xb.append([i, median([int(w[1]) for w in data[i:(i+w_size)]])] )
        return cov_sliding_w_xb

    """
    def avg_window(self,data, w_size):
        cov_sliding_w_xb = [] #[pos, soverage for sliding window of length 10]
        data_length = len(data)
        for i in range(0, data_length-w_size):
            cov_sliding_w_xb.append([i, np.average([int(w[1]) for w in data[i:(i+w_size)]])] )
        return cov_sliding_w_xb
    """

####################################################################################
######### BARCODES ###########
####################################################################################

    def get_barcodes(self, force_reload=False):
        if len(self.sorted_barcodes)==0 or force_reload:
            self.read_barcode_summary()
        return self.sorted_barcodes

    def get_colors(self, force_reload=False):
        if not self.barcode_colors or force_reload:
            self.read_barcode_summary()
        return self.barcode_colors

    def get_processed_reads(self, force_reload=False):
        if not self.processed_reads or force_reload:
            self.read_barcode_summary()
        return self.processed_reads

    def get_mapped_reads(self, force_reload=False):
        if not self.mapped_reads or force_reload:
            self.read_mapped_reads()
        return self.mapped_reads

    def get_avg_read_lengths(self, force_reload=False):
        if not self.avg_read_lengths or force_reload:
            self.read_barcode_summary()
        return self.avg_read_lengths

    def read_barcode_summary(self):
        with open(self.summary_path+"barcode_summary/summary.csv", 'r') as f:
            for line in f:
                line = line.strip()
                l = line.split(",") # num of reads, num of bases, barcode
                barcode = l[2]
                self.add_barcode(barcode)
                self.processed_reads[barcode]=int(l[0])
                self.barcode_bases[barcode]=int(l[1])
                if self.barcode_bases[barcode]!= 0:
                    self.avg_read_lengths[barcode]=round(float(l[1])/float(l[0]), 1)
                else:
                    self.avg_read_lengths[barcode]=0

    def read_mapped_reads(self):
        with open(self.summary_path+"aligned/mapped.csv", 'r') as f:
            for line in f:
                line = line.strip()
                l = line.split(",") # barcode, num of reads
                barcode = l[0]
                self.add_barcode(barcode)
                self.mapped_reads[barcode]=int(l[1])

    def add_barcode(self,barcode):
        if not barcode in self.sorted_barcodes: # update color and add barcode only if it is not already there
            self.sorted_barcodes.append(barcode)
            self.barcode_colors[barcode] = self.choose_color()#we dont want the colors to be completely random

    def choose_color(self):
        if self.color_picker_pos>=len(self.colors)-1:
            color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
        else:
            color = self.colors[self.color_picker_pos]
            self.color_picker_pos += 1 
        return color


#######################################################################

    def reload_barcode_data(self, barcode):
        self.read_base_counts_for_barcode(barcode)
        self.read_variant_SNPs_for_barcode(barcode)

    def read_variant_SNPs_for_barcode(self, barcode):
        self.variant_SNPs[barcode] = []
        self.recursive_variant_SNPs(self.variants[barcode], barcode)

####################################################################################
###### VARIANTS ######
####################################################################################

    def get_variant_summary(self):
        if not self.variant_summary:
            self.variant_counts()
        return self.variant_summary

    def get_variant_SNPs(self, barcode):
        if not self.variant_SNPs:
            self.read_variant_SNPs()
        return self.variant_SNPs[barcode]

    def variant_counts(self):
        self.variant_summary = {}
        if not self.variants:
            self.read_variants()
        for barcode in self.variants:
            self.recursive_variant_summary(self.variants[barcode], self.variant_summary, barcode)

    def recursive_variant_summary(self, variants, variant_summary_parent, barcode):
        for variant in variants:
            vname = variant["name"]
            if not vname in variant_summary_parent:
                variant_summary_parent[vname]={}
                variant_summary_parent[vname]["count"] = 1
                variant_summary_parent[vname]["barcodes"] = [barcode]
                variant_summary_parent[vname]["children"] = {}
            else:
                variant_summary_parent[vname]["count"]+=1
                variant_summary_parent[vname]["barcodes"].append(barcode)
            
            self.recursive_variant_summary(variant["subs"], variant_summary_parent[vname]["children"], barcode)

    def read_variant_SNPs(self):
        for barcode in self.sorted_barcodes:
            self.variant_SNPs[barcode] = []
            self.recursive_variant_SNPs(self.variants[barcode], barcode)

    
    def recursive_variant_SNPs(self, variants, barcode):
        for variant in variants:
            vmut = variant["mutations"]
            for snp in vmut:
                self.variant_SNPs[barcode].append(snp)
            self.recursive_variant_SNPs(variant["subs"], barcode)

    def read_variants(self):
        self.variants = {}
        for barcode in self.sorted_barcodes:
            self.variants[barcode] = self.read_variants_for_barcode(barcode)

    def read_variants_for_barcode(self, barcode):
        with open(self.summary_path+"variants/"+barcode+".json", "r") as read_file:
            data = json.load(read_file)
        return data

#######   get_SNPs (unused)  ################################################################
    """
    def get_SNPs(self,barcode, force_reload=False):
        if not self.SNPs:
            self.read_SNPs()
        return self.SNPs[barcode]
    
    def read_SNPs(self):
        self.SNPs = {}
        for barcode in self.sorted_barcodes:
            self.variants[barcode] = self.read_SNPs_for_barcode(barcode)

    def read_SNPs_for_barcode(self, barcode):
        with open(self.summary_path+"SNPs/"+barcode+".csv", 'r') as f:
            self.SNPs[barcode]=[]
            f.readline()
            for line in f:
                l = line.split(",")
                self.SNPs[barcode].append([l[0], int(l[1]), int(l[2])])
    """

####################################################################################
###### AMPLICONS #####
####################################################################################

    def get_amplicons(self):
        if len(self.amplicons)==0:
            return self.read_amplicons()
        return self.amplicons

    def read_amplicons(self):
        f = open(self.ampliconsjson)
        data = json.load(f)
        int_start_end = [ [int(a[0]), int(a[1])] for a in data["amplicons"] ]
        self.amplicons = int_start_end
        return self.amplicons

    def int_before_dash(amp):
        amp_split = amp.split("-")
        return int(amp_split[0])

    def create_amplicon_data(self,trim=500, ignore_unclassified=True):
        f = open(self.ampliconsjson)
        data = json.load(f)
        amplicon_captions = [] 
        amplicons = sorted(data['amplicons'])
        for i in range(len(amplicons)):
            amplicon_captions.append(str(amplicons[i][0])+"-"+str(amplicons[i][1]))
        if not ignore_unclassified:
            barcode_captions = self.sorted_barcodes
        else:
            barcode_captions = [b for b in self.sorted_barcodes if (not b=="unclassified" or len(self.sorted_barcodes)==1)]
        barcode_amplicon = np.array([self.avgs_for_amplicons(amplicons, self.pos_coverage[barcode],trim, barcode) for barcode in barcode_captions])
        f.close()
        return amplicon_captions, barcode_captions, barcode_amplicon

    def avgs_for_amplicons(self,amplicons, barcode, trim, barcode_name):
        avg_cov = []
        for amplicon in amplicons:
            amp_arr = [v[1] for v in barcode[int(amplicon[0])-1+trim : amplicon[1]-1-trim]]
            avg_cov.append(sum(amp_arr)/len(amp_arr))
        return avg_cov

    def divide_amplicons_into_pools(self):
        if len(self.amplicons) == 0:
            self.read_amplicons()
        if len(self.amplicon_pools) == 0:
            pool = 0
            prev = [-1,-1]
            pools = [[]]
            for amplicon in self.amplicons:
                if prev[1]<amplicon[0]:
                    pools[pool].append(amplicon)
                else:
                    pools.append([])
                    pool+=1
                    pools[pool].append(amplicon)
                prev = amplicon
            self.amplicon_pools = pools
            return pools
        else:
            return self.amplicon_pools
