# Copyright 2019 by Kent Kawashima.  All rights reserved.

import os
import numpy as np
from .utils import fasta_to_seqobj_list
from .mafft import mafft_align_codons
from .conspos import encode_codon_aln_pos, add_conspos_marker
from .mafft import mafft_ginsi_align, mafft_linsi_align, mafft_einsi_align

def codon_mafft_conspos(
    nucl_path, aln_path, use_aln='einsi', aln_iterations=1000,
    consistent_marker='C', inconsistent_marker='N', 
    consistent_threshold=1.0, # every body has to have template pattern -> 1. Threshold for consistency
    marker_name='conspos_marker', marker_location=0
    ):
    # Perform mafft alignment for global, local, and affine gap methods
    # This is for CDS
    # Paths for codon alignments
    out_suffix='.fna.aln'
    ginsi_cod_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + f'.ginsi.cod{out_suffix}')
    linsi_cod_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + f'.linsi.cod{out_suffix}')
    einsi_cod_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + f'.einsi.cod{out_suffix}')
    # Paths for amino acid alignments
    ginsi_aa_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + '.ginsi.faa.aln')
    linsi_aa_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + '.linsi.faa.aln')
    einsi_aa_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + '.einsi.faa.aln')

    mafft_align_codons(
        nucl_path, ginsi_aa_aln_path, ginsi_cod_aln_path,
        aln_method='ginsi', aln_iterations=aln_iterations)
    mafft_align_codons(
        nucl_path, linsi_aa_aln_path, linsi_cod_aln_path,
        translate=False, aln_method='linsi', aln_iterations=aln_iterations)
    mafft_align_codons(
        nucl_path, einsi_aa_aln_path, einsi_cod_aln_path,
        translate=False, aln_method='einsi', aln_iterations=aln_iterations)

    return add_conspos_marker(
        ginsi_cod_aln_path, linsi_cod_aln_path, einsi_cod_aln_path, aln_path,
        marker_name, marker_location, use_aln,
        consistent_marker, inconsistent_marker, 
        consistent_threshold,
    )
    
#----------------------------------------------------------------------------------------
def intron_mafft_conspos(
    nucl_path, aln_path, aln_iterations=1000,
    marker_name='conspos_marker', marker_location=0, use_aln='einsi',
    consistent_marker='C', inconsistent_marker='N', 
    consistent_threshold=1.0,
    ):
    out_suffix='.fna.aln'
    # Perform mafft alignment for global, local, and affine gap methods
    ginsi_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + f'.ginsi{out_suffix}')
    linsi_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + f'.linsi{out_suffix}')
    einsi_aln_path = os.path.join(
        os.path.dirname(aln_path), 
        os.path.basename(nucl_path) + f'.einsi{out_suffix}')

    mafft_ginsi_align(
        infile = nucl_path, output_aln = ginsi_aln_path,
        iterations=aln_iterations)
    mafft_linsi_align(
        infile = nucl_path, output_aln = linsi_aln_path, 
        iterations=aln_iterations)
    mafft_einsi_align(
        infile = nucl_path, output_aln = einsi_aln_path, 
        iterations=aln_iterations)

    return add_conspos_marker(
        ginsi_aln_path, linsi_aln_path, einsi_aln_path, aln_path,
        marker_name, marker_location, use_aln,
        consistent_marker, inconsistent_marker, 
        consistent_threshold,
    )
