import pandas as pd
import numpy as np

def merge_dexseq_rmats(dexseqResults, MATS, fromGTF):
    event_Dexseq = fromGTF[['ID', 'DexseqFragment']]
    DXFrag_MATS = MATS.merge(event_Dexseq, how="left", on="ID")

    DXFrag_MATS['DexseqFragment'] = DXFrag_MATS["DexseqFragment"].str.split(',')
    DXFrag_MATS = DXFrag_MATS.explode('DexseqFragment')

    dexseqResults_event = dexseqResults.merge(DXFrag_MATS, how="inner", left_on=['groupID', 'featureID'], right_on=['GeneID', 'DexseqFragment'])

    return dexseqResults_event


dexseqResults = pd.read_table("/home/dwito/Raji_Jurkat_Cells/dexseq/Raji_Jurkat_dexseq_corrected_padj.txt")


A3SS_MATS = pd.read_table("/home/dwito/Raji_Jurkat_Cells/rmats/A3SS.MATS.JCEC.txt")
fromGTF_A3SS = pd.read_table("/home/dwito/merging_rmats_dexseq/genome_wide_analysis_raji_jurkat/combined_fromGTFs/combined_fromGTF.A3SS.txt")
A3SS_rmats_dexseq = merge_dexseq_rmats(dexseqResults, A3SS_MATS, fromGTF_A3SS)
A3SS_rmats_dexseq['Event'] = 'A3SS'


A5SS_MATS = pd.read_table("/home/dwito/Raji_Jurkat_Cells/rmats/A5SS.MATS.JCEC.txt")
fromGTF_A5SS = pd.read_table("/home/dwito/merging_rmats_dexseq/genome_wide_analysis_raji_jurkat/combined_fromGTFs/combined_fromGTF.A5SS.txt")
A5SS_rmats_dexseq = merge_dexseq_rmats(dexseqResults, A5SS_MATS, fromGTF_A5SS)
A5SS_rmats_dexseq['Event'] = 'A5SS'


SE_MATS = pd.read_table("/home/dwito/Raji_Jurkat_Cells/rmats/SE.MATS.JCEC.txt")
fromGTF_SE = pd.read_table("/home/dwito/merging_rmats_dexseq/genome_wide_analysis_raji_jurkat/combined_fromGTFs/combined_fromGTF.SE.txt")
SE_rmats_dexseq = merge_dexseq_rmats(dexseqResults, SE_MATS, fromGTF_SE)
SE_rmats_dexseq['Event'] = 'SE'


RI_MATS = pd.read_table("/home/dwito/Raji_Jurkat_Cells/rmats/RI.MATS.JCEC.txt")
fromGTF_RI = pd.read_table("/home/dwito/merging_rmats_dexseq/genome_wide_analysis_raji_jurkat/combined_fromGTFs/combined_fromGTF.RI.txt")
RI_rmats_dexseq = merge_dexseq_rmats(dexseqResults, RI_MATS, fromGTF_RI)
RI_rmats_dexseq['Event'] = 'RI'


frames = [A3SS_rmats_dexseq, A5SS_rmats_dexseq, SE_rmats_dexseq, RI_rmats_dexseq]
rmats_dexseq_all = pd.concat(frames)
rmats_dexseq_all.to_csv("/home/dwito/Raji_Jurkat_Cells/rmats_dexseq_all.txt", index=False, sep='\t')