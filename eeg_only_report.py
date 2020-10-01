#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:02:58 2019

@author: miguelrebelo
"""
#This script takes ~60min to run

#CAUTION! change MAF file, Linear, WALD, and VCF
#Important files: Annotaton file, DisGeNET variants file, Targets file, and AD meta-analysis variants file

import sys
old_stdout = sys.stdout

log_file = open("frequency_bands_report.log","w")

sys.stdout = log_file


print('start report')

import pandas as pd
#read linear regression files
p1 = pd.read_csv('linear_add_eaad_all_05.alpha.assoc.linear', delim_whitespace=True)
p2 = pd.read_csv('linear_add_eaad_all_05.beta1.assoc.linear', delim_whitespace=True)
p3 = pd.read_csv('linear_add_eaad_all_05.beta2.assoc.linear', delim_whitespace=True)
p4 = pd.read_csv('linear_add_eaad_all_05.delta.assoc.linear', delim_whitespace=True)
p6 = pd.read_csv('linear_add_eaad_all_05.theta.assoc.linear', delim_whitespace=True)
#get beta1 coefficients
alpha = p1[p1.TEST=='ADD']
beta1 = p2[p2.TEST=='ADD']
beta2 = p3[p3.TEST=='ADD']
delta = p4[p4.TEST=='ADD']
theta = p6[p6.TEST=='ADD']
#write to new file
alpha.to_csv('alpha.txt', sep='\t', index=False)
beta1.to_csv('beta1.txt', sep='\t', index=False)
beta2.to_csv('beta2.txt', sep='\t', index=False)
delta.to_csv('delta.txt', sep='\t', index=False)
theta.to_csv('theta.txt', sep='\t', index=False)
print('linear regression files written as "-band-.txt"')
#rename P column for future reference
alpha.rename({'P':'P_alpha'}, axis=1, inplace=True)
beta1.rename({'P':'P_beta1'}, axis=1, inplace=True)
beta2.rename({'P':'P_beta2'}, axis=1, inplace=True)
delta.rename({'P':'P_delta'}, axis=1, inplace=True)
theta.rename({'P':'P_theta'}, axis=1, inplace=True)
#subset p<0.05
subset_t = 0.01
snps = alpha[alpha.P_alpha < subset_t].append(beta1[beta1.P_beta1 < subset_t]).append(beta2[beta2.P_beta2 < subset_t]).append(delta[delta.P_delta < subset_t]).append(theta[theta.P_theta < subset_t])
#print how many variants were selected
print(len(snps.drop_duplicates(['SNP'])), ' variants selected at p < 0.01')
snps.SNP.drop_duplicates().to_csv('snps_01_cut.txt', index=False, header=None)
print('SNP subset at P<0.01 exported as "snps_01_cut.txt"')
#annot file
annot = pd.read_csv('multianno_AD.hg19_multianno.txt', sep='\t')
annot.rename({'Otherinfo':'SNP', 'Func.wgEncodeGencodeBasicV19': 'consequence', 'Gene.wgEncodeGencodeBasicV19': 'gene', 'CLNSIG': 'clinsig', 'GeneDetail.wgEncodeGencodeBasicV19': 'detail', 'ExonicFunc.wgEncodeGencodeBasicV19': 'exonic_func', 'AAChange.wgEncodeGencodeBasicV19': 'aa_change'}, axis=1, inplace=True)
annot = annot.filter(['SNP', 'start', 'end', 'Ref', 'Alt', 'consequence', 'gene', 'detail', 'exonic_func', 'aa_change', 'CADD13_PHRED', 'clinsig'])
#How many variants on disgenet database (drop duplicates!).
disgenet = pd.read_table('C0002395_disease_vda_summary.tsv')
disgenet.rename({'Variant':'SNP'}, axis='columns', inplace=True)
contain = snps.merge(disgenet, on='SNP', how='left')
eeg_disgnet_05 = contain[contain.Gene.notna()]
print(len(eeg_disgnet_05.drop_duplicates(['SNP'])), 'variants from subset in disgenet')
print('from ', len(disgenet), ' variants present in disgenet dataset')

#How many variants have log regression data (drop duplicates!) and check significance for logistic regression.
log = pd.read_csv('log_add_eaa_all_05.assoc.logistic', delim_whitespace=True)
log_in_eeg_05 = log[log.TEST == 'ADD'][log.SNP.isin(contain.SNP[contain.Gene.notna()])].drop_duplicates(['SNP'])
print('log regression info for ', len(log_in_eeg_05), ' variants')
log_in_eeg_05.to_csv('eeg_log.txt', sep='\t', index=False)
print('p_values from logistic regression written to "eeg_log.txt"')
print(len(log_in_eeg_05[log_in_eeg_05.P < 0.05]), ' variants with p_log < 0.05')
print(log_in_eeg_05[log_in_eeg_05.P < 0.05])
print(len(log_in_eeg_05[log_in_eeg_05.P < 0.001]), ' variants with p_log < 0.001')
print(log_in_eeg_05[log_in_eeg_05.P < 0.001])

#Check if there’s some variant (or proxy) previously associated to AD by meta-analysis (European Alzheimer’s Disease Initiative (EADI) et al., 2013; Jansen et al., 2019).
meta = pd.read_table('AD-meta.txt')
print(meta[meta.SNP.isin(eeg_disgnet_05.SNP)])
print(len(meta[meta.SNP.isin(eeg_disgnet_05.SNP)]), ' snps previously associated with AD from meta-analysis')

#Incorporate MAF information.
freq = pd.read_csv('Common_SNPs_eaa_all_05.frq.strat', delim_whitespace=True)
freq.drop(['A1','A2', 'NCHROBS', 'CHR', 'MAC'], axis=1, inplace=True)
freq2 = freq.pivot_table(index=['SNP'], columns='CLST', aggfunc= lambda x: x)
freq2.columns = freq2.columns.droplevel().rename(None)
freq2 = freq2.reset_index()
freq2.rename({'case':'MAF_cases', 'control':'MAF_controls'}, axis=1, inplace=True)

eeg_disgenet_freq = eeg_disgnet_05.merge(freq2, on='SNP', how='left')
eeg_disgenet_freq = eeg_disgenet_freq.merge(annot, on='SNP', how='left')
eeg_disgenet_freq.to_csv('snp_in_disgenet_freq.txt', sep='\t', index=False)
print('variants present in disgenet with MAF info for cases and controls written to "snp_in_disgenet_freq.txt"')

#Most significant from logistic regression among the nominally significant for linear regression
log_add = log[log.TEST == 'ADD']
log_add.to_csv('log.txt', sep='\t', index=False) #for manhattan
snps2 = snps.drop_duplicates(['SNP'])
login=log_add[log_add.SNP.isin(snps2.SNP)]
print(len(login[login.P < 0.05]), ' variants with log P < 0.05')
print(len(login[login.P < 0.001]), ' variants with log P < 0.001')
print(len(login[login.P < 0.0001]), ' variants with log P < 0.0001')
snps_log=snps.merge(login, on='SNP', how='left').merge(freq2, on='SNP', how='left').merge(annot, on='SNP', how='left')
snps_log = snps_log.filter(['SNP', 'A1', 'P', 'OR', 'P_alpha', 'P_beta1', 'P_beta2', 'P_delta', 'P_theta', 'BETA', 'SE', 'gene', 'MAF_cases', 'MAF_controls'])
snps_log_sig = snps_log[snps_log.P < 0.0001]
snps_log_sig.to_csv('sig_log_in_lin_10-4.txt', sep='\t', index=False)
print('log significant snps at 10e-4 exported as "sig_log_in_lin_10-4.txt"')


#add OR and logistic P values
summary = eeg_disgenet_freq.merge(log_in_eeg_05, on='SNP', how='left', suffixes=('','_'))
summary = summary[['SNP', 'A1', 'CHR', 'BP', 'P_alpha', 'P_beta1', 'P_beta2', 'P_delta', 'P_theta', 'BETA', 'SE', 'Gene', 'Consequence', 'Alt_Ref_Genome', 'OR', 'P', 'MAF_cases', 'MAF_controls']]

import numpy as np
summary['effect'] = pd.np.where(summary.OR > 1, 'RISK', 'PROTECT')
summary['duplicated'] = summary.SNP.duplicated()


summary = summary.merge(annot, on='SNP', how='left')
summary.to_csv('summary.txt', index=False, sep='\t')
print('summary info with linear and logistic p_values as well as MAF and duplicates info written to "summary.txt" file')

print('getting genotypes and frequency bands')

# Get genotypes from the VCF file. subset VCF exported for performance reasons.
import io

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

vcf = read_vcf('Common_SNPs_eaa_all_05.vcf')
vcf.drop(['CHROM', 'POS', 'FILTER', 'REF', 'ALT', 'QUAL', 'INFO', 'FORMAT'], 
         axis='columns', inplace=True)
vcf2 = vcf.set_index('ID').T
vcf2.reset_index(level=0, inplace=True)
vcf2.rename({'index':'IID'}, axis='columns', inplace=True)
vcf2.IID = vcf2.IID.str.replace('\_ALZ\_\d*\.CEL','')
df=vcf2.replace('0/1','1').replace('1/0', '1').replace('0/0',
                 '0').replace('1/1', '2').replace('./.', '')
#Read Frequency band info
eeg = pd.read_csv('/Users/miguelrebelo/Desktop/IPATIMUP/AD-EEGWA/Entregable 2.5 AnexoI.txt', sep='\t')
eeg = pd.DataFrame.dropna(eeg)
eeg.rename({'Unnamed: 0': 'id'}, axis='columns', inplace=True)
#Read sample codes to identify samples
samples = pd.read_table('SAMPLE_codes.txt', header=None)
samples.rename({0:'IID', 1:'id'}, axis='columns', inplace=True)
samples.id = samples.id.str.replace('\-','_')
#identify eeg samples
eeg2 = samples.merge(eeg, on='id', how='left')
eeg2 = eeg2[['IID', 'RP-alpha', 'RP_beta1', 'RP_beta2', 'RP_delta', 'RP_theta']]
#merge with genotypes
df2 = df.merge(eeg2, how='left', on='IID')

#automaticly subsets selected variants
lista = ['IID', 'RP-alpha', 'RP_beta1', 'RP_beta2', 'RP_delta', 'RP_theta'] + list(eeg_disgnet_05.SNP.drop_duplicates())
df3 = df2[df2.columns.intersection(lista)]
#keeps selected variants as individual columns and melts waves
lista2 = ['IID'] + list(eeg_disgnet_05.SNP.drop_duplicates())
df4 = df3.melt(id_vars=lista2)
df4.rename({'variable':'Frequency band', 'value':'Relative Power'}, axis='columns', inplace=True)
#melt snps to keep genotypes as values
df5 = df4.melt(id_vars=['IID', 'Frequency band', 'Relative Power'])
df5['Frequency band'] = df5['Frequency band'].str.replace('^RP\_', '').str.replace('^RP\-', '')
df5.rename({'variable':'variant', 'value':'n altered alleles'}, axis='columns', inplace=True)

#exclude wave if needed
#df6 = df5[df5['Frequency band'] != 'gamma']
df6 = df5.copy()
df6['presence of altered alleles'] = df6['n altered alleles'].replace('2','1')

#Strat ploting
print('plots')
import seaborn as sns
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

plt.rc('font', **font)

import pathlib
import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
#list selected SNPs
lista_rs = list(eeg_disgnet_05.SNP.drop_duplicates())
#create folder to dump plots
####pathlib.Path('genotype/').mkdir(parents=True, exist_ok=True)
#plot with genotype information
####with cd('genotype/'):
####    for rs in lista_rs:
####        sns.catplot(data=df6[df6['n altered alleles'] != ''][df6.variant == rs], x='Frequency band', y='Relative Power', hue='n altered alleles', palette='pastel', ci='sd', errcolor='grey', kind='bar', errwidth=0.5)
####        plt.subplots_adjust(top=0.88)
####        plt.title(rs)
####        plt.savefig(rs)
####print('RP vs genotype dumped in genotype folder')

####pathlib.Path('allele/').mkdir(parents=True, exist_ok=True)
#plot with allele info
####with cd('allele/'):
####    for rs in lista_rs:
####        sns.catplot(data=df6[df6['presence of altered alleles'] != ''][df6.variant == rs], x='Frequency band', y='Relative Power', hue='presence of altered alleles', palette='pastel', ci='sd', errcolor='grey', kind='bar', errwidth=0.5)
####        plt.subplots_adjust(top=0.88)
####        plt.title(rs)
####        plt.savefig(rs)
####print('RP vs allele dumped in allele folder')

####print('ploting RP vs genotype distributions for all samples and cases and controls separatelly')

#import case/control info
fam = pd.read_csv('Common_SNPs_PLINK.fam', delim_whitespace=True, header=None)
fam.drop([1,2,3], axis='columns', inplace=True)
fam.rename({0:'IID', 4:'sex', 5:'phenotype'}, axis='columns', inplace=True)
df7 = df6.merge(fam, on='IID', how='left')
df7.phenotype = df7.phenotype.replace(2,'case').replace(1, 'control')
df7_cases = df7[df7.phenotype == 'case']
df7_controls = df7[df7.phenotype == 'control']
#plot ditributions for cases, controls and all

####pathlib.Path('all_dist/').mkdir(parents=True, exist_ok=True)
####with cd('all_dist/'):
####    for rs in lista_rs:
####        sns.catplot(data=df7[df7['n altered alleles'] != ''][df7.variant == rs], x='n altered alleles', y='Relative Power', hue='phenotype', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype ditribution for cases and controls dumped in all_dist folder')

####pathlib.Path('cases_dist/').mkdir(parents=True, exist_ok=True)
####with cd('cases_dist/'):
####    for rs in lista_rs:
####        sns.catplot(data=df7_cases[df7_cases['n altered alleles'] != ''][df7_cases.variant == rs], x='n altered alleles', y='Relative Power', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype distribution for cases dumped in cases_dist folder')

####pathlib.Path('controls_dist/').mkdir(parents=True, exist_ok=True)
####with cd('controls_dist/'):
####    for rs in lista_rs:
####        sns.catplot(data=df7_controls[df7_controls['n altered alleles'] != ''][df7_controls.variant == rs], x='n altered alleles', y='Relative Power', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype distribution for controls dumped in controls_dist folder')


#Start analysis for most significant snps from linear regression
print('start most significant')
tsld = 0.05/len(log[log.TEST == 'ADD'])
print('genome-wide threshold: ', tsld)
print('number of variants being tested: ', len(log[log.TEST == 'ADD']))

snps = alpha[alpha.P_alpha < tsld].append(beta1[beta1.P_beta1 < tsld]).append(beta2[beta2.P_beta2 < tsld]).append(delta[delta.P_delta < tsld]).append(theta[theta.P_theta < tsld])
sum2 = snps.merge(log[log.TEST == 'ADD'], on='SNP', how='left', suffixes=('','_')).merge(freq2, on='SNP', how='left', suffixes=('','_'))
sum2 = sum2[['SNP', 'A1', 'CHR', 'BP', 'P_alpha', 'P_beta1', 'P_beta2', 'P_delta', 'P_theta', 'BETA', 'SE', 'OR', 'P', 'MAF_cases', 'MAF_controls']]
sum2['effect'] = pd.np.where(sum2.OR > 1, 'RISK', 'PROTECT')
sum2['duplicated'] = sum2.SNP.duplicated()
sum2 = sum2.merge(annot, on='SNP', how='left')
sum2.to_csv('summary_sig.txt', index=False, sep='\t')
print('summary info with linear and logistic p_values as well as MAF info written to "summary_sig.txt" file')
print('found ', len(sum2.drop_duplicates(['SNP'])),' highly significant SNPs:')
print(sum2)

#list selected variants and subset
lista3 = ['IID', 'RP-alpha', 'RP_beta1', 'RP_beta2', 'RP_delta', 'RP_theta'] + list(sum2.SNP.drop_duplicates())
df8 = df2[df2.columns.intersection(lista3)]

#melt waves
lista4 = ['IID'] + list(sum2.SNP.drop_duplicates())
df9 = df8.melt(id_vars=lista4)
#melt snps
df9.rename({'variable':'Frequency band', 'value':'Relative Power'}, axis='columns', inplace=True)
df10 = df9.melt(id_vars=['IID', 'Frequency band', 'Relative Power'])
df10['Frequency band'] = df10['Frequency band'].str.replace('^RP\_', '').str.replace('^RP\-', '')
df10.rename({'variable':'variant', 'value':'n altered alleles'}, axis='columns', inplace=True)
#convert genotype to allelic info
df10['presence of altered alleles'] = df10['n altered alleles'].replace('2','1')

#list selected snps
lista_rs2 = list(sum2.SNP.drop_duplicates())

pathlib.Path('genotype_sig/').mkdir(parents=True, exist_ok=True)
with cd('genotype_sig/'):
    for rs in lista_rs2:
        sns.catplot(data=df10[df10['n altered alleles'] != ''][df10.variant == rs], x='Frequency band', y='Relative Power', hue='n altered alleles', palette='pastel', ci='sd', errcolor='grey', kind='bar', errwidth=0.5)
        plt.subplots_adjust(top=0.88)
        plt.title(rs)
        plt.savefig(rs)
print('RP vs genotype for significant SNPs dumped in genotype_sig folder')

pathlib.Path('allele_sig/').mkdir(parents=True, exist_ok=True)

with cd('allele_sig/'):
    for rs in lista_rs2:
        sns.catplot(data=df10[df10['presence of altered alleles'] != ''][df10.variant == rs], x='Frequency band', y='Relative Power', hue='presence of altered alleles', palette='pastel', ci='sd', errcolor='grey', kind='bar', errwidth=0.5)
        plt.subplots_adjust(top=0.88)
        plt.title(rs)
        plt.savefig(rs)
print('RP vs allele for significant SNPs dumped in allele_sig folder')

print('ploting RP vs genotype distributions for all samples and cases and controls separatelly')

df11 = df10.merge(fam, on='IID', how='left')
df11.phenotype = df11.phenotype.replace(2,'case').replace(1, 'control')
df11_cases = df11[df11.phenotype == 'case']
df11_controls = df11[df11.phenotype == 'control']

pathlib.Path('all_dist_sig/').mkdir(parents=True, exist_ok=True)
with cd('all_dist_sig/'):
    for rs in lista_rs2:
        sns.catplot(data=df11[df11['n altered alleles'] != ''][df11.variant == rs], x='n altered alleles', y='Relative Power', hue='phenotype', palette='pastel', col='Frequency band')
        plt.subplots_adjust(top=0.88)
        plt.savefig(rs)
print('RP vs genotype for significant SNPs distribution for all samples dumped in all_dist_sig folder')

####pathlib.Path('cases_dist_sig/').mkdir(parents=True, exist_ok=True)
####with cd('cases_dist_sig/'):
####    for rs in lista_rs2:
####        sns.catplot(data=df11_cases[df11_cases['n altered alleles'] != ''][df11_cases.variant == rs], x='n altered alleles', y='Relative Power', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype for significant SNPs distribution for cases dumped in cases_dist_sig folder')

####pathlib.Path('controls_dist_sig/').mkdir(parents=True, exist_ok=True)
####with cd('controls_dist_sig/'):
####    for rs in lista_rs2:
####        sns.catplot(data=df11_controls[df11_controls['n altered alleles'] != ''][df11_controls.variant == rs], x='n altered alleles', y='Relative Power', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype for significant SNPs distribution for controls dumped in controls_dist_sig folder')


tsld3 = 0.000005

#other significant ones
print('Start analysis for most significant snps from linear regression, with low SE and significant beta')
snps = alpha[alpha.P_alpha < tsld3].append(beta1[beta1.P_beta1 < tsld3]).append(beta2[beta2.P_beta2 < tsld3]).append(delta[delta.P_delta < tsld3]).append(theta[theta.P_theta < tsld3])
#snps = snps[abs(snps.SE) < 0.1]
sum3 = snps.merge(log[log.TEST == 'ADD'], on='SNP', how='left', suffixes=('','_')).merge(freq2, on='SNP', how='left', suffixes=('','_'))
sum3 = sum3[['SNP', 'A1', 'CHR', 'BP', 'P_alpha', 'P_beta1', 'P_beta2', 'P_delta', 'P_theta', 'BETA', 'SE', 'OR', 'P', 'MAF_cases', 'MAF_controls']]
sum3['effect'] = pd.np.where(sum3.OR > 1, 'RISK', 'PROTECT')
sum3['duplicated'] = sum3.SNP.duplicated()
sum3 = sum3.merge(annot, on='SNP', how='left')
sum3.to_csv('summary_sig_se.txt', index=False, sep='\t')
print('summary info with linear p < ', tsld3, ' and logistic p_values as well as MAF info written to "summary_sig_se.txt" file')
print('found ', len(sum3.drop_duplicates(['SNP'])),' significant SNPs at p < ', tsld3)
print(sum3)

#list selected variants and subset
lista3 = ['IID', 'RP-alpha', 'RP_beta1', 'RP_beta2', 'RP_delta', 'RP_theta'] + list(sum3.SNP.drop_duplicates())
df8 = df2[df2.columns.intersection(lista3)]

#melt waves
lista4 = ['IID'] + list(sum3.SNP.drop_duplicates())
df9 = df8.melt(id_vars=lista4)
#melt snps
df9.rename({'variable':'Frequency band', 'value':'Relative Power'}, axis='columns', inplace=True)
df10 = df9.melt(id_vars=['IID', 'Frequency band', 'Relative Power'])
df10['Frequency band'] = df10['Frequency band'].str.replace('^RP\_', '').str.replace('^RP\-', '')
df10.rename({'variable':'variant', 'value':'n altered alleles'}, axis='columns', inplace=True)
#convert genotype to allelic info
df10['presence of altered alleles'] = df10['n altered alleles'].replace('2','1')

#list selected snps
lista_rs3 = list(sum3.SNP.drop_duplicates())

pathlib.Path('genotype_sig_se/').mkdir(parents=True, exist_ok=True)
with cd('genotype_sig_se/'):
    for rs in lista_rs3:
        sns.catplot(data=df10[df10['n altered alleles'] != ''][df10.variant == rs], x='Frequency band', y='Relative Power', hue='n altered alleles', palette='pastel', ci='sd', errcolor='grey', kind='bar', errwidth=0.5)
        plt.subplots_adjust(top=0.88)
        plt.title(rs)
        plt.savefig(rs)
print('RP vs genotype for significant SNPs filtered by SE dumped in genotype_sig_se folder')

pathlib.Path('allele_sig_se/').mkdir(parents=True, exist_ok=True)

with cd('allele_sig_se/'):
    for rs in lista_rs3:
        sns.catplot(data=df10[df10['presence of altered alleles'] != ''][df10.variant == rs], x='Frequency band', y='Relative Power', hue='presence of altered alleles', palette='pastel', ci='sd', errcolor='grey', kind='bar', errwidth=0.5)
        plt.subplots_adjust(top=0.88)
        plt.title(rs)
        plt.savefig(rs)
print('RP vs allele for significant SNPs filterd by SE dumped in allele_sig_se folder')

print('ploting RP vs genotype distributions for all samples and cases and controls separatelly')

df11 = df10.merge(fam, on='IID', how='left')
df11.phenotype = df11.phenotype.replace(2,'case').replace(1, 'control')
df11_cases = df11[df11.phenotype == 'case']
df11_controls = df11[df11.phenotype == 'control']

pathlib.Path('all_dist_sig_se/').mkdir(parents=True, exist_ok=True)
with cd('all_dist_sig_se/'):
    for rs in lista_rs3:
        sns.catplot(data=df11[df11['n altered alleles'] != ''][df11.variant == rs], x='n altered alleles', y='Relative Power', hue='phenotype', palette='pastel', col='Frequency band')
        plt.subplots_adjust(top=0.88)
        plt.savefig(rs)
print('RP vs genotype for significant SNPs filtered by SE distribution for all samples dumped in all_dist_sig_se folder')

####pathlib.Path('cases_dist_sig_se/').mkdir(parents=True, exist_ok=True)
####with cd('cases_dist_sig_se/'):
####    for rs in lista_rs3:
####        sns.catplot(data=df11_cases[df11_cases['n altered alleles'] != ''][df11_cases.variant == rs], x='n altered alleles', y='Relative Power', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype for significant SNPs filtered by SE distribution for cases dumped in cases_dist_sig_se folder')

####pathlib.Path('controls_dist_sig_se/').mkdir(parents=True, exist_ok=True)
####with cd('controls_dist_sig_se/'):
####    for rs in lista_rs3:
####        sns.catplot(data=df11_controls[df11_controls['n altered alleles'] != ''][df11_controls.variant == rs], x='n altered alleles', y='Relative Power', palette='pastel', col='Frequency band')
####        plt.subplots_adjust(top=0.88)
####        plt.savefig(rs)
####print('RP vs genotype for significant SNPs filtered by SE distribution for controls dumped in controls_dist_sig_se folder')




#Manhattan plots

def FormatData(path, sep = '\t', chromosome = 'CHR', p_value = 'P'):
    data = pd.read_table(path, sep = sep)
    data['-log10(p_value)'] = -np.log10(data[p_value])
    data[chromosome] = data[chromosome].astype('category')
    data['ind'] = range(len(data))
    data_grouped = data.groupby((chromosome))
    return data, data_grouped

def GenerateManhattan(pyhattan_object, export_path = None, significance = 6, colors = ['#E24E42', '#008F95'], refSNP = False):
    data = pyhattan_object[0]
    data_grouped = pyhattan_object[1]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(p_value)', color=colors[num % len(colors)], ax=ax, s= 1)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(data)])
    ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
    ax.set_xlabel('Chromosome')
    plt.axhline(y=significance, color='gray', linestyle='-', linewidth = 0.5)
    plt.xticks(fontsize=8, rotation=60)
    plt.yticks(fontsize=8)

    if refSNP:
        for index, row in data.iterrows():
            if row['-log10(p_value)'] >= significance:
                ax.annotate(str(row[refSNP]), xy = (index, row['-log10(p_value)']), xytext=(index + 7, row['-log10(p_value)']+0.2), size=7, arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))

    if export_path:
        plt.savefig(export_path)

    plt.show()

GenerateManhattan(FormatData('alpha.txt'), significance=-np.log10(tsld3), refSNP='SNP', colors=['red', 'blue'], export_path='alpha_man')
GenerateManhattan(FormatData('beta1.txt'), significance=-np.log10(tsld3), refSNP='SNP', colors=['red', 'blue'], export_path='beta1_man')
GenerateManhattan(FormatData('beta2.txt'), significance=-np.log10(tsld3), refSNP='SNP', colors=['red', 'blue'], export_path='beta2_man')
GenerateManhattan(FormatData('delta.txt'), significance=-np.log10(tsld3), refSNP='SNP', colors=['red', 'blue'], export_path='delta_man')
GenerateManhattan(FormatData('theta.txt'), significance=-np.log10(tsld3), refSNP='SNP', colors=['red', 'blue'], export_path='theta_man')

print('manhattan plots exported to current folder')


#Q-Q plots

import decimal
def drange(x, y, jump):
  while x < y:
    yield float(x)
    x += decimal.Decimal(jump)

x=sorted(-np.log10(list(drange(0,1,1/len(alpha.dropna())))))

alpha = p1[p1.TEST=='ADD']
beta1 = p2[p2.TEST=='ADD']
beta2 = p3[p3.TEST=='ADD']
delta = p4[p4.TEST=='ADD']
theta = p6[p6.TEST=='ADD']

lista_w = [alpha, beta1, beta2, delta, theta]


for w in lista_w:
    y = -np.log10(sorted(w.P.dropna(), key=float, reverse=True))
    name =[x for x in globals() if globals()[x] is w][0]
    scat = plt.scatter(x,y, label=name)
    plt.plot([0, 7], [0, 7], 'orange', 1)
    plt.xlabel('Theoretical quantiles')
    plt.ylabel('Ordered -log(p-values)')
    plt.title('Probability Plot')
plt.legend()
plt.savefig('QQ')

print('qq plot saved to QQ-plot folder')


#OPEN targets: Some of them are from experiments using nucleotide sequencing. Others are drugs from clinical trials.
targets = pd.read_csv('targets_AD.csv')
targets.rename({'target.gene_info.symbol':'gene'}, axis=1, inplace=True)
lista = list(sum3.gene)
lista = [i.split(';') for i in lista]
lista = [item for sublist in lista for item in sublist]
lista = [i.split('-', 1)[0] for i in lista]
print(len(targets[targets.gene.isin(lista)]), ' genes from sig_se in OpenTargets platform')
pd.options.display.max_columns = 15
print(targets[targets.gene.isin(lista)])



#manhattan plot
GenerateManhattan(FormatData('log.txt'), significance=-np.log10(0.0005), refSNP='SNP', colors=['darkgreen', 'brown'], export_path='log_man')

print('job finished')


sys.stdout = old_stdout

log_file.close()