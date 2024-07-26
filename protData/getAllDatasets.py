'''
Get all cptac data
downloads data into docker image
'''

import cptac
import os
os.mkdir('/data/')

import pandas as pd
import synapseclient

syn = synapseclient.Synapse()
syn.login()

def getCancerObj(cancertype):
   # cptac.download(dataset=cancertype,source='harmonized',)
    if cancertype == 'brca':
        dat = cptac.Brca()
    elif cancertype == 'ccrcc':
        dat = cptac.Ccrcc()
    elif cancertype == 'coad':
        dat = cptac.Coad()
    elif cancertype == 'gbm':
        dat = cptac.Gbm()
    elif cancertype == 'hnscc':
        dat = cptac.Hnscc()
    elif cancertype == 'lscc':
        dat = cptac.Lscc()
    elif cancertype == 'luad':
        dat = cptac.Luad()
    elif cancertype == 'ovarian':
        dat = cptac.Ov()
    elif cancertype =='pdac':
        dat = cptac.Pdac()
    elif cancertype =='ucec':
        dat = cptac.Ucec()
    elif cancertype =='aml':
        aml_syn = 'syn25714248'
        dat = pd.read_csv(syn.get(aml_syn).path, delimiter='\t')
        # make sure gene names are rownames?
    elif cancertype =='AML_Monocyte_DIA': # AML sorted proteomics
        aml_syn = 'syn58914135'
        dat = pd.read_csv(syn.get(aml_syn).path)
        
        # extract columns with CD14+ samples
        dat = dat.filter(regex="CD14")
    elif cancertype =='AML_Progenitor_DIA': # AML sorted proteomics
        aml_syn = 'syn58914135'
        dat = pd.read_csv(syn.get(aml_syn).path)
        
        # extract columns with CD34+ samples
        dat = dat.filter(regex="CD34")
    elif cancertype =='AML_MSC_DIA': # AML sorted proteomics
        aml_syn = 'syn58914135'
        dat = pd.read_csv(syn.get(aml_syn).path)
        
        # extract columns with MSC samples
        dat = dat.filter(regex="MSC")
    elif cancertype =='AML_Monocyte_TMT': # AML sorted proteomics
        aml_syn = 'syn53493076'
        dat = pd.read_csv(syn.get(aml_syn).path, delimiter='\t')
        
        # load metadata
        meta_syn = 'syn53461075'
        meta = pd.read_excel(
            syn.get(meta_syn).path,
            engine='openpyxl'
        )
        meta['id'] = meta.index + 1
        meta = meta.dropna()
        meta['id2'] = meta['patient'] + '_' + meta['SampleType']

        # extract columns with CD14+ samples
        dat.columns = meta['id2']
        dat = dat.filter(regex="CD14")
    elif cancertype =='AML_Progenitor_TMT': # AML sorted proteomics
        aml_syn = 'syn53493076'
        dat = pd.read_csv(syn.get(aml_syn).path, delimiter='\t')
        
        # load metadata
        meta_syn = 'syn53461075'
        meta = pd.read_excel(
            syn.get(meta_syn).path,
            engine='openpyxl'
        )
        meta['id'] = meta.index + 1
        meta = meta.dropna()
        meta['id2'] = meta['patient'] + '_' + meta['SampleType']

        # extract columns with CD34+ samples
        dat.columns = meta['id2']
        dat = dat.filter(regex="CD34")
    elif cancertype =='AML_MSC_TMT': # AML sorted proteomics
        aml_syn = 'syn53493076'
        dat = pd.read_csv(syn.get(aml_syn).path, delimiter='\t')
        
        # load metadata
        meta_syn = 'syn53461075'
        meta = pd.read_excel(
            syn.get(meta_syn).path,
            engine='openpyxl'
        )
        meta['id'] = meta.index + 1
        meta = meta.dropna()
        meta['id2'] = meta['patient'] + '_' + meta['SampleType']

        # extract columns with MSC samples
        dat.columns = meta['id2']
        dat = dat.filter(regex="MSC")
    else:
        print('Wrong cancer type: '+cancertype)
        exit()
    return dat


for ds in ['brca', 'ccrcc', 'ucec', 'coad','pdac', 'ovarian', 'luad', 'hnscc', 'gbm','lscc']:
    dat=getCancerObj(ds)

    #this call changed in recent version
    dat_list = dat.list_data_sources().set_index('Data type').to_dict()['Available sources']
    clinsource = dat_list['clinical']
    if 'harmonized' in clinsource:
        cs = 'harmonized'
    else:
        cs = clinsource[0]
    dat.get_clinical(cs)
    tsource = dat_list['proteomics']
    res = dat.get_proteomics(tsource[0])
    if res.columns.nlevels == 2:
        res.columns = res.columns.droplevel(1)
    
    print(ds+':',res.shape)
    res.to_csv('/data/'+ds+'.csv')
