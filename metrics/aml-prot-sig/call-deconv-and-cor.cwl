#!/usr/bin/env cwltool
class: Workflow
label: call-deconv-and-cor
id: call-deconv-and-cor
cwlVersion: v1.2

requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
  
inputs:
   signature1: string
   signature2: string
   prot-alg: string
   cancerType: string
   tissueType: string

outputs:
  sig-cor-file:
     type: File
     outputSource: sig-cor/corr
  prot-file1:
     type: File
     outputSource: deconv-prot1/deconvoluted
  prot-file2:
     type: File
     outputSource: deconv-prot2/deconvoluted

steps:
  deconv-prot1:
     run: ../prot-deconv.cwl
     in:
       cancerType: cancerType
       protAlg: prot-alg
       signature: signature1
       sampleType: tissueType
     out: [deconvoluted]
  deconv-prot2:
     run: ../prot-deconv.cwl
     in:
       cancerType: cancerType
       protAlg: prot-alg
       signature: signature2
       sampleType: tissueType
     out: [deconvoluted]
  sig-cor:
     run: ../correlations/deconv-corr-cwl-tool-sigs.cwl
     in:
       cancerType: cancerType
       protAlg: prot-alg
       signature1: signature1
       signature2: signature2
       sampleType: tissueType
       proteomics1:
         source: deconv-prot1/deconvoluted
       proteomics2:
         source: deconv-prot2/deconvoluted
     out: [corr]

