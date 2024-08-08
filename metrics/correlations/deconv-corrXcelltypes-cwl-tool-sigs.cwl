#!/usr/bin/env cwltool

label: deconv-corrXcelltypes-cwl-tool-sigs
id:  deconv-corrXcelltypes-cwl-tool-sigs
cwlVersion: v1.0
class: CommandLineTool
baseCommand: python

arguments:
  - /bin/correlationXcelltypes-sigs.py

requirements:
  - class: DockerRequirement
    dockerPull: correlations:latest
  - class: InlineJavascriptRequirement
  
inputs:
  proteomics1:
    type: File
    inputBinding:
      prefix: --proteomics1
  proteomics2:
    type: File
    inputBinding:
      prefix: --proteomics2
  protAlg:
    type: string
  cancerType:
    type: string
  spearOrPears:
    type: string
    inputBinding:
      prefix: --spearOrPears
    default: "spearman"
  signature1:
    type: string
  signature2:
    type: string
  sampleVal:
    type: int
    default: 100
  sampleType:
    type: string
  sampleRep:
    type: int
    default: 0

outputs:
  corr:
    type: File
    outputBinding:
      glob: "corrXcelltypes.tsv" #$(inputs.output)
      outputEval: |
        ${
          var name = inputs.sampleType + '-' + inputs.cancerType + '-' + inputs.signature1 + '-to-' + inputs.signature2 +'-'+ inputs.sampleVal + '-sample-' + inputs.sampleRep +'-cellTypecorr.tsv'
          self[0].basename = name;
          return self[0]
         }
