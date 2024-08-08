#!/usr/bin/env cwltool

label: deconv-corr-cwl-tool-sigs
id:  deconv-corr-cwl-tool-sigs
cwlVersion: v1.0
class: CommandLineTool
baseCommand: python

arguments:
  - /bin/correlation-sigs.py

requirements:
  - class: DockerRequirement
    dockerPull: correlations:latest
  - class: StepInputExpressionRequirement
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
  spearOrPears:
    type: string
    inputBinding:
      prefix: --spearOrPears
    default: "spearman"
  cancerType:
    type: string
  protAlg:
    type: string
  signature1:
    type: string
  signature2:
    type: string
  sampleType:
    type: string

outputs:
  corr:
    type: File
    outputBinding:
      glob: "corr.tsv" 
      outputEval: |
        ${
          var name = inputs.sampleType + '-' + inputs.cancerType + '-' + inputs.signature1 + '-to-' + inputs.signature2 + '-corr.tsv'
          self[0].basename = name;
          return self[0]
         }
