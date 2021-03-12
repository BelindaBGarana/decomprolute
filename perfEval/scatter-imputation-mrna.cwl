#!/usr/bin/env cwltool
class: Workflow
label: scatter-imputation-mrna
id: scatter-imputation-mrna
cwlVersion: v1.2


requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement

   
inputs:
  tissueTypes:
    type: string[]
  cancerTypes:
    type: string[]
  algorithms:
    type: string[]
  signatures:
    type: File[]
      
outputs:
  raw-data:
    type: File[]
    outputSource: download-mrna/matrix
  imputed-data:
    type: File[]
    outputSource: impute-mrna/matrix
  # deconvoluted:
  #   type: File[]
  #   outputSource: deconvolution/deconvoluted
  imputed-deconvoluted:
    type: File[]
    outputSource: imputed-deconvolution/deconvoluted

steps:
  download-raw:
    run: ../protData/mrna-data-cwl-tool.cwl
    scatter: [cancerType, sampleType]
    scatterMethod: flat_crossproduct
    in:
      cancerType: cancerTypes
      sampleType: tissueTypes
    out:
      [matrix]

  impute-prot:
    run: ../imputation/imputation-tool.cwl
    scatter: [input_f]
    scatterMethod: flat_crossproduct
    in:
      input_f:
        source: download-prot/matrix
      use_missForest:
        valueFrom: 'false'
    out:
      [matrix]

  deconv-prot:
    run: prot-deconv-imputed.cwl
    scatter: [expression, protAlg, signature]
    scatterMethod: flat_crossproduct
    in:
      expression: 
        source: download-prot/matrix
      protAlg: prot-algorithms
      signature: signatures
    out: [deconvoluted]

  deconv-prot-imputed:
    run: prot-deconv-imputed.cwl
    scatter: [expression, protAlg, signature]
    scatterMethod: flat_crossproduct
    in:
      expression:
        source: impute-prot/matrix
      protAlg: prot-algorithms
      signature: signatures
    out: [deconvoluted]