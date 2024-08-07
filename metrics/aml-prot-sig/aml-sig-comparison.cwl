#!/usr/bin/env cwltool
class: Workflow
label: scatter-test
id: scatter-test
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
   prot-algorithms:
      type: string[]
   signatures:
      type: string[]
      
outputs:
   cell-cor-tab:
      type: File
      outputSource: get-celltype-cors/table
   cell-fig:
      type: File[]
      outputSource: get-celltype-cors/fig
   prot-files1:
      type: File[]
      outputSource: run-all-algs-by-sig/prot-file1
   prot-files2:
      type: File[]
      outputSource: run-all-algs-by-sig/prot-file2
    

steps:
   run-all-algs-by-sig:
      run: call-deconv-and-cor.cwl
      scatter: [signature,prot-alg,tissueType,cancerType]
      scatterMethod: flat_crossproduct
      in:
        signature: signatures
        prot-alg: prot-algorithms
        cancerType: cancerTypes
        tissueType: tissueTypes
      out:
        [sig-cor-file,prot-file1,prot-file2]
   get-celltype-cors:
      run: ../figures/plot-figs.cwl
      in:
        metricType:
            valueFrom: "cellType"
        files:
            source: run-all-algs-by-sig/sig-cor-file
      out:
        [table,fig]
