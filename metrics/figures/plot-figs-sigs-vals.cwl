#!/usr/bin/env cwltool

class: CommandLineTool
label: plot-figs-sigs-vals
id: plot-figs-sigs-vals
cwlVersion: v1.2
baseCommand: Rscript

arguments:
   - /bin/compare_sigs.R

requirements:
   - class: DockerRequirement
     dockerPull: figures:latest
   - class: MultipleInputFeatureRequirement

inputs:
   metricType:
     type: string
     inputBinding:
       position: 1
   metric:
     type: string
     default: "correlation"
     inputBinding:
       position: 2
   files:
     type:
       type: array
       items: File
     inputBinding:
        position: 4
   repNumber:
     type: int
     default: 0
     inputBinding:
        position: 3

outputs:
   table:
     type:
       type: array
       items: File
     outputBinding:
        glob: "*.tsv"
   fig:
     type:
       type: array
       items: File
     outputBinding:
        glob: ["*.pdf","*.html"]
    
