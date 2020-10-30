label: run-mcpcounter-on-mrna-tool
id: run-mcpcounter-on-mrna-tool
cwlversion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/r-mcpcounter:1.1.0--r40_1

arguments:
    - mcpcounter.r

inputs:
    rnaseq:
        type: string
        inputBinding:
            position: 1
    