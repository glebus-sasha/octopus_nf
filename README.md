# Your Pipeline Name

This repository contains a Nextflow pipeline for analyzing Next-Generation Sequencing (NGS) data using octopus (https://github.com/luntergroup/octopus).

<!--
  ~ Copyright 2013-2023, Seqera Labs
  ~
  ~ Licensed under the Apache License, Version 2.0 (the "License");
  ~ you may not use this file except in compliance with the License.
  ~ You may obtain a copy of the License at
  ~
  ~     http://www.apache.org/licenses/LICENSE-2.0
  ~
  ~ Unless required by applicable law or agreed to in writing, software
  ~ distributed under the License is distributed on an "AS IS" BASIS,
  ~ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ~ See the License for the specific language governing permissions and
  ~ limitations under the License.
  -->
<html>
<head>
<meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">
</head>
<body>
<pre class="mermaid" style="text-align: center;">
flowchart TB
    subgraph " "
    v0["Channel.fromFilePairs"]
    v1["reference"]
    v5["reference"]
    v7["reference"]
    v9["reference"]
    end
    v2([REFINDEX])
    v3([QCONTROL])
    subgraph " "
    v4[" "]
    v11[" "]
    end
    v6([ALIGN])
    v8([PREPARE])
    v10([VARCALL])
    v0 --> v3
    v1 --> v2
    v2 --> v6
    v3 --> v6
    v3 --> v4
    v5 --> v6
    v6 --> v8
    v6 --> v10
    v7 --> v8
    v8 --> v10
    v9 --> v10
    v10 --> v11

</pre>
<script type="module">
  import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
  mermaid.initialize({ startOnLoad: true });
</script>
</body>
</html>


## Description

The pipeline is implemented in Nextflow and includes several stages for NGS data analysis:

1. **REFINDEX:** Index creation using BWA (Burrows-Wheeler Aligner).
2. **QCONTROL:** Data preprocessing using Fastp.
3. **ALIGN:** Sequence alignment using BWA mem.
4. **PREPARE:** File processing and preparation using Samtools.
5. **VARCALL:** Variant calling using Octopus.

## Usage

### Quick Start

To quickly run the pipeline, use the following command:

```bash
nextflow run glebus-sasha/octopus
```

### Requirements

- Nextflow (https://www.nextflow.io/)
- Singularity (https://sylabs.io/singularity/)

### Running the Pipeline

1. Install all the necessary dependencies such as Nextflow, Singularity.
3. Clone this repository: `git clone https://github.com/glebus-sasha/octopus.git`
4. Navigate to the pipeline directory: `cd octopus`
5. Edit the `nextflow.config` file to set the required parameters, if necessary.
6. Run the pipeline, setting the required parameters, for example:

```bash
nextflow run main.nf
```
