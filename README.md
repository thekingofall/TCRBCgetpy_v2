
# TCRBCgetpy_v2

TCRBCgetpy_v2 是一个用于分析 TCR（T细胞受体）序列的工具。

## 安装

TCRBCgetpy_v2 需要 Python 3.8 环境和一些特定的依赖项。请按照以下步骤进行安装：

1. 创建并激活 conda 环境：

   ```bash
   conda create -n tcr_env python=3.8
   conda activate tcr_env
   ```

2. 安装 BLAST：

   ```bash
   conda install -c bioconda blast
   ```

3. 克隆 TCRBCgetpy_v2 仓库：

   ```bash
   git clone https://github.com/thekingofall/TCRBCgetpy_v2.git
   cd TCRBCgetpy_v2
   ```

4. 安装 Python 依赖项：

   ```bash
   pip install pandas
   ```

   注意：其他依赖项（如 re, os, sys, datetime, argparse, gzip, configparser）通常是 Python 标准库的一部分，不需要单独安装。



## 使用方法

### 单个样本使用方法

对于单个样本的分析，使用以下命令：

```bash
python /path/to/TCRBCgetpy_v2/TCRGetpy/PXTCR01_main.py \
    --FQ1 /path/to/sample_R1.fastq.gz \
    --FQ2 /path/to/sample_R2.fastq.gz \
    --Module fsm \
    --ini /path/to/TCRBCgetpy_v2/TCRGetpy/congfig.ini \
    --FN SampleName
```

#### 参数说明

- `--FQ1`: 第一个 FASTQ 文件（R1）的路径
- `--FQ2`: 第二个 FASTQ 文件（R2）的路径
- `--Module`: 使用的模块，这里是 `fsm`
- `--ini`: 配置文件的路径
- `--FN`: 输出文件名前缀

### 批量处理多个样本

对于批量处理多个样本，可以使用以下脚本：

```bash
for i in *; do
    if [ -d "$i" ]; then
        echo "Processing sample: $i"
        python /path/to/TCRBCgetpy_v2/TCRGetpy/PXTCR01_main.py \
            --FQ1 $i/*R1.fastq.gz \
            --FQ2 $i/*R2.fastq.gz \
            --Module fsm \
            --ini /path/to/TCRBCgetpy_v2/TCRGetpy/congfig.ini \
            --FN $i
    fi
done
```

#### 批量处理步骤

1. 确保每个样本都有自己的目录。
2. 每个样本目录中应包含 R1 和 R2 的 FASTQ 文件（gzip 压缩格式）。如果没有请改名
3. 修改命令中的路径，使其指向您的 TCRBCgetpy_v2 安装目录。
4. 在包含所有样本目录的父目录中运行上述命令。

注意：请根据您的实际安装路径和文件结构调整命令中的路径。

