# 仓库维护

## 文件放在哪里

- 首页只放项目目的、入口与常用操作。
- `docs/` 放跨项目说明、学习路线与评估；科学脚本、报告、数据留在各实验原目录。
- `tools/` 放仓库维护与运行准备工具，`tests/` 检验这些工具。
- `runs/` 放新运行副本，已加入忽略规则。结果审查后再选择必要的输入、报告、日志和检查点纳入对应实验目录。
- `_archive/` 保留历史来源。不要仅因同名或文件哈希相同就删除研究证据。

## 更新台账

台账工具需 NumPy，在仓库根目录执行：

```bash
python -m pip install numpy
# 先暂存本次要提交的文件，台账哈希按 Git 暂存区规范化内容计算
# 用具体文件路径替换 <files>
git add -- <files>
python tools/artifact_ledger.py
git add ARTIFACT_LEDGER.md ARTIFACT_LEDGER.csv
git diff --cached --stat
git commit
```

`ARTIFACT_LEDGER.md/csv` 是派生文件，不手工编辑。已跟踪文件的大小和 SHA256 来自暂存区 blob；内容分析读取工作区，所以生成前要暂存完整的文件，避免混用两个版本。台账会覆盖所有未忽略的项目、文档与工具目录，新增实验不再需要手工注册才能被收录。

Windows 可用 `refresh_ledger.cmd`。可选钩子安装方式见 [tools/README](../tools/README.md)。该钩子会自动暂存台账，生成失败时仍放行；重要提交应手工检查台账。

## 缺失的大文件

全新克隆不含轨迹、SOAP 矩阵、文献 PDF 等被忽略文件。恢复清单见 [DATA_MANIFEST.md](../DATA_MANIFEST.md)，完整哈希见同名 CSV。

只在保存了完整研究数据的原机器上执行 `python tools/artifact_ledger.py --data-manifest`。不要在缺数据的新克隆里重建该清单，以免把可恢复项误写成当前机器的缓存或残缺数据。

`.restart` 保存续算状态，已纳入 Git；不能用只含部分状态的结构文件随意替代。

## 提交前检查

```bash
python -m unittest discover -s tests -v
git diff --check
git diff --cached --check
```

工具测试检查 SOAP 副本准备、目录完整覆盖及台账稳定性，不代表 MD 或 SOAP 科学结果已复现。科学代码有独立测试和运行环境要求，见[工作流](WORKFLOWS.md)。

`.gitattributes` 规定普通文本和 shell 脚本使用 LF，Windows 批处理与 PowerShell 使用 CRLF，结构文件和二进制数据保持原始字节。
