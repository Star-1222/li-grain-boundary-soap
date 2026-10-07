# 环境配置指南

> **主要运行平台：Linux**

## 原始环境信息

| 项目 | 值 |
|------|-----|
| Python 版本 | 3.14.4 |
| 创建平台 | Linux (Ubuntu, 路径 /home/star_22/) |
| 虚拟环境路径 | venv/（Linux 格式） |
| 包管理器 | pip 25.1.1 |

> 原始 venv 仍保留在 Linux 机器 `/home/star_22/li-snap-test/venv`，
> 若该路径未删除，可直接复用（见下方"复用已有 venv"）。

---

## Linux 配置（主环境）

### 方式一：复用已有 venv

如果 Linux 上原始 venv 还在，直接激活即可：

```bash
cd ~/li-snap-test          # 或你的原始项目路径
source venv/bin/activate
cd ~/Li_GB_Project         # 进入整合后的项目目录

# 验证关键依赖
python -c "import ase; print('ASE', ase.__version__)"
python -c "import numpy; print('NumPy', numpy.__version__)"
python -c "import scipy; print('SciPy', scipy.__version__)"
```

### 方式二：新建 venv

```bash
# 1. 进入项目目录
cd ~/Li_GB_Project

# 2. 创建虚拟环境（需要 Python 3.12+）
python3 -m venv venv

# 3. 激活
source venv/bin/activate

# 4. 安装依赖
pip install -r requirements.txt

# 5. 验证
python -c "import ase; print('ASE', ase.__version__)"
python -c "import numpy; print('NumPy', numpy.__version__)"
python -c "import scipy; print('SciPy', scipy.__version__)"
```

### LAMMPS + OpenKIM

LAMMPS 和 OpenKIM 是系统级安装，不在 venv 内。确认已安装：

```bash
lmp -h                          # 检查 LAMMPS
kim-api-collections-manager list  # 检查 OpenKIM 模型
# 需要的模型：SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

---

## Windows 配置（备用）

仅在需要在 Windows 上测试 Python 部分时使用：

```powershell
cd C:\Users\1\Documents\Li_GB_Project
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

> 注意：LAMMPS/OpenKIM 在 Windows 上配置较复杂，建议仅用 Windows 做
> Python 开发与调试，实际 MD 模拟在 Linux 上运行。

---

## 依赖说明

| 包 | 版本 | 用途 |
|----|------|------|
| **ase** | 3.29.0 | 双晶模型生成（Li_Bicrystal_Generator 核心） |
| **numpy** | 2.5.1 | 矩阵运算、数组操作 |
| **scipy** | 1.18.0 | 结构优化、空间分析 |
| **matplotlib** | 3.11.1 | 数据可视化、绘图 |

其余 7 个包为上述包的间接依赖，pip 会自动安装。

---

## 注意事项

1. **Python 版本**：原始环境为 3.14.4（较新），建议使用 3.12+ 以确保兼容性
2. **LAMMPS 不在 venv 内**：LAMMPS + OpenKIM 是独立的系统级安装，不在 Python 依赖范围内
3. **SOAP/DScribe**：如后续需要 SOAP 描述符计算，需额外安装 `dscribe`：
   ```bash
   pip install dscribe
   ```
