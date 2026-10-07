#!/bin/bash
# Li_GB_Project 环境激活脚本
# 用法：bash setup_env.sh 或 source setup_env.sh

# 项目路径
PROJECT_DIR="$HOME/Li_GB_Project"
VENV_DIR="$HOME/li-snap-test/venv"

# 进入项目目录
cd "$PROJECT_DIR" || { echo "错误：项目目录不存在 $PROJECT_DIR"; return 1 2>/dev/null || exit 1; }

# 激活 venv
if [ -f "$VENV_DIR/bin/activate" ]; then
    source "$VENV_DIR/bin/activate"
    echo "✓ venv 已激活：$VENV_DIR"
else
    echo "✗ venv 不存在，正在创建..."
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    echo "✓ 新 venv 已创建并安装依赖"
fi

# 验证关键依赖
echo ""
echo "--- 环境验证 ---"
python -c "import ase; print(f'  ASE     {ase.__version__}')" 2>/dev/null || echo "  ASE     未安装"
python -c "import numpy; print(f'  NumPy   {numpy.__version__}')" 2>/dev/null || echo "  NumPy   未安装"
python -c "import scipy; print(f'  SciPy   {scipy.__version__}')" 2>/dev/null || echo "  SciPy   未安装"
python -c "import matplotlib; print(f'  MPL     {matplotlib.__version__}')" 2>/dev/null || echo "  MPL     未安装"

# 检查 LAMMPS
if command -v lmp &>/dev/null; then
    echo "  LAMMPS  $(lmp -h 2>&1 | head -1)"
else
    echo "  LAMMPS  未找到（需系统级安装）"
fi

echo ""
echo "✓ 环境就绪，当前目录：$(pwd)"
echo "  生成双晶模型：python -m Li_Bicrystal_Generator --list"
echo "  退出虚拟环境：deactivate"
