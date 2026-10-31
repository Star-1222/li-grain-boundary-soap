#!/bin/bash
# Li-SNAP OpenKIM 验证运行脚本
# 使用方法: bash run.sh 或 chmod +x run.sh && ./run.sh

LMP_EXE="${LMP:-lmp}"

echo "================================================"
echo "  Li-SNAP OpenKIM 验证测试"
echo "================================================"
echo "  LAMMPS executable: ${LMP_EXE}"
echo "================================================"
echo ""

${LMP_EXE} -in in.test_li

if [ $? -eq 0 ]; then
    echo ""
    echo "  Li-SNAP OpenKIM 验证成功完成"
    echo "  检查输出中的 Total Energy 和 Pressure 值"
fi