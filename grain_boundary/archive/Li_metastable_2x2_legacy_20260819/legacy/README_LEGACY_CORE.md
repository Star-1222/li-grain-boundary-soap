# README_LEGACY_CORE

旧版权威项目 `v2_production` 是一个已验证的 BCC Li Σ5 [001](310) 对称
倾转晶界工作流：

```text
3600 atoms, K=6 J=3 L=10
63 automated tests passed
Stage 0-4 completed: minimize -> 300 K NVT 20 ps anneal -> 300->1 K
50 ps quench -> final 0 K minimization
a = 3.506186475312 A
potential = SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

旧版流程的权威说明见 `README_v2_final.md`，环境与最终状态见
`reports_v2/ENVIRONMENT.txt` 与 `reports_v2/FINAL_STATUS.md`。

本目录只作为参考与溯源，不参与新实验的计算；新实验使用
`Li_Sigma5_310_metastable_2x2` 根目录下自己的生成器副本、LAMMPS 输入
和输出目录。
