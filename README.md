# 晶界构造与 SOAP 方法

研究目标：建立锂金属晶界的结构无序表征，检验其与亚稳状态、微观变形机制和临界剪切强度的关联。

## 中期进展

汇报日期：2026-11-05；仓库整理更新：2026-11-15。此次更新整理已有记录，没有新增模拟或参数扫描。

| 内容 | 已有进展与限制 | 入口 |
|---|---|---|
| Li晶界构造与优化 | 已实现周期性双晶与热处理流程；既有短轨迹对照支持采用aniso降低残余偏应力 | [晶界构造](grain_boundary/README.md) |
| Ni SOAP参数探索 | 完成240组扫描；PPT口径下角向多数相邻比较通过，径向0/180通过，参数尚未冻结 | [扫描结果](soap/experiments/2026-11-01_midterm/REPORT.md) |
| 径向补充诊断 | 扩到14/16阶仍未建立收敛证据；polynomial对照存在高阶正交化及本机l=10旋转偏差 | [补充对照](soap/experiments/2026-11-01_radial_followup/REPORT.md) |
| 老师建议的新试验 | 扩大x、y尺寸，高温NPT保持2 ps后按10 K/ps淬火；10 K、零压稳定超过1 ps可继续研究 | [待实施方案](grain_boundary/current/TEACHER_PLAN_2026-11-05.md) |

[中期进展与证据索引](maintenance/docs/MIDTERM_2026-11.md) · [方法说明与答辩问答](maintenance/docs/MIDTERM_QA_2026-11.md) · [当前研究状态](maintenance/docs/RESEARCH_STATUS.md)

SOAP距离目前是相对体相参考的局域结构差异候选表征；它与亚稳程度、力学响应的关系仍需验证。老师给出的低温稳定判据用于决定构型是否可以继续研究，不要求先完成相变问题研究。

## 两条研究主线

本项目有两条研究主线：

| 主线 | 内容 | 入口 |
|---|---|---|
| 晶界构造 | Li Σ5(310)[001] 晶界的构建、热处理与结构验证 | [grain_boundary/](grain_boundary/README.md) |
| SOAP 方法 | 用 Ni 双晶探索局域结构描述、距离与无序度分析 | [soap/](soap/README.md) |

晶界构造中的基准、实验和历史版本统一收在 `grain_boundary/`；SOAP 的六阶段脚本、报告与输出统一收在 `soap/`。

台账、学习资料及维护工具见 [maintenance/](maintenance/README.md)。
