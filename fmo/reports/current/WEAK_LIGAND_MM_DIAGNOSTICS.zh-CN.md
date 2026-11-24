# 弱集合配体内部几何代价：匹配 MM 诊断

OpenFF 2.2.0 / NAGL；每个配体相同原子、参数和电荷，只换坐标。能量单位 kcal/mol。
这是无受体、无溶剂的诊断，不是结合自由能；绝不把不同分子的绝对总能量直接作排序。

| 集合 | n | B−基线内部能差 中位数 | 最小 / 最大差 | 局部最小化未达力阈值例数 |
|---|---:|---:|---|---:|
| CHK1 | 13 | -2.496 | -11.494 / -0.961 | 1 |
| BACE | 36 | -4.573 | -5.790 / 48.127 | 1 |
| CDK8 | 10 | 0.119 | -13.938 / 0.976 | 0 |
| HIV1PROT | 13 | 0.000 | 0.000 / 0.000 | 4 |
| Galectin3 | 26 | 0.036 | -2.453 / 0.508 | 0 |

## CHK1：最明显的内部能上升例

17 torsion adjust：B−基线 = -0.961 kcal/mol。
分项差：{'0:NonbondedForce': -0.44, '1:HarmonicBondForce': 0.098, '2:PeriodicTorsionForce': -0.427, '3:CMMotionRemover': 0.0, '4:HarmonicAngleForce': -0.192}。
最大角度变化：原子索引 [0, 3, 33]（['C', 'C', 'H']），117.10° → 118.03°。
孤立配体参照未达力阈值：12 docked；不能称其为已收敛局部最低点。

## BACE：最明显的内部能上升例

CAT-13i：B−基线 = 48.127 kcal/mol。
分项差：{'0:NonbondedForce': 0.845, '1:HarmonicBondForce': 0.678, '2:PeriodicTorsionForce': 4.914, '3:CMMotionRemover': 0.0, '4:HarmonicAngleForce': 41.69}。
最大角度变化：原子索引 [35, 9, 36]（['H', 'C', 'H']），113.70° → 66.13°。
孤立配体参照未达力阈值：CAT-17f；不能称其为已收敛局部最低点。

## CDK8：最明显的内部能上升例

13：B−基线 = 0.976 kcal/mol。
分项差：{'0:NonbondedForce': 0.612, '1:HarmonicBondForce': -0.045, '2:PeriodicTorsionForce': 0.455, '3:CMMotionRemover': 0.0, '4:HarmonicAngleForce': -0.047}。
最大角度变化：原子索引 [27, 7, 28]（['H', 'N', 'H']），119.43° → 116.68°。

## HIV1PROT：最明显的内部能上升例

1：B−基线 = 0.000 kcal/mol。
分项差：{'0:NonbondedForce': 0.0, '1:HarmonicBondForce': 0.0, '2:PeriodicTorsionForce': 0.0, '3:CMMotionRemover': 0.0, '4:HarmonicAngleForce': 0.0}。
最大角度变化：原子索引 [1, 0, 14]（['C', 'C', 'C']），108.41° → 108.41°。
孤立配体参照未达力阈值：1, 14, 15, 20；不能称其为已收敛局部最低点。

## Galectin3：最明显的内部能上升例

18a：B−基线 = 0.508 kcal/mol。
分项差：{'0:NonbondedForce': 0.061, '1:HarmonicBondForce': -0.012, '2:PeriodicTorsionForce': 0.475, '3:CMMotionRemover': 0.0, '4:HarmonicAngleForce': -0.017}。
最大角度变化：原子索引 [2, 51, 52]（['C', 'O', 'H']），112.32° → 112.13°。

## 解读

BACE CAT-13i 出现明显内部角度代价，说明仅审计重原子碰撞和酰胺平面仍可能漏掉含氢键角异常。
该异常不能独自解释整个 BACE 集合错排；需要完整姿态保持与氢约束对照，不删除异常配体。
CHK1 内部代价下降也没有同步带来排序提升：局部模型几何改善不是亲和力排序的充分条件。
六例孤立参照未达力阈值时仅保留近似参照状态；本表直接比较两个绑定几何，差值不依赖该参照。
不按内部 MM 最低点选择 tPIE，更不利用实验标签筛分支。
