"""Stable task identities; no inferred cross-model parameter transfer."""
PROJECTS = (
    "task01_lead_developability", "task02_tpd", "task03_covalent_kinetics",
    "task04_pbpk", "task05_cyp_ddi", "task06_asd", "task07_qsp",
    "task08_adc", "task09_cryoem_allostery", "task10_rna_splicing",
)
DRIVERS = (
    "run_task1_mpo_admet_developability.py", "run_task2_tpd_ternary_cooperativity.py",
    "run_task3_covalent_kinetics_residence_time.py", "run_task4_pbpk_pharmacokinetics_dose_prediction.py",
    "run_task5_cyp_ddi_mechanism_based_inhibition.py", "run_task6_asd_formulation_supersaturation_kinetics.py",
    "run_task7_qsp_tumor_immune_pkpd_synergy.py", "run_task8_adc_dar_cleavage_bystander_dynamics.py",
    "run_task9_cryoem_cryptic_pocket_allostery.py", "run_task10_rna_targeted_small_molecule_dynamics.py",
)
TITLES = (
    "先导可开发性 · MPO / ADMET", "靶向蛋白降解 · TPD", "共价抑制动力学 · Covalent kinetics",
    "生理药代动力学 · PBPK", "CYP 药物相互作用 · DDI", "无定形固体分散体 · ASD",
    "肿瘤免疫系统药理 · QSP", "抗体偶联药物 · ADC", "构象与变构 · Structural allostery",
    "RNA 结合与剪接 · RNA targeting",
)
FIGURES = (
    "fig_task1_developability_mpo.png", "fig_task2_tpd_hook_effect.png",
    "fig_task3_covalent_kinetics.png", "fig_task4_pbpk_pharmacokinetics.png",
    "fig_task5_cyp_ddi_mbi.png", "fig_task6_asd_supersaturation.png",
    "fig_task7_qsp_immuno_oncology.png", "fig_task8_adc_multiscale.png",
    "fig_task9_cryoem_allostery.png", "fig_task10_rna_targeted_cadd.png",
)
REPORTS = (
    "DEVELOPABILITY_MPO_REPORT", "TPD_TERNARY_COOPERATIVITY_REPORT", "COVALENT_DRUG_KINETICS_REPORT",
    "PBPK_DOSE_PREDICTION_REPORT", "CYP_DDI_KINETICS_REPORT", "ASD_FORMULATION_KINETICS_REPORT",
    "QSP_IMMUNO_ONCOLOGY_REPORT", "ADC_TRANSLATIONAL_ENGINEERING_REPORT",
    "CRYOEM_CRYPTIC_POCKET_REPORT", "RNA_TARGETED_CADD_REPORT",
)
BOUNDARIES = (
    "结构身份与描述符有记录；ADMET 为未校准代理，气相构象不证明变色龙性。",
    "平衡与合成降解情景；峰值不是临床最优剂量，连接子尚未包含真实蛋白出口向量。",
    "实际执行 EHT；速率参数为情景假设，GSH 筛查不等于 DILI 安全评价。",
    "简化分配模型与假设输入；第 7 天不自动等于稳态，剂量不是临床建议。",
    "实际调用 Task 4；M12 形式筛查基于假设参数，不构成临床禁忌判定。",
    "假设聚合物与析晶参数；吸收增强为代理，VFT 适用域之外保留缺失。",
    "50 名配对虚拟个体；质量进展终点不等于临床 RECIST / PFS。",
    "合成 DAR/HIC 与反应扩散情景；阈值半径不是经实验验证的杀伤半径。",
    "公开晶体结构 + 合成轨迹；状态自由能不是活化势垒，ROI 截断保留未知。",
    "ViennaRNA 计算 + 公开 NMR；工程化 RNA 和假设结合参数，U1 占有率不是临床剪接率。",
)
STAGES = ("discovery", "discovery", "discovery", "exposure", "exposure", "exposure",
          "response", "response", "structure", "structure")
