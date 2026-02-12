from __future__ import annotations

from pathlib import Path
from typing import Optional

import streamlit as st

from mco_mnox.model import load_model, predict, train_model
from mco_mnox.utils import load_config, set_seed, setup_logging


def _safe_path(p: str) -> Optional[Path]:
    p = (p or "").strip()
    if not p:
        return None
    return Path(p)


def _file_ok(path: Optional[Path], required: bool = True) -> bool:
    if path is None:
        return not required
    return path.exists() and path.is_file()


def run_app() -> None:
    st.set_page_config(page_title="MCO MnOx PU GUI", layout="wide")
    st.title("MCO Mn(II)-oxidizing 倾向预测 GUI")
    st.caption("PU 学习：只把 positives 当已知阳性，unlabeled 视为混合分布，不直接当负例。")

    tab_train, tab_predict, tab_about = st.tabs(["训练", "预测", "说明"])

    with tab_train:
        st.subheader("训练模型")
        c1, c2 = st.columns(2)
        with c1:
            pos = st.text_input("positives.fasta 路径", "data/toy/positives.fasta")
            unl = st.text_input("unlabeled.fasta 路径", "data/toy/unlabeled.fasta")
            neg = st.text_input("negatives.fasta 路径（可选）", "data/toy/negatives.fasta")
        with c2:
            outdir = st.text_input("输出目录", "runs/gui_run")
            config_path = st.text_input("配置文件", "configs/default.yaml")
            seed = st.number_input("随机种子", min_value=0, value=42, step=1)

        strategy = st.selectbox("PU 策略", ["bagging", "rn"], index=0)
        embedder = st.selectbox("Embedding", ["none", "esm2_t12_35M", "esm2_t33_650M", "protT5"], index=0)

        if st.button("开始训练", type="primary"):
            pos_p, unl_p, neg_p = _safe_path(pos), _safe_path(unl), _safe_path(neg)
            outdir_p, cfg_p = _safe_path(outdir), _safe_path(config_path)

            if not _file_ok(pos_p, True):
                st.error("positives.fasta 路径不存在。")
            elif not _file_ok(unl_p, True):
                st.error("unlabeled.fasta 路径不存在。")
            elif not _file_ok(cfg_p, True):
                st.error("配置文件路径不存在。")
            elif neg_p is not None and not _file_ok(neg_p, False):
                st.error("negatives.fasta 路径不存在。")
            else:
                cfg = load_config(cfg_p)
                cfg["seed"] = int(seed)
                cfg["pu"]["strategy"] = strategy
                cfg["features"]["embedder"] = embedder
                outdir_p.mkdir(parents=True, exist_ok=True)

                set_seed(int(seed))
                setup_logging(outdir_p)

                with st.spinner("训练中，请稍候..."):
                    model = train_model(
                        pos_fasta=str(pos_p),
                        unl_fasta=str(unl_p),
                        neg_fasta=str(neg_p) if neg_p and neg_p.exists() else None,
                        outdir=str(outdir_p),
                        config=cfg,
                    )
                st.success(f"训练完成。模型保存至: {outdir_p / 'model.pkl'}")
                st.json(
                    {
                        "strategy": model.strategy,
                        "n_positives": len(model.positives),
                        "feature_dim": len(model.feature_names),
                        "embedder": embedder,
                    }
                )

    with tab_predict:
        st.subheader("加载模型并预测")
        c1, c2 = st.columns(2)
        with c1:
            model_path = st.text_input("model.pkl 路径", "runs/gui_run/model.pkl")
            fasta_path = st.text_input("待预测 FASTA", "data/toy/unlabeled.fasta")
        with c2:
            out_tsv = st.text_input("输出 TSV", "runs/gui_run/results.tsv")

        if st.button("开始预测", type="primary"):
            m_p = _safe_path(model_path)
            f_p = _safe_path(fasta_path)
            o_p = _safe_path(out_tsv)
            if not _file_ok(m_p, True):
                st.error("model.pkl 路径不存在。")
            elif not _file_ok(f_p, True):
                st.error("FASTA 路径不存在。")
            else:
                o_p.parent.mkdir(parents=True, exist_ok=True)
                with st.spinner("预测中..."):
                    model = load_model(str(m_p))
                    df, rec = predict(model=model, fasta=str(f_p), out_tsv=str(o_p))
                st.success(f"预测完成。结果写入: {o_p} 和 {o_p.with_suffix('.json')}")
                st.dataframe(df, use_container_width=True)
                st.subheader("实验优先级建议")
                st.json(rec)

    with tab_about:
        st.markdown(
            """
            ### 这个 GUI 做了什么
            - 训练：读取 positives / unlabeled / (可选) negatives，执行 PU 学习。
            - 预测：对未知序列输出 score_mean / score_std / motif / 定位代理等字段。
            - 输出：`results.tsv` + `results.json`（候选优先级）。

            ### 注意
            - 这是 PU 排序任务，不是传统全监督二分类。
            - 推荐按 top-N 候选做实验验证，并将新标签回流迭代训练。
            """
        )


def main() -> None:
    run_app()


if __name__ == "__main__":
    main()
