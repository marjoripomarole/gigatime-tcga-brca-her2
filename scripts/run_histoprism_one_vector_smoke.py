#!/usr/bin/env python3
"""Run a minimal HistoPrism smoke test on one precomputed UNI feature vector."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import torch
import torch.nn as nn
import yaml
from huggingface_hub import hf_hub_download


MARKER_GENES = [
    "ERBB2",
    "ESR1",
    "PGR",
    "MKI67",
    "EPCAM",
    "KRT8",
    "KRT18",
    "PTPRC",
    "CD3D",
    "CD8A",
    "MS4A1",
    "COL1A1",
]


class CrossAttentionLayer(nn.Module):
    def __init__(self, query_dim: int, context_dim: int, n_heads: int, dropout_rate: float):
        super().__init__()
        self.attention = nn.MultiheadAttention(
            embed_dim=query_dim,
            kdim=context_dim,
            vdim=context_dim,
            num_heads=n_heads,
            dropout=dropout_rate,
            batch_first=True,
        )
        self.norm = nn.LayerNorm(query_dim)

    def forward(self, query: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        attn_output, _ = self.attention(query=query, key=context, value=context)
        return self.norm(query + attn_output)


class MinimalHistoPrism(nn.Module):
    """Subset of HistoPrism needed for split-0 inference with zero graph layers.

    The official split-0 config sets graph_config.n_conv_layers to 0, so graph
    convolution weights are present in the checkpoint but not used in forward.
    Keeping this smoke model minimal avoids adding torch-geometric just to inspect
    one prediction vector.
    """

    def __init__(self, config: dict, gene_dim: int, emb_dim: int = 1024):
        super().__init__()
        hidden_dim = int(config.get("hidden_dim", 256))
        dropout_rate = float(config.get("dropout_rate", 0.1))
        self.input_proj = nn.Linear(emb_dim, hidden_dim)
        self.onco_embedder = nn.Linear(int(config.get("onco_onehot_dim", 29)), emb_dim)
        self.conditioner = CrossAttentionLayer(
            query_dim=emb_dim,
            context_dim=emb_dim,
            n_heads=int(config.get("cross_attn_heads", 4)),
            dropout_rate=dropout_rate,
        )
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=hidden_dim,
            nhead=int(config.get("transformer_heads", 8)),
            batch_first=True,
            dropout=dropout_rate,
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer,
            num_layers=int(config.get("transformer_layers", 2)),
        )
        self.gene_regression_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, gene_dim),
        )

    def forward(self, emb: torch.Tensor, onco_onehot: torch.Tensor) -> torch.Tensor:
        onco_embedding = self.onco_embedder(onco_onehot).unsqueeze(0)
        query = emb.unsqueeze(0)
        conditioned_emb = self.conditioner(query=query, context=onco_embedding).squeeze(0)
        h = self.input_proj(conditioned_emb)
        h = self.transformer(h.unsqueeze(0)).squeeze(0)
        return self.gene_regression_head(h)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histoprism-repo", default="external/HistoPrism")
    parser.add_argument(
        "--embedding-file",
        default="external/HistoPrism/data/sample_processed_HEST1K/uni/ZEN48_embeddings.pt",
        help="Official HistoPrism-format dict of barcode -> 1024-dim UNI feature.",
    )
    parser.add_argument("--barcode", default="", help="Barcode key to use. Defaults to the first vector.")
    parser.add_argument("--oncotree-code", default="auto", help="Use auto to read the sample split table when possible.")
    parser.add_argument("--checkpoint", default="", help="Local checkpoint path. Defaults to HuSusu/HistoPrism split 0.")
    parser.add_argument("--out-dir", default="results/histoprism_one_vector_smoke")
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--device", choices=["auto", "cpu", "mps", "cuda"], default="auto")
    return parser.parse_args()


def resolve_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_yaml(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def sample_id_from_embedding_file(path: Path) -> str:
    name = path.name
    return name.removesuffix("_embeddings.pt")


def infer_oncotree_code(repo: Path, sample_id: str, requested: str) -> str:
    if requested != "auto":
        return requested
    for split_dir in [repo / "data/split_0", repo / "data/split_1"]:
        for split_name in ["train_split.csv", "val_split.csv", "test_split.csv"]:
            split_path = split_dir / split_name
            if not split_path.exists():
                continue
            with split_path.open("r", newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    if row.get("sample_id") == sample_id:
                        return row.get("oncotree_code") or row.get("oncotree") or "Unknown"
    return "Unknown"


def onehot_oncotree(code: str, codes: list[str], device: torch.device) -> torch.Tensor:
    onehot = torch.zeros((1, len(codes)), dtype=torch.float32, device=device)
    if code in codes:
        onehot[0, codes.index(code)] = 1.0
    return onehot


def symbol_lookup(repo: Path) -> dict[str, str]:
    mapping_path = repo / "data/gene_panel/STPath_genes.json"
    if not mapping_path.exists():
        return {}
    raw = json.loads(mapping_path.read_text(encoding="utf-8"))
    by_ensembl: dict[str, str] = {}
    for symbol, ensembl in raw.items():
        if not isinstance(symbol, str) or not isinstance(ensembl, str):
            continue
        if symbol.startswith("ENSG") or "__" in symbol or "." in symbol:
            continue
        by_ensembl.setdefault(ensembl, symbol)
    return by_ensembl


def rows_for_predictions(genes: list[str], predictions: torch.Tensor, symbols: dict[str, str]) -> list[dict[str, object]]:
    values = predictions.detach().cpu().float().numpy().tolist()
    rows = []
    for idx, (gene_id, value) in enumerate(zip(genes, values, strict=True)):
        rows.append(
            {
                "rank_input_order": idx,
                "gene_id": gene_id,
                "gene_symbol": symbols.get(gene_id, ""),
                "prediction": float(value),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    repo = Path(args.histoprism_repo)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = load_yaml(repo / "checkpoints_paper/split_0_paper/config.yml")
    genes = load_yaml(repo / "data/gene_panel/stpath_genes.yaml")["genes"]
    oncotree_codes = (repo / "data/split_0/oncotree_code_list.txt").read_text(encoding="utf-8").splitlines()
    symbols = symbol_lookup(repo)

    embedding_file = Path(args.embedding_file)
    sample_id = sample_id_from_embedding_file(embedding_file)
    oncotree_code = infer_oncotree_code(repo, sample_id, args.oncotree_code)
    if oncotree_code not in oncotree_codes:
        raise ValueError(f"Oncotree code {oncotree_code!r} is not in split_0 oncotree list.")

    embedding_dict = torch.load(embedding_file, map_location="cpu", weights_only=True)
    barcode = args.barcode or next(iter(embedding_dict))
    emb = embedding_dict[barcode].detach().float()
    if emb.ndim != 1 or emb.shape[0] != 1024:
        raise ValueError(f"Expected one 1024-dim UNI vector, got shape {tuple(emb.shape)}")

    checkpoint = args.checkpoint or hf_hub_download("HuSusu/HistoPrism", "HistoPrism_split0.ckpt")
    ckpt = torch.load(checkpoint, map_location="cpu", weights_only=False)
    state_dict = ckpt["model_state"]

    device = resolve_device(args.device)
    model = MinimalHistoPrism(config=config, gene_dim=len(genes), emb_dim=1024)
    load_result = model.load_state_dict(state_dict, strict=False)
    model.to(device)
    model.eval()

    with torch.inference_mode():
        predictions = model(
            emb=emb.unsqueeze(0).to(device),
            onco_onehot=onehot_oncotree(oncotree_code, oncotree_codes, device),
        ).squeeze(0)

    rows = rows_for_predictions(genes, predictions, symbols)
    high_rows = sorted(rows, key=lambda row: float(row["prediction"]), reverse=True)[: args.top_n]
    low_rows = sorted(rows, key=lambda row: float(row["prediction"]))[: args.top_n]
    marker_set = set(MARKER_GENES)
    marker_rows = [row for row in rows if row["gene_symbol"] in marker_set or row["gene_id"] in marker_set]

    write_csv(out_dir / "all_gene_predictions.csv", rows)
    write_csv(out_dir / "top_predicted_genes.csv", high_rows)
    write_csv(out_dir / "bottom_predicted_genes.csv", low_rows)
    write_csv(out_dir / "selected_marker_predictions.csv", marker_rows)

    erbb2_row = next((row for row in rows if row["gene_id"] == "ENSG00000141736"), None)
    summary = {
        "model": "HistoPrism split 0 minimal smoke",
        "checkpoint": str(checkpoint),
        "input_kind": "precomputed UNI feature vector from official HistoPrism sample data",
        "embedding_file": str(embedding_file),
        "sample_id": sample_id,
        "barcode": barcode,
        "oncotree_code": oncotree_code,
        "device": str(device),
        "n_genes": len(rows),
        "erbb2_prediction": erbb2_row,
        "missing_keys": list(load_result.missing_keys),
        "unexpected_keys": list(load_result.unexpected_keys),
        "outputs": {
            "all_gene_predictions": str(out_dir / "all_gene_predictions.csv"),
            "top_predicted_genes": str(out_dir / "top_predicted_genes.csv"),
            "bottom_predicted_genes": str(out_dir / "bottom_predicted_genes.csv"),
            "selected_marker_predictions": str(out_dir / "selected_marker_predictions.csv"),
        },
    }
    (out_dir / "histoprism_one_vector_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
