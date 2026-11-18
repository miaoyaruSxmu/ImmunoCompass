from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn


class BetaVAE(nn.Module):
    def __init__(self, input_dim, latent_dim, encoder_layers, decoder_layers):
        super().__init__()
        encoder = []
        current = input_dim
        for units in encoder_layers:
            encoder.extend([nn.Linear(current, units), nn.ReLU()])
            current = units
        self.encoder = nn.Sequential(*encoder)
        self.mu_layer = nn.Linear(current, latent_dim)
        self.logvar_layer = nn.Linear(current, latent_dim)

        decoder = []
        current = latent_dim
        for units in decoder_layers:
            decoder.extend([nn.Linear(current, units), nn.ReLU()])
            current = units
        decoder.append(nn.Linear(current, input_dim))
        self.decoder = nn.Sequential(*decoder)

    def encode(self, x):
        encoded = self.encoder(x)
        return self.mu_layer(encoded), self.logvar_layer(encoded)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        return mu + torch.randn_like(std) * std

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_load_model(path: Path):
    expected_globals = {
        "__main__.BetaVAE",
        "torch.nn.modules.container.Sequential",
        "torch.nn.modules.linear.Linear",
        "torch.nn.modules.activation.ReLU",
    }
    observed_globals = set(torch.serialization.get_unsafe_globals_in_checkpoint(path))
    if observed_globals != expected_globals:
        raise RuntimeError(
            "Unexpected serialized globals in model file: "
            + "; ".join(sorted(observed_globals - expected_globals))
        )
    allowlist = [BetaVAE, nn.Sequential, nn.Linear, nn.ReLU]
    with torch.serialization.safe_globals(allowlist):
        model = torch.load(path, map_location="cpu", weights_only=True)
    model.eval()
    return model, sorted(observed_globals)


def standardize(frame: pd.DataFrame, feature_names, mean, scale) -> tuple[np.ndarray, np.ndarray]:
    values = frame.loc[:, feature_names].to_numpy(dtype=np.float64, copy=True)
    if np.isinf(values).any():
        raise ValueError("Infinite values remain in the ssGSEA matrix")
    missing = np.isnan(values)
    if missing.any():
        row_index, column_index = np.where(missing)
        values[row_index, column_index] = mean[column_index]
    return (values - mean) / scale, missing


def encode_mu(model, values: np.ndarray, batch_size: int = 256) -> np.ndarray:
    outputs = []
    with torch.no_grad():
        for start in range(0, values.shape[0], batch_size):
            x = torch.as_tensor(values[start : start + batch_size], dtype=torch.float32)
            mu, _ = model.encode(x)
            outputs.append(mu.cpu().numpy())
    return np.vstack(outputs)


def scale_latent(values: np.ndarray, minimum: np.ndarray, span: np.ndarray) -> np.ndarray:
    return (values - minimum) / span


def defense_score(latent: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denominator = np.abs(reference).sum()
    if denominator <= 0:
        raise ValueError("The fixed real-control reference has a zero norm")
    distance = np.abs(latent - reference).sum(axis=1) / denominator
    return np.where(latent[:, 1] >= reference[1], 1.0 + distance, 1.0 - distance)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--ssgsea", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    parameter_dir = args.repo / "Data" / "params"
    data_dir = args.repo / "Data" / "data"
    model_path = (
        args.repo
        / "Data"
        / "anneal_betaFalse_anneal_steps20000enc512_128_32_lr0.0001_bs16_beta0.1.pth"
    )
    feature_path = parameter_dir / "feature_vector_1d.json"
    selected_features = json.loads(feature_path.read_text(encoding="utf-8"))
    if len(selected_features) != 619:
        raise ValueError("The fixed 619-feature panel is incomplete")

    group_files = sorted(data_dir.glob("ssgsea_*.csv.bz2"))
    if len(group_files) != 8:
        raise ValueError(f"Expected eight public ImmuDef training groups, found {len(group_files)}")
    training_frames = []
    training_groups = []
    for path in group_files:
        group = re.sub(
            r"_[Cc]7\.csv\.bz2$", "", path.name.removeprefix("ssgsea_")
        )
        frame = pd.read_csv(path, index_col=0).T
        if frame.shape[1] != 5219 or len(set(frame.columns)) != 5219:
            raise ValueError(f"{path.name} does not contain the expected 5219 C7 pathways")
        training_frames.append(frame)
        training_groups.extend([group] * frame.shape[0])
    training = pd.concat(training_frames, axis=0)
    feature_names = list(training.columns)
    if len(feature_names) != 5219 or len(set(feature_names)) != 5219:
        raise ValueError("Unexpected ImmuDef C7 training matrix structure")
    mean = training.loc[:, feature_names].to_numpy(dtype=np.float64).mean(axis=0)
    scale = training.loc[:, feature_names].to_numpy(dtype=np.float64).std(axis=0, ddof=0)
    if np.any(scale <= 0):
        raise ValueError("A public ImmuDef training pathway has zero variance")
    if not set(selected_features).issubset(feature_names):
        raise ValueError("The fixed 619-feature panel does not match public training pathways")
    selected_index = np.asarray([feature_names.index(x) for x in selected_features])

    external = pd.read_csv(args.ssgsea, compression="gzip")
    metadata_columns = ["sample_key", "sample_id", "cohort"]
    missing_pathways = sorted(set(feature_names) - set(external.columns))
    if missing_pathways:
        raise ValueError(f"External ssGSEA matrix lacks {len(missing_pathways)} pathways")
    external_scaled_all, external_missing = standardize(
        external, feature_names, mean, scale
    )
    external_scaled = external_scaled_all[:, selected_index]
    training_scaled_all, training_missing = standardize(
        training, feature_names, mean, scale
    )
    if training_missing.any():
        raise ValueError("Bundled ImmuDef training profiles unexpectedly contain missing pathways")
    training_scaled = training_scaled_all[:, selected_index]

    model, serialized_globals = safe_load_model(model_path)
    training_mu = encode_mu(model, training_scaled)
    external_mu = encode_mu(model, external_scaled)

    latent_minimum = training_mu.min(axis=0)
    latent_span = training_mu.max(axis=0) - latent_minimum
    if np.any(latent_span <= 0):
        raise ValueError("A training latent dimension has zero range")
    training_latent = scale_latent(training_mu, latent_minimum, latent_span)
    external_latent = scale_latent(external_mu, latent_minimum, latent_span)

    training_groups = np.asarray(training_groups)
    control_mask = training_groups == "Control"
    if control_mask.sum() != 368:
        raise ValueError(f"Expected 368 public controls, found {control_mask.sum()}")
    reference = np.median(training_latent[control_mask], axis=0)

    training_scores = defense_score(training_latent, reference)
    external_scores = defense_score(external_latent, reference)
    output = external.loc[:, metadata_columns].copy()
    output["ImmuDef_compatible_score"] = external_scores
    output["latent_z0"] = external_latent[:, 0]
    output["latent_z1"] = external_latent[:, 1]
    output["L1_distance_from_real_control_reference"] = np.abs(
        external_latent - reference
    ).sum(axis=1)
    output["missing_C7_features_imputed"] = external_missing.sum(axis=1)
    output["missing_VAE619_features_imputed"] = external_missing[:, selected_index].sum(axis=1)
    output.to_csv(args.outdir / "immudef_compatible_external_scores.csv", index=False)

    training_output = pd.DataFrame(
        {
            "training_group": training_groups,
            "ImmuDef_compatible_score": training_scores,
            "latent_z0": training_latent[:, 0],
            "latent_z1": training_latent[:, 1],
        },
        index=training.index,
    )
    training_output.index.name = "sample_id"
    training_output.to_csv(args.outdir / "immudef_public_training_scores.csv")
    summary = (
        training_output.groupby("training_group")["ImmuDef_compatible_score"]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    summary.to_csv(args.outdir / "immudef_public_training_group_summary.csv", index=False)

    provenance = {
        "score_name": "ImmuDef-compatible real-control-reference score",
        "official_DImmuScore_reproduced": False,
        "reason": (
            "The public repository does not include the generated high-confidence "
            "control reference or the three ensemble classifiers required by run_immudef.py."
        ),
        "reference_definition": "median VAE posterior mean of 368 bundled real controls",
        "latent_scaling": "min-max parameters fitted only on eight bundled ImmuDef training groups",
        "external_labels_used": False,
        "vae_representation": "posterior mean (mu), deterministic",
        "training_samples": int(training.shape[0]),
        "external_samples": int(external.shape[0]),
        "external_missing_C7_values_imputed": int(external_missing.sum()),
        "external_missing_VAE619_values_imputed": int(
            external_missing[:, selected_index].sum()
        ),
        "missing_value_rule": "pretrained scaler mean; standardized value equals zero",
        "scaler_features": len(feature_names),
        "vae_features": len(selected_features),
        "reference": reference.tolist(),
        "latent_minimum": latent_minimum.tolist(),
        "latent_span": latent_span.tolist(),
        "model_sha256": sha256(model_path),
        "scaler_sha256": sha256(parameter_dir / "scaler.pkl"),
        "scaler_reconstruction": (
            "mean and population standard deviation recomputed from the eight bundled"
            " public ImmuDef training ssGSEA matrices; no external labels used"
        ),
        "feature_list_sha256": sha256(feature_path),
        "allowed_serialized_globals": serialized_globals,
        "torch_version": torch.__version__,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
    }
    (args.outdir / "immudef_compatible_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
