#!/usr/bin/env python3
"""Reproduce MaleCNS mesh-proximity calculations used for Fig. 7M/O.

This version adds explicit progress/timing output and diagnostic modes so that
slow remote mesh retrievals can be localized. Distances are screening metrics
derived from segmentation meshes. They are not exact membrane-to-membrane
distances and do not establish synaptic connectivity.
"""

from __future__ import annotations

import argparse
import inspect
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from cloudvolume import CloudVolume
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
META = ROOT / "metadata" / "fig7_representative_sites.csv"

MESH_SOURCE = "precomputed://gs://flyem-male-cns/v1.0/segmentation"
VOXEL_SIZE_NM = np.array([8.0, 8.0, 8.0])
EXPECTED_FIG7M_NM = 19.595918


def log(message: str) -> None:
    print(message, flush=True)


def elapsed(start: float) -> str:
    return f"{time.perf_counter() - start:.1f} s"


def mesh_get_kwargs(vol: CloudVolume, lod: int) -> dict:
    """Build only keyword arguments supported by the installed CloudVolume."""
    try:
        params = inspect.signature(vol.mesh.get).parameters
    except (TypeError, ValueError):
        params = {}
    kwargs = {}
    if "lod" in params:
        kwargs["lod"] = lod
    if "progress" in params:
        kwargs["progress"] = True
    return kwargs


def fetch_vertices(vol: CloudVolume, body_id: int, lod: int = 0) -> np.ndarray:
    log(f"[mesh] Fetching body {body_id} at LOD {lod} ...")
    start = time.perf_counter()
    kwargs = mesh_get_kwargs(vol, lod)
    if "lod" not in kwargs and lod != 0:
        raise RuntimeError(
            "This CloudVolume mesh backend does not expose an LOD argument; "
            "diagnostic LOD > 0 cannot be requested."
        )
    mesh = vol.mesh.get(int(body_id), **kwargs)
    log(f"[mesh] Download/decode finished for {body_id} ({elapsed(start)}).")

    if isinstance(mesh, dict):
        if int(body_id) in mesh:
            mesh = mesh[int(body_id)]
        else:
            mesh = next(iter(mesh.values()))

    vertices = np.asarray(mesh.vertices, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise RuntimeError(f"Unexpected mesh vertex shape for {body_id}: {vertices.shape}")

    log(
        f"[mesh] body {body_id}: {len(vertices):,} vertices; "
        f"array size ~{vertices.nbytes / (1024**2):.1f} MiB"
    )
    return vertices


def probe_manifest(vol: CloudVolume, body_id: int) -> None:
    log(f"[probe] Checking mesh manifest for body {body_id} ...")
    start = time.perf_counter()
    try:
        if hasattr(vol.mesh, "get_manifest"):
            manifest = vol.mesh.get_manifest(int(body_id))
            if manifest is None:
                log(f"[probe] body {body_id}: manifest not found ({elapsed(start)}).")
                return
            attrs = []
            for name in ("num_lods", "segment_id", "chunk_shape", "grid_origin"):
                if hasattr(manifest, name):
                    attrs.append(f"{name}={getattr(manifest, name)}")
            details = ", ".join(attrs) if attrs else type(manifest).__name__
            log(f"[probe] body {body_id}: manifest OK; {details} ({elapsed(start)}).")
        elif hasattr(vol.mesh, "exists"):
            exists = vol.mesh.exists([int(body_id)], progress=True)
            log(f"[probe] body {body_id}: mesh existence result={exists} ({elapsed(start)}).")
        else:
            log("[probe] Mesh backend exposes neither get_manifest() nor exists().")
    except Exception as exc:
        log(f"[probe] body {body_id}: ERROR {type(exc).__name__}: {exc}")


def min_mesh_distance_nm(a_nm: np.ndarray, b_nm: np.ndarray) -> dict:
    """Nearest-neighbor distance between two vertex sets using cKDTree."""
    log(
        f"[calc] Nearest-neighbor search: A={len(a_nm):,} vertices, "
        f"B={len(b_nm):,} vertices ..."
    )
    start = time.perf_counter()
    if len(a_nm) <= len(b_nm):
        tree = cKDTree(a_nm)
        d, idx = tree.query(b_nm, k=1, workers=-1)
        j = int(np.argmin(d))
        i = int(idx[j])
        pa, pb = a_nm[i], b_nm[j]
    else:
        tree = cKDTree(b_nm)
        d, idx = tree.query(a_nm, k=1, workers=-1)
        i = int(np.argmin(d))
        j = int(idx[i])
        pa, pb = a_nm[i], b_nm[j]
    distance = float(np.linalg.norm(pa - pb))
    log(f"[calc] Nearest-neighbor search finished ({elapsed(start)}).")
    return {
        "distance_nm": distance,
        "distance_um": distance / 1000.0,
        "point_a_nm": pa.tolist(),
        "point_b_nm": pb.tolist(),
    }


def row_coord(meta: pd.DataFrame, panel: str, role: str) -> np.ndarray:
    row = meta[(meta["panel"] == panel) & (meta["coordinate_role"] == role)]
    if len(row) != 1:
        raise RuntimeError(f"Expected one metadata row for panel={panel}, role={role}")
    r = row.iloc[0]
    return np.array([r.x_voxel, r.y_voxel, r.z_voxel], dtype=float)


def output_name(base: str, lod: int) -> str:
    if lod == 0:
        return base
    stem, suffix = base.rsplit(".", 1)
    return f"{stem}_LOD{lod}_DIAGNOSTIC.{suffix}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--probe",
        action="store_true",
        help="Only inspect mesh availability/manifests; do not download full meshes.",
    )
    parser.add_argument(
        "--panel",
        choices=("M", "O", "all"),
        default="all",
        help="Run only Fig. 7M, only Fig. 7O, or both (default: all).",
    )
    parser.add_argument(
        "--lod",
        type=int,
        default=0,
        help="Mesh level of detail. LOD 0 is the manuscript/final analysis. Higher LODs are diagnostic only.",
    )
    args = parser.parse_args()

    if args.lod < 0:
        raise SystemExit("--lod must be >= 0")

    OUTPUTS.mkdir(exist_ok=True)
    log("[init] Starting MaleCNS mesh-proximity analysis.")
    log(f"[init] Mesh source: {MESH_SOURCE}")
    log(f"[init] Requested panel: {args.panel}; LOD: {args.lod}")
    if args.lod != 0:
        log("[WARNING] LOD > 0 is DIAGNOSTIC ONLY and must not be used for manuscript distances.")

    start = time.perf_counter()
    vol = CloudVolume(
        MESH_SOURCE,
        use_https=True,
        progress=True,
        parallel=1,
        cache=False,
    )
    log(f"[init] CloudVolume initialized ({elapsed(start)}).")
    log(f"[init] Mesh backend: {type(vol.mesh).__name__}")
    log(f"[init] info['mesh']: {vol.info.get('mesh', '<not listed>')}")

    if args.probe:
        for body_id in (512245, 535798, 10342):
            probe_manifest(vol, body_id)
        log("[probe] Probe finished. No full mesh was downloaded.")
        return

    meta = pd.read_csv(META)
    kc = None

    if args.panel in ("M", "all"):
        log("[M] Fig. 7M: fetching SIFa 512245 and KC 535798.")
        sifa_m = fetch_vertices(vol, 512245, args.lod)
        kc = fetch_vertices(vol, 535798, args.lod)
        m_result = min_mesh_distance_nm(sifa_m, kc)
        m_result.update(
            {
                "panel": "M",
                "sifa_body_id": 512245,
                "kc_body_id": 535798,
                "lod": args.lod,
                "interpretation": "mesh-screened close apposition; not classified as a synapse",
            }
        )
        log(
            f"[M] Minimum distance: {m_result['distance_nm']:.6f} nm "
            f"({m_result['distance_um']:.6f} um)"
        )
        if args.lod == 0:
            delta = abs(m_result["distance_nm"] - EXPECTED_FIG7M_NM)
            log(
                f"[M] Expected prior value: {EXPECTED_FIG7M_NM:.6f} nm; "
                f"absolute difference: {delta:.6f} nm"
            )
        with open(
            OUTPUTS / output_name("fig7M_mesh_proximity.json", args.lod),
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(m_result, f, indent=2)
        log("[M] Output written.")

    if args.panel in ("O", "all"):
        log("[O] Fig. 7O: preparing soma-associated proximity analysis.")
        if kc is None:
            kc = fetch_vertices(vol, 535798, args.lod)

        soma_vox = row_coord(meta, "O", "KC_soma")
        sifa_center_vox = row_coord(meta, "O", "SIFa_selected_center")
        soma_nm = soma_vox * VOXEL_SIZE_NM
        sifa_center_nm = sifa_center_vox * VOXEL_SIZE_NM

        soma_radius_nm = 3000.0
        log("[O] Selecting KC mesh vertices within 3 um of annotated soma coordinate ...")
        kc_soma_mask = np.linalg.norm(kc - soma_nm[None, :], axis=1) <= soma_radius_nm
        kc_soma_vertices = kc[kc_soma_mask]
        log(f"[O] KC soma-associated vertices: {len(kc_soma_vertices):,}")
        if len(kc_soma_vertices) == 0:
            raise RuntimeError("No KC mesh vertices found within 3 um of the annotated soma coordinate")

        sifa_o = fetch_vertices(vol, 10342, args.lod)
        local_rows = []
        for radius_um in (0.5, 1.0, 1.5):
            radius_nm = radius_um * 1000.0
            log(f"[O] Testing SIFa local radius {radius_um:.1f} um ...")
            local_mask = np.linalg.norm(sifa_o - sifa_center_nm[None, :], axis=1) <= radius_nm
            local_sifa = sifa_o[local_mask]
            log(f"[O]   local SIFa vertices: {len(local_sifa):,}")
            if len(local_sifa) == 0:
                local_rows.append(
                    {
                        "panel": "O",
                        "radius_um": radius_um,
                        "n_sifa_vertices": 0,
                        "n_kc_soma_vertices": int(len(kc_soma_vertices)),
                        "distance_nm": np.nan,
                        "distance_um": np.nan,
                    }
                )
                continue
            result = min_mesh_distance_nm(local_sifa, kc_soma_vertices)
            local_rows.append(
                {
                    "panel": "O",
                    "radius_um": radius_um,
                    "n_sifa_vertices": int(len(local_sifa)),
                    "n_kc_soma_vertices": int(len(kc_soma_vertices)),
                    "distance_nm": result["distance_nm"],
                    "distance_um": result["distance_um"],
                }
            )
            log(f"[O]   minimum distance: {result['distance_um']:.6f} um")

        center_to_soma_nm = float(np.linalg.norm(sifa_center_nm - soma_nm))
        log(
            f"[O] Selected-center to annotated-soma coordinate: "
            f"{center_to_soma_nm:.3f} nm ({center_to_soma_nm / 1000.0:.6f} um)"
        )

        pd.DataFrame(local_rows).to_csv(
            OUTPUTS / output_name("fig7O_soma_associated_proximity.csv", args.lod),
            index=False,
        )
        pd.DataFrame(
            [
                {
                    "panel": "O",
                    "sifa_body_id": 10342,
                    "kc_body_id": 535798,
                    "lod": args.lod,
                    "sifa_selected_center_x_voxel": sifa_center_vox[0],
                    "sifa_selected_center_y_voxel": sifa_center_vox[1],
                    "sifa_selected_center_z_voxel": sifa_center_vox[2],
                    "kc_soma_x_voxel": soma_vox[0],
                    "kc_soma_y_voxel": soma_vox[1],
                    "kc_soma_z_voxel": soma_vox[2],
                    "center_to_soma_distance_nm": center_to_soma_nm,
                    "center_to_soma_distance_um": center_to_soma_nm / 1000.0,
                    "interpretation": "vicinity of KC soma-associated region; not demonstrated direct soma contact",
                }
            ]
        ).to_csv(
            OUTPUTS / output_name("fig7O_site_metadata_reproduced.csv", args.lod),
            index=False,
        )
        log("[O] Outputs written.")

    log(f"[done] Wrote outputs to: {OUTPUTS}")


if __name__ == "__main__":
    main()

