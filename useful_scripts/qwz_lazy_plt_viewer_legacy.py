#!/usr/bin/env python3
"""
Lazy interactive viewer for SpecLoc plots.

Workflow:
1) Build packet:
   julia simulations/data_processing/lazy_packet_builder.jl <input_root> <output_root>
2) Launch viewer on packet folder or packet records file:
   python3 useful_scripts/qwz_lazy_plt_viewer.py <.../lazy_packet or .../lazy_packet_records.tsv>

The viewer requests only visible plots and caches rendered images locally.
Background worker continuously generates plots, spreading from default parameters.
"""

import csv
import json
import math
import os
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path
from queue import PriorityQueue

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (
    QApplication,
    QGridLayout,
    QHBoxLayout,
    QLineEdit,
    QLabel,
    QMainWindow,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QVBoxLayout,
    QWidget,
)


PREFERRED_PLOT_ORDER = [
    "ribbon_ipr_fixed",
    "ribbon_dcdE_fixed",
    "ribbon_chern_acc_fixed",
    "ribbon_ipr_auto",
    "ribbon_dcdE_auto",
    "ribbon_chern_acc_auto",
    "specloc_spectrum_vs_gamma",
    "specloc_signature_vs_E",
    "specloc_loggap_vs_E",
    "band3d",
    "ldos_target",
    "dos",
    "specloc_gap_gamma_E",
    "specloc_sig_gamma_E",
    "specloc_gap_gamma_kappa",
    "specloc_gap_gamma_W",
    "specloc_gap_W_kappa",
    "specloc_sig_gamma_kappa",
    "specloc_sig_gamma_W",
    "specloc_sig_W_kappa",
    "specloc_spectrum_vs_W",
    "specloc_spectrum_vs_kappa",
]

PLOT_TITLES = {
    "ribbon_ipr_fixed": "Ribbon Spectrum (IPR, fixed)",
    "ribbon_dcdE_fixed": "Ribbon Spectrum (dC/dE, fixed)",
    "ribbon_chern_acc_fixed": "Accumulated Chern (fixed)",
    "ribbon_ipr_auto": "Ribbon Spectrum (IPR, auto)",
    "ribbon_dcdE_auto": "Ribbon Spectrum (dC/dE, auto)",
    "ribbon_chern_acc_auto": "Accumulated Chern (auto)",
    "band3d": "3D Bandstructure",
    "ldos_target": "LDOS at Target Energy",
    "ldos_lowest": "LDOS of Lowest-|E| States",
    "dos": "DOS vs Energy",
    "specloc_signature_vs_E": "Speclocaliser Signature vs Energy",
    "specloc_loggap_vs_E": "Speclocaliser Gap vs Energy",
    "specloc_gap_gamma_E": "Speclocaliser Gap: gamma vs E",
    "specloc_sig_gamma_E": "Speclocaliser Signature: gamma vs E",
    "specloc_gap_gamma_kappa": "Speclocaliser Gap: gamma vs kappa",
    "specloc_gap_gamma_W": "Speclocaliser Gap: gamma vs W",
    "specloc_gap_W_kappa": "Speclocaliser Gap: W vs kappa",
    "specloc_sig_gamma_kappa": "Speclocaliser Signature: gamma vs kappa",
    "specloc_sig_gamma_W": "Speclocaliser Signature: gamma vs W",
    "specloc_sig_W_kappa": "Speclocaliser Signature: W vs kappa",
    "specloc_spectrum_vs_gamma": "Speclocaliser Spectrum vs gamma",
    "specloc_spectrum_vs_W": "Speclocaliser Spectrum vs W",
    "specloc_spectrum_vs_kappa": "Speclocaliser Spectrum vs kappa",
}

RIBBON_BASE_TYPES = ["ribbon_ipr", "ribbon_dcdE", "ribbon_chern_acc"]

# For large packets, keep only essential panels enabled by default.
DISABLED_PLOT_TYPES = {
    "ldos_lowest",
    "specloc_gap_gamma_kappa",
    "specloc_gap_gamma_W",
    "specloc_gap_W_kappa",
    "specloc_sig_gamma_kappa",
    "specloc_sig_gamma_W",
    "specloc_sig_W_kappa",
    "specloc_spectrum_vs_W",
    "specloc_spectrum_vs_kappa",
}

VIRTUAL_SPEClOC_PLOT_TYPES = {
    "specloc_gap_gamma_E",
    "specloc_sig_gamma_E",
}


class PlotGenerationWorker(threading.Thread):
    """Background worker thread that systematically generates and caches plots."""

    def __init__(self, viewer):
        super().__init__(daemon=True)
        self.viewer = viewer
        self.request_queue = PriorityQueue()  # (priority, record_id, plot_type, target_params)
        self.running = True
        self.generated_index = self._load_generated_index()
        self.lock = threading.Lock()
        self.exploration_depth = 0
        self.exploration_visited = set()

    def _load_generated_index(self):
        """Load index of already-generated plots from persistent storage."""
        index_file = self.viewer.cache_dir / "background_generated_index.json"
        if index_file.exists():
            try:
                with open(index_file) as f:
                    data = json.load(f)
                    return set(data.get("generated", []))
            except Exception:
                return set()
        return set()

    def _save_generated_index(self):
        """Persist the index of generated plots."""
        index_file = self.viewer.cache_dir / "background_generated_index.json"
        try:
            with open(index_file, 'w') as f:
                json.dump({"generated": sorted(list(self.generated_index))}, f)
        except Exception:
            pass

    def request_plot(self, record_id, plot_type, priority=0):
        """Queue a plot generation request. priority=0 is user (high), higher=background (low)."""
        self.request_queue.put((priority, record_id, plot_type))

    def run(self):
        """Main worker loop: process user requests first, then explore background."""
        while self.running:
            try:
                # Check for user-prioritized requests (blocking, 0.1s timeout)
                try:
                    priority, record_id, plot_type = self.request_queue.get(timeout=0.1)
                    key = f"{record_id}:{plot_type}"

                    if key not in self.generated_index:
                        pairs = [(record_id, plot_type)]
                        self.viewer._render_pairs(pairs)
                        with self.lock:
                            self.generated_index.add(key)
                            self._save_generated_index()
                    continue
                except:
                    pass

                # No user requests; do background exploration
                self._explore_background()

            except Exception as e:
                if not self.running:
                    break
                time.sleep(0.5)

    def _explore_background(self):
        """Systematically explore parameter space from default parameters."""
        defaults = {
            "A": 1.0,
            "B": 1.0,
            "m": -1.0,
            "E": 0.0,
            "gamma": 0.0,
        }

        # Get current exploration state
        state_key = ("explore", self.exploration_depth, str(sorted(self.exploration_visited)))

        if state_key not in self.exploration_visited:
            # Generate plots at current neighbors
            candidates = []

            if self.exploration_depth == 0:
                # Start at defaults
                target = dict(defaults)
                for plot_type in self.viewer.display_plot_types:
                    rec = self.viewer._best_record_for_type(plot_type, target)
                    if rec is not None:
                        candidates.append((rec["record_id"], plot_type))
                self.exploration_visited.add(state_key)
            else:
                # Explore neighbors at current depth
                base_idx = {}
                for pname in sorted(self.viewer.param_ranges.keys()):
                    vals = self.viewer.param_ranges.get(pname, [])
                    if not vals:
                        continue
                    target_val = defaults.get(pname, vals[0])
                    idx = min(range(len(vals)), key=lambda i: abs(float(vals[i]) - target_val))
                    base_idx[pname] = idx

                # Generate neighbor targets
                neighbor_targets = []
                for pname in sorted(base_idx.keys()):
                    max_i = len(self.viewer.param_ranges.get(pname, [])) - 1
                    if max_i < 1:
                        continue
                    for delta in range(1, self.exploration_depth + 1):
                        for sign in (-1, 1):
                            ni = base_idx[pname] + sign * delta
                            if 0 <= ni <= max_i:
                                idx_copy = dict(base_idx)
                                idx_copy[pname] = ni
                                neighbor_targets.append(self.viewer._target_from_indices(idx_copy))

                # Sample and generate
                for target in neighbor_targets[:100]:  # Batch of 100 per depth level
                    for plot_type in self.viewer.display_plot_types:
                        rec = self.viewer._best_record_for_type(plot_type, target)
                        if rec is not None:
                            candidates.append((rec["record_id"], plot_type))

                self.exploration_visited.add(state_key)

            # Render candidates, skip already-generated
            to_render = []
            for record_id, plot_type in candidates:
                key = f"{record_id}:{plot_type}"
                if key not in self.generated_index:
                    to_render.append((record_id, plot_type))

            if to_render:
                self.viewer._render_pairs(to_render[:50])  # Batch render up to 50
                with self.lock:
                    for record_id, plot_type in to_render[:50]:
                        self.generated_index.add(f"{record_id}:{plot_type}")
                    self._save_generated_index()

            # Move to next depth after current batch
            if len(to_render) == 0 or len(candidates) < 20:
                self.exploration_depth += 1
                if self.exploration_depth > 5:
                    self.exploration_depth = 0
                    time.sleep(2)  # Pause before restarting exploration

    def stop(self):
        """Gracefully stop the worker."""
        self.running = False
        self._save_generated_index()


class LazyPlotViewer(QMainWindow):
    def __init__(self, source_path):
        super().__init__()
        self.source_path = Path(source_path)
        self.records_file, self.groups_file = self._resolve_packet_paths(self.source_path)
        self.packet_dir = self.records_file.parent
        self.cache_dir = self.packet_dir / "lazy_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.records = []
        self.rec_by_id = {}
        self.records_by_plot_type = {}
        self.current_params = {}
        self.param_ranges = {}
        self.available_plot_types = []
        self.display_plot_types = []
        self.sliders = {}
        self.labels = {}
        self.plot_cards = {}
        self.rendered_path_cache = {}
        self._manifest_cache = {}
        self.selection_cache = {}
        self.perturbation_values = []
        self.selected_perturbation = None
        self.range_mode = os.environ.get("SPECLOC_LAZY_RIBBON_RANGE", "fixed").lower()
        if self.range_mode not in ("fixed", "auto"):
            self.range_mode = "fixed"
        self.perturb_button = None
        self.range_button = None
        self.plots_container = None
        self.spectrum_ymin = None
        self.spectrum_ymax = None
        self._spectrum_ylim_version = 0
        self._current_target_E = None

        self.prefetch_enabled = True
        self.prefetch_radius = max(1, int(os.environ.get("SPECLOC_LAZY_PREFETCH_RADIUS", "1")))
        self.prefetch_limit = max(8, int(os.environ.get("SPECLOC_LAZY_PREFETCH_LIMIT", "48")))
        self.cache_limit_bytes = int(float(os.environ.get("SPECLOC_LAZY_CACHE_GB", "8")) * (1024 ** 3))

        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_display)
        self.update_timer.setSingleShot(True)

        self.prefetch_timer = QTimer()
        self.prefetch_timer.timeout.connect(self._prefetch_neighbors)
        self.prefetch_timer.setSingleShot(True)

        self._last_target_params = None

        self._load_records()
        default_prefetch = "0" if len(self.records) > 100000 else "1"
        self.prefetch_enabled = os.environ.get("SPECLOC_LAZY_PREFETCH", default_prefetch).lower() not in ("0", "false", "no")
        self._init_ui()

        # Start background plot generation worker
        self.plot_worker = PlotGenerationWorker(self)
        self.plot_worker.start()
        print("Background plot worker started")

    def _resolve_packet_paths(self, source):
        if source.is_file() and source.name == "lazy_packet_records.tsv":
            records = source
            groups = source.with_name("lazy_packet_groups.tsv")
            if not groups.exists():
                raise RuntimeError(f"Missing groups file: {groups}")
            return records, groups

        if source.is_dir():
            direct_records = source / "lazy_packet_records.tsv"
            direct_groups = source / "lazy_packet_groups.tsv"
            if direct_records.exists() and direct_groups.exists():
                return direct_records, direct_groups

            rec_matches = sorted(source.rglob("lazy_packet_records.tsv"), key=lambda p: p.stat().st_mtime)
            if rec_matches:
                records = rec_matches[-1]
                groups = records.with_name("lazy_packet_groups.tsv")
                if groups.exists():
                    return records, groups

        raise RuntimeError("Could not resolve lazy packet records/groups from source path")

    @staticmethod
    def _try_float(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return float("nan")

    def _ordered_plot_types(self, available_types):
        preferred_present = [p for p in PREFERRED_PLOT_ORDER if p in available_types]
        extras = sorted([p for p in available_types if p not in PREFERRED_PLOT_ORDER])
        return preferred_present + extras

    def _refresh_display_plot_types(self):
        available = set(self.available_plot_types)

        ribbon_types = []
        for base in RIBBON_BASE_TYPES:
            pt = f"{base}_{self.range_mode}"
            if pt in available and pt not in DISABLED_PLOT_TYPES:
                ribbon_types.append(pt)

        other_types = []
        for pt in self._ordered_plot_types(self.available_plot_types):
            if pt in DISABLED_PLOT_TYPES:
                continue
            if any(pt == f"{base}_fixed" or pt == f"{base}_auto" for base in RIBBON_BASE_TYPES):
                continue
            other_types.append(pt)

        self.display_plot_types = ribbon_types + other_types

    def _relevant_target_params(self, plot_type, target_params):
        base_keys = (
            "A", "B", "m", "B_y", "winding", "Lx_ribbon", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y",
        )

        def pick(keys):
            return {k: target_params[k] for k in keys if k in target_params}

        if plot_type in ("specloc_signature_vs_E", "specloc_loggap_vs_E"):
            return pick(base_keys + ("gamma", "W", "kappa"))
        if plot_type in ("specloc_gap_gamma_E", "specloc_sig_gamma_E"):
            return pick(base_keys + ("W", "kappa"))
        if plot_type == "specloc_spectrum_vs_gamma":
            return pick(base_keys + ("W", "kappa", "E"))
        if plot_type == "specloc_spectrum_vs_W":
            return pick(base_keys + ("gamma", "kappa", "E"))
        if plot_type == "specloc_spectrum_vs_kappa":
            return pick(base_keys + ("gamma", "W", "E"))
        if plot_type in ("dos", "ldos_target", "ldos_lowest"):
            keys = base_keys + ("gamma",)
            if plot_type == "ldos_target":
                keys = keys + ("E",)
            return pick(keys)
        return dict(target_params)

    def _load_records(self):
        numeric_keys = {
            "A", "B", "m", "B_y", "winding", "Lx_ribbon", "Lx_obc", "Ly_obc",
            "specloc_x", "specloc_y", "gamma", "W", "kappa", "E",
        }

        self.records = []
        with open(self.records_file, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rec = dict(row)
                rec["record_id"] = str(rec["record_id"])
                for k in numeric_keys:
                    if k in rec:
                        rec[k] = self._try_float(rec[k])
                self.records.append(rec)

        self.rec_by_id = {r["record_id"]: r for r in self.records}

        for k in sorted(numeric_keys):
            vals = sorted({r[k] for r in self.records if k in r and not math.isnan(r[k])})
            if len(vals) > 1:
                self.param_ranges[k] = vals

        self.available_plot_types = sorted({r["plot_type"] for r in self.records})
        self.perturbation_values = sorted({str(r.get("perturbation_type", "")) for r in self.records if str(r.get("perturbation_type", ""))})
        if self.selected_perturbation not in self.perturbation_values:
            if "symmetric" in self.perturbation_values:
                self.selected_perturbation = "symmetric"
            elif self.perturbation_values:
                self.selected_perturbation = self.perturbation_values[0]
            else:
                self.selected_perturbation = None
        self._refresh_display_plot_types()
        self.records_by_plot_type = {ptype: [] for ptype in self.available_plot_types}
        for rec in self.records:
            self.records_by_plot_type[rec["plot_type"]].append(rec)

        # Backward compatibility: allow gamma-vs-E heatmap cards to render even
        # when older packets do not explicitly contain these plot types.
        if self._has_specloc_group_records():
            for ptype in sorted(VIRTUAL_SPEClOC_PLOT_TYPES):
                if ptype not in self.available_plot_types:
                    self.available_plot_types.append(ptype)

        print(f"Loaded records: {self.records_file}")
        print(f"Loaded groups:  {self.groups_file}")
        print(f"Records: {len(self.records)}")
        print(f"Plot types: {len(self.available_plot_types)}")

    def _distance(self, rec, target_params):
        d = 0.0
        used = 0
        missing_penalty = 1e6
        for k, v in target_params.items():
            if k not in rec:
                continue
            rv = rec[k]
            if isinstance(rv, float) and math.isnan(rv):
                if isinstance(v, float) and not math.isnan(v):
                    # Prefer records that actually carry the requested axis value,
                    # e.g. gamma-resolved DOS/LDOS over legacy NaN-gamma rows.
                    d += missing_penalty
                continue
            d += (float(rv) - float(v)) ** 2
            used += 1
        if used == 0:
            return float("inf")
        return d

    @staticmethod
    def _is_finite_number(v):
        return isinstance(v, float) and not math.isnan(v)

    def _has_specloc_group_records(self):
        for rec in self.records:
            if rec.get("render_kind") == "specloc_group":
                return True
        return False

    def _fallback_candidates_for_virtual_type(self, plot_type):
        if plot_type not in VIRTUAL_SPEClOC_PLOT_TYPES:
            return []

        # Reuse any specloc-group record id; renderer can still resolve the
        # requested plot_type from the generated bundle manifest.
        candidates = []
        for anchor_type in ("specloc_signature_vs_E", "specloc_loggap_vs_E", "specloc_spectrum_vs_gamma"):
            candidates.extend(self.records_by_plot_type.get(anchor_type, []))

        if candidates:
            return candidates

        return [r for r in self.records if r.get("render_kind") == "specloc_group"]

    def _best_record_for_type(self, plot_type, target_params):
        relevant = self._relevant_target_params(plot_type, target_params)
        key = (plot_type, self.selected_perturbation, tuple(sorted(relevant.items())))
        cached = self.selection_cache.get(key)
        if cached is not None:
            return cached

        candidates = self.records_by_plot_type.get(plot_type, [])
        if not candidates:
            candidates = self._fallback_candidates_for_virtual_type(plot_type)

        if self.selected_perturbation is not None:
            candidates = [r for r in candidates if str(r.get("perturbation_type", "")) == self.selected_perturbation]

        # For DOS/LDOS, prefer gamma-resolved rows when available.
        if plot_type in ("dos", "ldos_target", "ldos_lowest") and "gamma" in relevant:
            finite_gamma = [r for r in candidates if self._is_finite_number(r.get("gamma"))]
            if finite_gamma:
                candidates = finite_gamma

        # For target-LDOS, prefer rows with finite target energy when available.
        if plot_type == "ldos_target" and "E" in relevant:
            finite_e = [r for r in candidates if self._is_finite_number(r.get("E"))]
            if finite_e:
                candidates = finite_e

        if not candidates:
            return None
        best = min(candidates, key=lambda r: self._distance(r, relevant))
        self.selection_cache[key] = best
        return best

    @staticmethod
    def _path_exists_and_nonempty(path):
        if not path:
            return False
        try:
            return Path(path).is_file() and Path(path).stat().st_size > 0
        except OSError:
            return False

    def _render_variant(self, plot_type, target_params=None):
        if target_params is None:
            target_params = self.current_params

        gamma = target_params.get("gamma", float("nan"))
        gamma_tag = "nan" if math.isnan(float(gamma)) else f"{float(gamma):.12g}"

        if plot_type == "specloc_spectrum_vs_gamma":
            return f"{self._spectrum_ylim_version}:{self.spectrum_ymin}:{self.spectrum_ymax}"
        if plot_type in ("specloc_signature_vs_E", "specloc_loggap_vs_E"):
            return f"Eref:{self._current_target_E}"
        if plot_type in ("dos", "ldos_target"):
            return f"gamma:{gamma_tag}"
        return "base"

    def _read_manifest_rows(self, image_path):
        p = Path(image_path)
        try:
            bundle_dir = p.parents[1]
        except IndexError:
            return []
        manifest = bundle_dir / "manifest.tsv"
        if not manifest.is_file():
            return []

        try:
            mtime = manifest.stat().st_mtime
        except OSError:
            return []

        key = str(manifest)
        cached = self._manifest_cache.get(key)
        if cached is not None and cached.get("mtime") == mtime:
            return cached.get("rows", [])

        rows = []
        try:
            with open(manifest, "r", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    rows.append({
                        "plot_type": str(row.get("plot_type", "")),
                        "gamma": self._try_float(row.get("gamma", "")),
                        "W": self._try_float(row.get("W", "")),
                        "kappa": self._try_float(row.get("kappa", "")),
                        "E": self._try_float(row.get("E", "")),
                        "plot_path": str(row.get("plot_path", "")),
                    })
        except OSError:
            return []

        self._manifest_cache[key] = {"mtime": mtime, "rows": rows}
        return rows

    def _distance_manifest_row(self, row, target_params):
        dist = 0.0
        used = 0
        for k in ("gamma", "W", "kappa", "E"):
            rv = row.get(k, float("nan"))
            tv = target_params.get(k, float("nan"))
            if not math.isnan(float(rv)) and not math.isnan(float(tv)):
                dist += (float(rv) - float(tv)) ** 2
                used += 1
        return dist if used > 0 else float("inf")

    def _gamma_filtered_path(self, image_path, plot_type, target_params):
        if plot_type not in ("dos", "ldos_target"):
            return image_path

        rows = self._read_manifest_rows(image_path)
        if not rows:
            return image_path

        candidates = [r for r in rows if r.get("plot_type") == plot_type and self._path_exists_and_nonempty(r.get("plot_path", ""))]
        if not candidates:
            return image_path

        req_gamma = target_params.get("gamma", float("nan"))
        if not math.isnan(float(req_gamma)):
            finite_gamma = [r for r in candidates if not math.isnan(float(r.get("gamma", float("nan"))))]
            if finite_gamma:
                candidates = finite_gamma

        best = min(candidates, key=lambda r: self._distance_manifest_row(r, target_params))
        return best.get("plot_path", image_path)

    def _render_pairs(self, pairs, target_params=None):
        if not pairs:
            return {}

        if target_params is None:
            target_params = dict(self.current_params)

        # Deduplicate while preserving order.
        seen = set()
        deduped = []
        for p in pairs:
            if p in seen:
                continue
            seen.add(p)
            deduped.append(p)
        pairs = deduped

        out = {}
        missing = []
        missing_has_spectrum = False
        for pair in pairs:
            rid, ptype = pair
            key = (rid, ptype, self._render_variant(ptype, target_params))
            cached_path = self.rendered_path_cache.get(key, "")
            if self._path_exists_and_nonempty(cached_path):
                out[pair] = self._gamma_filtered_path(cached_path, ptype, target_params)
            else:
                missing.append(pair)
                if ptype == "specloc_spectrum_vs_gamma":
                    missing_has_spectrum = True

        if not missing:
            return out

        # Build request file for batch rendering.
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False, dir=str(self.cache_dir)) as tf:
            tf.write("record_id\tplot_type\n")
            for rid, ptype in missing:
                tf.write(f"{rid}\t{ptype}\n")
            req_file = Path(tf.name)

        script_path = Path(__file__).resolve().parents[1] / "simulations" / "data_processing" / "lazy_render.jl"
        cmd = [
            "julia",
            str(script_path),
            str(self.records_file),
            str(self.groups_file),
            str(self.cache_dir),
            str(req_file),
        ]

        env = dict(os.environ)
        env.setdefault("GKSwstype", "png")
        env.setdefault("GKS_WSTYPE", env["GKSwstype"])
        env.setdefault("GKS_NO_GUI", "1")
        env.setdefault("SPECLOC_SUPPRESS_GR_WARNINGS", "true")
        if self._current_target_E is not None:
            env["SPECLOC_SIGGAP_E_REF"] = str(self._current_target_E)
        if self.spectrum_ymin is not None and self.spectrum_ymax is not None:
            env["SPECLOC_SPECTRUM_YMIN"] = str(self.spectrum_ymin)
            env["SPECLOC_SPECTRUM_YMAX"] = str(self.spectrum_ymax)
        if missing_has_spectrum:
            env["SPECLOC_LAZY_FORCE_REBUILD"] = "true"

        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
            if proc.returncode != 0:
                print(proc.stderr)
                return out
            for line in proc.stdout.splitlines():
                parts = line.split("\t")
                if len(parts) != 3:
                    continue
                rid, ptype, path = parts
                pair = (rid, ptype)
                key = (rid, ptype, self._render_variant(ptype, target_params))
                self.rendered_path_cache[key] = path
                out[pair] = self._gamma_filtered_path(path, ptype, target_params)

            self._enforce_cache_limit()
            # Drop stale cache entries after eviction.
            stale = [k for k, p in self.rendered_path_cache.items() if not self._path_exists_and_nonempty(p)]
            for k in stale:
                self.rendered_path_cache.pop(k, None)
            return out
        finally:
            try:
                req_file.unlink(missing_ok=True)
            except Exception:
                pass

    def _cache_files(self):
        files = []
        if not self.cache_dir.exists():
            return files
        for p in self.cache_dir.rglob("*"):
            if not p.is_file():
                continue
            try:
                st = p.stat()
            except OSError:
                continue
            files.append((p, st.st_mtime, st.st_size))
        return files

    def _enforce_cache_limit(self):
        files = self._cache_files()
        total = sum(s for _, _, s in files)
        if total <= self.cache_limit_bytes:
            return

        files.sort(key=lambda x: x[1])
        for p, _, s in files:
            if total <= self.cache_limit_bytes:
                break
            try:
                p.unlink()
                total -= s
            except OSError:
                pass

        # Cleanup empty directories.
        for d in sorted(self.cache_dir.rglob("*"), reverse=True):
            if d.is_dir():
                try:
                    d.rmdir()
                except OSError:
                    pass

    def _current_slider_indices(self):
        idx = {}
        for pname, slider in self.sliders.items():
            idx[pname] = slider.value()
        return idx

    def _target_from_indices(self, idx_map):
        t = {}
        for pname, i in idx_map.items():
            vals = self.param_ranges.get(pname, [])
            if not vals:
                continue
            i2 = max(0, min(i, len(vals) - 1))
            t[pname] = vals[i2]
        return t

    def _prefetch_neighbors(self):
        if not self.prefetch_enabled or self._last_target_params is None:
            return

        base_idx = self._current_slider_indices()
        neighbor_targets = []

        for pname in sorted(base_idx.keys()):
            max_i = len(self.param_ranges.get(pname, [])) - 1
            if max_i < 1:
                continue
            for delta in range(1, self.prefetch_radius + 1):
                for sign in (-1, 1):
                    ni = base_idx[pname] + sign * delta
                    if ni < 0 or ni > max_i:
                        continue
                    idx_copy = dict(base_idx)
                    idx_copy[pname] = ni
                    neighbor_targets.append(self._target_from_indices(idx_copy))

        pairs = []
        for t in neighbor_targets:
            for plot_type in self.display_plot_types:
                rec = self._best_record_for_type(plot_type, t)
                if rec is None:
                    continue
                pairs.append((rec["record_id"], plot_type))
                if len(pairs) >= self.prefetch_limit:
                    break
            if len(pairs) >= self.prefetch_limit:
                break

        if pairs:
            self._render_pairs(pairs)

    def _init_ui(self):
        self.setWindowTitle("SpecLoc Lazy Plot Viewer")
        self.setGeometry(80, 80, 1800, 1040)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout()

        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        control_layout.setSpacing(6)
        control_layout.setAlignment(Qt.AlignTop)
        control_panel.setMinimumWidth(290)

        for param_name in sorted(self.param_ranges.keys()):
            if not self.param_ranges[param_name]:
                continue
            param_layout = QVBoxLayout()

            label = QLabel(f"{param_name}: {self.param_ranges[param_name][0]:.3f}")
            self.labels[param_name] = label
            param_layout.addWidget(label)

            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(0)
            slider.setMaximum(len(self.param_ranges[param_name]) - 1)
            slider.setValue(0)
            slider.setTickPosition(QSlider.TicksBelow)
            slider.setTickInterval(max(1, len(self.param_ranges[param_name]) // 10))
            slider.valueChanged.connect(lambda value, p=param_name: self._on_slider_move(p, value))

            self.sliders[param_name] = slider
            param_layout.addWidget(slider)
            control_layout.addLayout(param_layout)

        self._set_default_slider_values()

        reset_button = QPushButton("Reset")
        reset_button.clicked.connect(self._reset_params)
        control_layout.addWidget(reset_button)

        self.perturb_button = QPushButton()
        self.perturb_button.clicked.connect(self._toggle_perturbation)
        control_layout.addWidget(self.perturb_button)

        self.range_button = QPushButton()
        self.range_button.clicked.connect(self._toggle_range_mode)
        control_layout.addWidget(self.range_button)

        ylim_row = QHBoxLayout()
        self.ymin_input = QLineEdit()
        self.ymin_input.setPlaceholderText("spectrum y min")
        self.ymin_input.setFixedWidth(120)
        self.ymax_input = QLineEdit()
        self.ymax_input.setPlaceholderText("spectrum y max")
        self.ymax_input.setFixedWidth(120)
        ylim_apply = QPushButton("Apply y-lims")
        ylim_apply.clicked.connect(self._apply_spectrum_ylims)
        ylim_row.addWidget(self.ymin_input)
        ylim_row.addWidget(self.ymax_input)
        ylim_row.addWidget(ylim_apply)
        control_layout.addLayout(ylim_row)

        self._update_toggle_labels()

        control_layout.addStretch()

        control_scroll = QScrollArea()
        control_scroll.setWidgetResizable(True)
        control_scroll.setMinimumWidth(280)
        control_scroll.setMaximumWidth(320)
        control_scroll.setWidget(control_panel)

        self.plots_container = QWidget()
        self.image_layout = QGridLayout(self.plots_container)
        self.image_layout.setSpacing(6)
        self.image_layout.setContentsMargins(4, 4, 4, 4)
        self._rebuild_plot_grid()

        plots_scroll = QScrollArea()
        plots_scroll.setWidgetResizable(True)
        plots_scroll.setWidget(self.plots_container)

        main_layout.addWidget(control_scroll, 0)
        main_layout.addWidget(plots_scroll, 1)
        main_layout.setSpacing(8)
        main_layout.setContentsMargins(6, 6, 6, 6)

        central_widget.setLayout(main_layout)
        self.update_display()

    def _update_toggle_labels(self):
        if self.perturb_button is not None:
            if self.selected_perturbation is None:
                self.perturb_button.setText("Distortion: all")
            else:
                self.perturb_button.setText(f"Distortion: {self.selected_perturbation}")
        if self.range_button is not None:
            self.range_button.setText(f"Ribbon range: {self.range_mode}")

    def _toggle_perturbation(self):
        if not self.perturbation_values:
            return
        if self.selected_perturbation not in self.perturbation_values:
            self.selected_perturbation = self.perturbation_values[0]
        else:
            i = self.perturbation_values.index(self.selected_perturbation)
            self.selected_perturbation = self.perturbation_values[(i + 1) % len(self.perturbation_values)]
        self.selection_cache.clear()
        self._update_toggle_labels()
        self.update_display()

    def _toggle_range_mode(self):
        self.range_mode = "auto" if self.range_mode == "fixed" else "fixed"
        self._refresh_display_plot_types()
        self.selection_cache.clear()
        self._rebuild_plot_grid()
        self._update_toggle_labels()
        self.update_display()

    def _apply_spectrum_ylims(self):
        tmin = self.ymin_input.text().strip() if self.ymin_input is not None else ""
        tmax = self.ymax_input.text().strip() if self.ymax_input is not None else ""

        if tmin == "" and tmax == "":
            self.spectrum_ymin = None
            self.spectrum_ymax = None
        else:
            try:
                ymin = float(tmin)
                ymax = float(tmax)
                if ymin >= ymax:
                    raise ValueError("ymin must be less than ymax")
                self.spectrum_ymin = ymin
                self.spectrum_ymax = ymax
            except ValueError:
                return

        self._spectrum_ylim_version += 1
        self.selection_cache.clear()
        self.rendered_path_cache.clear()
        self.update_display()

    def _rebuild_plot_grid(self):
        while self.image_layout.count():
            item = self.image_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

        self.plot_cards = {}
        for i, plot_type in enumerate(self.display_plot_types):
            meta = QLabel("Record: pending")
            meta.setAlignment(Qt.AlignCenter)
            meta.setWordWrap(True)
            meta.setStyleSheet("color: #444; font-size: 11px;")

            image = QLabel("Pending")
            image.setAlignment(Qt.AlignCenter)
            image.setStyleSheet("border: 1px solid #999; background: white;")
            image.setMinimumSize(420, 320)
            image.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

            card = QWidget()
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(3, 3, 3, 3)
            card_layout.setSpacing(4)
            card_layout.addWidget(meta)
            card_layout.addWidget(image)

            self.plot_cards[plot_type] = {"meta": meta, "image": image, "pixmap": None}
            r = i // 3
            c = i % 3
            self.image_layout.addWidget(card, r, c)

    def _format_record_metadata(self, rec, plot_type):
        if rec is None:
            return "Record: none"

        rid = str(rec.get("record_id", "?"))
        keys = ["A", "B", "m", "gamma", "W", "kappa", "E"]
        parts = []
        for k in keys:
            v = rec.get(k, None)
            if isinstance(v, float):
                if math.isnan(v):
                    continue
                parts.append(f"{k}={v:.3f}")
            elif v is not None and v != "":
                parts.append(f"{k}={v}")

        pt = rec.get("perturbation_type", "")
        dt = rec.get("disorder_type", "")
        if pt:
            parts.append(f"pt={pt}")
        if dt:
            parts.append(f"dt={dt}")

        return f"Record {rid} ({plot_type})\n" + ", ".join(parts)

    def _set_default_slider_values(self):
        defaults = {
            "A": 1.0,
            "B": 1.0,
            "m": -1.0,
            "E": 0.0,
            "gamma": 0.0,
        }
        for pname, target in defaults.items():
            if pname not in self.sliders:
                continue
            vals = self.param_ranges.get(pname, [])
            if not vals:
                continue
            idx = min(range(len(vals)), key=lambda i: abs(float(vals[i]) - target))
            slider = self.sliders[pname]
            slider.blockSignals(True)
            slider.setValue(idx)
            slider.blockSignals(False)
            val = vals[idx]
            self.current_params[pname] = val
            self.labels[pname].setText(f"{pname}: {val:.3f}")

    def _on_slider_move(self, param_name, value):
        param_value = self.param_ranges[param_name][value]
        self.current_params[param_name] = param_value
        self.labels[param_name].setText(f"{param_name}: {param_value:.3f}")
        self.update_timer.start(180)

    def _display_image(self, label, image_path, plot_type, target_params):
        if not image_path:
            label.setText("No rendered path")
            label.setPixmap(QPixmap())
            return

        pixmap = QPixmap(image_path)
        if pixmap.isNull():
            label.setText(f"Error loading image:\n{image_path}")
            label.setPixmap(QPixmap())
            return

        target_rect = label.contentsRect()
        scaled = pixmap.scaled(target_rect.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
        label.setPixmap(scaled)
        label.setText("")

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.update_timer.start(120)

    def update_display(self):
        if not self.param_ranges:
            return

        target_params = {}
        for param_name in self.param_ranges:
            if param_name in self.current_params:
                target_params[param_name] = self.current_params[param_name]
            else:
                target_params[param_name] = self.param_ranges[param_name][0]
                self.current_params[param_name] = target_params[param_name]

            self._current_target_E = target_params.get("E", None)

        param_str = ", ".join(f"{k}={v:.3f}" for k, v in sorted(target_params.items()))
        ptxt = self.selected_perturbation if self.selected_perturbation is not None else "all"
        self.setWindowTitle(f"SpecLoc Lazy Plot Viewer - distortion={ptxt}, range={self.range_mode}, {param_str}")

        if len(self.selection_cache) > 20000:
            self.selection_cache.clear()

        selected = []
        selected_by_type = {}
        for plot_type in self.display_plot_types:
            rec = self._best_record_for_type(plot_type, target_params)
            if rec is None:
                continue
            selected.append((rec["record_id"], plot_type))
            selected_by_type[plot_type] = rec

        rendered_map = self._render_pairs(selected, target_params=target_params)

        for plot_type, card in self.plot_cards.items():
            rec = selected_by_type.get(plot_type)
            if rec is None:
                card["meta"].setText("Record: no matching record")
                card["image"].setText("No matching record")
                card["image"].setPixmap(QPixmap())
                continue

            card["meta"].setText(self._format_record_metadata(rec, plot_type))
            image_path = rendered_map.get((rec["record_id"], plot_type), "")
            self._display_image(card["image"], image_path, plot_type, target_params)

        self._last_target_params = dict(target_params)
        if self.prefetch_enabled:
            # Start prefetch shortly after visible plots are updated.
            self.prefetch_timer.start(120)

    def _reset_params(self):
        for param_name in self.sliders:
            self.sliders[param_name].setValue(0)
            self.current_params[param_name] = self.param_ranges[param_name][0]
            self.labels[param_name].setText(f"{param_name}: {self.param_ranges[param_name][0]:.3f}")
        self.update_display()


def main():
    if len(sys.argv) > 1:
        source = Path(sys.argv[1])
    else:
        source = Path("simulations/data_processing/plots")
        if not source.exists():
            print("Usage: python3 useful_scripts/qwz_lazy_plt_viewer.py <lazy_packet dir or lazy_packet_records.tsv>")
            sys.exit(1)

    app = QApplication(sys.argv[:1])
    viewer = LazyPlotViewer(source)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
