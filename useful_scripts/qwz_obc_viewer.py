#!/usr/bin/env python3
"""
OBC Density of States Viewer

Displays DOS (density of states) and LDOS (local density of states) for OBC systems.
Auto-detects gamma-resolved vs legacy records per plot type.

Usage:
    python3 useful_scripts/qwz_obc_viewer.py <lazy_packet_dir>
"""

import csv
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (
    QApplication, QGridLayout, QHBoxLayout, QLabel, QMainWindow,
    QPushButton, QScrollArea, QSizePolicy, QSlider, QVBoxLayout, QWidget,
)


class OBCViewer(QMainWindow):
    def __init__(self, packet_dir):
        super().__init__()
        self.packet_dir = Path(packet_dir)
        self.records_file = self.packet_dir / "lazy_packet_records.tsv"
        self.groups_file = self.packet_dir / "lazy_packet_groups.tsv"
        self.cache_dir = self.packet_dir / "lazy_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.records = []
        self.param_ranges = {}
        self.current_params = {}
        self.sliders = {}
        self.labels = {}
        self.plot_cards = {}
        self.rendered_cache = {}
        self.perturbation_values = []
        self.selected_perturbation = "symmetric"
        self.has_gamma = {}  # Track per plot_type whether gamma-resolved records exist

        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_display)
        self.update_timer.setSingleShot(True)

        self._load_records()
        self._init_ui()

    def _load_records(self):
        numeric_keys = {"A", "B", "m", "B_y", "winding", "Lx_obc", "Ly_obc", "gamma", "E"}
        
        with open(self.records_file, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("render_kind") != "case":
                    continue
                ptype = row.get("plot_type", "")
                if ptype not in ("dos", "ldos_target", "ldos_lowest"):
                    continue
                
                rec = dict(row)
                for k in numeric_keys:
                    if k in rec:
                        try:
                            rec[k] = float(rec[k])
                        except (TypeError, ValueError):
                            rec[k] = float("nan")
                self.records.append(rec)

        # Detect gamma availability per plot type
        for ptype in ("dos", "ldos_target", "ldos_lowest"):
            has_finite = any(not math.isnan(r.get("gamma", float("nan"))) for r in self.records if r.get("plot_type") == ptype)
            self.has_gamma[ptype] = has_finite

        for k in sorted(numeric_keys):
            vals = sorted({r[k] for r in self.records if k in r and not math.isnan(r[k])})
            if len(vals) > 1:
                self.param_ranges[k] = vals

        self.perturbation_values = sorted({str(r.get("perturbation_type", "")) for r in self.records if str(r.get("perturbation_type", ""))})
        if "symmetric" not in self.perturbation_values and self.perturbation_values:
            self.selected_perturbation = self.perturbation_values[0]

        print(f"Loaded {len(self.records)} OBC records")
        print(f"Gamma-resolved: dos={self.has_gamma.get('dos')}, ldos_target={self.has_gamma.get('ldos_target')}, ldos_lowest={self.has_gamma.get('ldos_lowest')}")

    def _init_ui(self):
        self.setWindowTitle("OBC Density of States Viewer")
        self.setGeometry(80, 80, 1400, 600)

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout()

        # Control panel
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        control_layout.setAlignment(Qt.AlignTop)
        control_panel.setMinimumWidth(280)

        for param in sorted(self.param_ranges.keys()):
            if not self.param_ranges[param]:
                continue
            label = QLabel(f"{param}: {self.param_ranges[param][0]:.3f}")
            self.labels[param] = label
            control_layout.addWidget(label)

            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(0)
            slider.setMaximum(len(self.param_ranges[param]) - 1)
            slider.setValue(0)
            slider.valueChanged.connect(lambda v, p=param: self._on_slider(p, v))
            self.sliders[param] = slider
            control_layout.addWidget(slider)

        self._set_defaults()

        reset_btn = QPushButton("Reset")
        reset_btn.clicked.connect(self._reset)
        control_layout.addWidget(reset_btn)

        self.perturb_btn = QPushButton()
        self.perturb_btn.clicked.connect(self._toggle_perturbation)
        control_layout.addWidget(self.perturb_btn)

        self._update_button_labels()
        control_layout.addStretch()

        control_scroll = QScrollArea()
        control_scroll.setWidgetResizable(True)
        control_scroll.setWidget(control_panel)
        control_scroll.setMaximumWidth(320)

        # Plot grid (1 row × 3 columns)
        plots_widget = QWidget()
        grid = QGridLayout(plots_widget)
        grid.setSpacing(8)

        plot_types = [
            ("dos", "Total DOS vs Energy"),
            ("ldos_target", "LDOS at Target Energy"),
            ("ldos_lowest", "LDOS at Lowest |E| States"),
        ]

        for i, (ptype, title) in enumerate(plot_types):
            card = self._make_card(title)
            self.plot_cards[ptype] = card
            grid.addWidget(card["widget"], 0, i)

        plots_scroll = QScrollArea()
        plots_scroll.setWidgetResizable(True)
        plots_scroll.setWidget(plots_widget)

        main_layout.addWidget(control_scroll, 0)
        main_layout.addWidget(plots_scroll, 1)
        central.setLayout(main_layout)

        self.update_display()

    def _make_card(self, title):
        title_lbl = QLabel(title)
        title_lbl.setAlignment(Qt.AlignCenter)
        title_lbl.setStyleSheet("font-weight: bold;")

        meta_lbl = QLabel("Record: pending")
        meta_lbl.setAlignment(Qt.AlignCenter)
        meta_lbl.setWordWrap(True)
        meta_lbl.setStyleSheet("color: #444; font-size: 10px;")

        img_lbl = QLabel("Pending")
        img_lbl.setAlignment(Qt.AlignCenter)
        img_lbl.setStyleSheet("border: 1px solid #999; background: white;")
        img_lbl.setMinimumSize(380, 340)
        img_lbl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(3, 3, 3, 3)
        layout.addWidget(title_lbl)
        layout.addWidget(meta_lbl)
        layout.addWidget(img_lbl)

        return {"widget": widget, "title": title_lbl, "meta": meta_lbl, "image": img_lbl}

    def _set_defaults(self):
        defaults = {"A": 1.0, "B": 1.0, "m": -1.0, "gamma": 0.0, "E": 0.0}
        for p, target in defaults.items():
            if p not in self.sliders:
                continue
            vals = self.param_ranges[p]
            idx = min(range(len(vals)), key=lambda i: abs(vals[i] - target))
            self.sliders[p].blockSignals(True)
            self.sliders[p].setValue(idx)
            self.sliders[p].blockSignals(False)
            self.current_params[p] = vals[idx]
            self.labels[p].setText(f"{p}: {vals[idx]:.3f}")

    def _on_slider(self, param, idx):
        val = self.param_ranges[param][idx]
        self.current_params[param] = val
        self.labels[param].setText(f"{param}: {val:.3f}")
        self.update_timer.start(150)

    def _reset(self):
        for p in self.sliders:
            self.sliders[p].setValue(0)
        self.update_display()

    def _toggle_perturbation(self):
        if not self.perturbation_values:
            return
        if self.selected_perturbation not in self.perturbation_values:
            self.selected_perturbation = self.perturbation_values[0]
        else:
            i = self.perturbation_values.index(self.selected_perturbation)
            self.selected_perturbation = self.perturbation_values[(i + 1) % len(self.perturbation_values)]
        self._update_button_labels()
        self.rendered_cache.clear()
        self.update_display()

    def _update_button_labels(self):
        self.perturb_btn.setText(f"Distortion: {self.selected_perturbation}")

    def _best_record(self, plot_type):
        target = dict(self.current_params)
        candidates = [r for r in self.records 
                     if r.get("plot_type") == plot_type 
                     and str(r.get("perturbation_type", "")) == self.selected_perturbation]
        
        # Filter for gamma-resolved if available
        if self.has_gamma.get(plot_type) and "gamma" in target:
            gamma_candidates = [c for c in candidates if not math.isnan(c.get("gamma", float("nan")))]
            if gamma_candidates:
                candidates = gamma_candidates

        # For ldos_target, prefer finite E if available
        if plot_type == "ldos_target" and "E" in target:
            e_candidates = [c for c in candidates if not math.isnan(c.get("E", float("nan")))]
            if e_candidates:
                candidates = e_candidates

        if not candidates:
            return None

        def distance(rec):
            d = 0.0
            count = 0
            for k, v in target.items():
                rv = rec.get(k)
                if rv is None or math.isnan(rv):
                    continue
                d += (rv - v) ** 2
                count += 1
            return d if count > 0 else float("inf")

        return min(candidates, key=distance)

    def _render_images(self, records_needed):
        pairs = [(r["record_id"], r["plot_type"]) for r in records_needed]
        missing = [p for p in pairs if p not in self.rendered_cache]
        
        if not missing:
            return

        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            tf.write("record_id\tplot_type\n")
            for rid, ptype in missing:
                tf.write(f"{rid}\t{ptype}\n")
            req_file = tf.name

        script = Path(__file__).resolve().parents[1] / "simulations" / "data_processing" / "lazy_render.jl"
        cmd = ["julia", str(script), str(self.records_file), str(self.groups_file), str(self.cache_dir), req_file]
        
        env = dict(os.environ)
        env.setdefault("GKSwstype", "png")
        env["SPECLOC_SUPPRESS_GR_WARNINGS"] = "true"

        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
            if proc.returncode == 0:
                for line in proc.stdout.splitlines():
                    parts = line.split("\t")
                    if len(parts) == 3:
                        rid, ptype, path = parts
                        self.rendered_cache[(rid, ptype)] = path
        finally:
            Path(req_file).unlink(missing_ok=True)

    def update_display(self):
        target_str = ", ".join(f"{k}={v:.3f}" for k, v in sorted(self.current_params.items()))
        self.setWindowTitle(f"OBC Viewer - {self.selected_perturbation} - {target_str}")

        selected = {}
        for ptype in self.plot_cards.keys():
            rec = self._best_record(ptype)
            if rec:
                selected[ptype] = rec

        self._render_images(list(selected.values()))

        for ptype, card in self.plot_cards.items():
            rec = selected.get(ptype)
            if not rec:
                card["meta"].setText("Record: no match")
                card["image"].setText("No matching record")
                card["image"].setPixmap(QPixmap())
                continue

            rid = rec["record_id"]
            path = self.rendered_cache.get((rid, ptype), "")
            
            meta_parts = [f"ID={rid}"]
            for k in ["A", "B", "m", "gamma", "E"]:
                v = rec.get(k)
                if v is not None and not math.isnan(v):
                    meta_parts.append(f"{k}={v:.3f}")
            gamma_tag = " [gamma-resolved]" if self.has_gamma.get(ptype) else " [legacy]"
            card["meta"].setText(", ".join(meta_parts) + gamma_tag)

            if path and Path(path).is_file():
                pixmap = QPixmap(path)
                if not pixmap.isNull():
                    scaled = pixmap.scaled(card["image"].size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
                    card["image"].setPixmap(scaled)
                    card["image"].setText("")
                else:
                    card["image"].setText("Image load error")
                    card["image"].setPixmap(QPixmap())
            else:
                card["image"].setText("Rendering...")
                card["image"].setPixmap(QPixmap())


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 qwz_obc_viewer.py <lazy_packet_dir>")
        sys.exit(1)

    packet = Path(sys.argv[1])
    if not packet.exists():
        print(f"Error: {packet} does not exist")
        sys.exit(1)

    app = QApplication(sys.argv[:1])
    viewer = OBCViewer(packet)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
