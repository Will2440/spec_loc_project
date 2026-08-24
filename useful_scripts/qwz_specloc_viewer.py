#!/usr/bin/env python3
"""
Spectral Localiser Viewer

Displays aggregated spectral localiser analysis across multiple cases.
Handles signature/gap vs E, spectrum vs gamma/W/kappa, and various 2D heatmaps.

Usage:
    python3 useful_scripts/qwz_specloc_viewer.py <lazy_packet_dir>
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
    QApplication, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QMainWindow,
    QPushButton, QScrollArea, QSizePolicy, QSlider, QVBoxLayout, QWidget,
)


# Core plot types always shown
CORE_PLOT_TYPES = [
    ("specloc_signature_vs_E", "Signature vs E"),
    ("specloc_loggap_vs_E", "log(gap) vs E"),
    ("specloc_spectrum_vs_gamma", "Spectrum vs gamma"),
    ("specloc_gap_gamma_E", "Gap: gamma×E"),
    ("specloc_sig_gamma_E", "Signature: gamma×E"),
]

# Optional extended plot types (off by default for large packets)
EXTENDED_PLOT_TYPES = [
    ("specloc_gap_gamma_W", "Gap: gamma×W"),
    ("specloc_sig_gamma_W", "Signature: gamma×W"),
    ("specloc_gap_gamma_kappa", "Gap: gamma×kappa"),
    ("specloc_sig_gamma_kappa", "Signature: gamma×kappa"),
    ("specloc_gap_W_kappa", "Gap: W×kappa"),
    ("specloc_sig_W_kappa", "Signature: W×kappa"),
    ("specloc_spectrum_vs_W", "Spectrum vs W"),
    ("specloc_spectrum_vs_kappa", "Spectrum vs kappa"),
]


class SpeclocViewer(QMainWindow):
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
        self.show_extended = False
        self.spectrum_ymin = None
        self.spectrum_ymax = None
        self._spectrum_ylim_version = 0

        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_display)
        self.update_timer.setSingleShot(True)

        self._load_records()
        self._init_ui()

    def _load_records(self):
        numeric_keys = {"A", "B", "m", "B_y", "winding", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y", "gamma", "W", "kappa", "E"}
        
        with open(self.records_file, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("render_kind") != "specloc_group":
                    continue
                
                rec = dict(row)
                for k in numeric_keys:
                    if k in rec:
                        try:
                            rec[k] = float(rec[k])
                        except (TypeError, ValueError):
                            rec[k] = float("nan")
                self.records.append(rec)

        for k in sorted(numeric_keys):
            vals = sorted({r[k] for r in self.records if k in r and not math.isnan(r[k])})
            if len(vals) > 1:
                self.param_ranges[k] = vals

        self.perturbation_values = sorted({str(r.get("perturbation_type", "")) for r in self.records if str(r.get("perturbation_type", ""))})
        if "symmetric" not in self.perturbation_values and self.perturbation_values:
            self.selected_perturbation = self.perturbation_values[0]

        # Auto-disable extended plots for large packets
        self.show_extended = len(self.records) < 5000

        print(f"Loaded {len(self.records)} specloc group records")
        print(f"Parameters: {list(self.param_ranges.keys())}")
        print(f"Extended plots: {'enabled' if self.show_extended else 'disabled (large packet)'}")

    def _init_ui(self):
        self.setWindowTitle("Spectral Localiser Viewer")
        self.setGeometry(100, 100, 1800, 1000)

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

        self.extended_btn = QPushButton()
        self.extended_btn.clicked.connect(self._toggle_extended)
        control_layout.addWidget(self.extended_btn)

        # Spectrum y-limits
        ylim_label = QLabel("Spectrum Y-limits:")
        control_layout.addWidget(ylim_label)
        ylim_row = QHBoxLayout()
        self.ymin_input = QLineEdit()
        self.ymin_input.setPlaceholderText("ymin")
        self.ymin_input.setFixedWidth(80)
        self.ymax_input = QLineEdit()
        self.ymax_input.setPlaceholderText("ymax")
        self.ymax_input.setFixedWidth(80)
        ylim_apply = QPushButton("Apply")
        ylim_apply.clicked.connect(self._apply_ylims)
        ylim_row.addWidget(self.ymin_input)
        ylim_row.addWidget(self.ymax_input)
        ylim_row.addWidget(ylim_apply)
        control_layout.addLayout(ylim_row)

        self._update_button_labels()
        control_layout.addStretch()

        control_scroll = QScrollArea()
        control_scroll.setWidgetResizable(True)
        control_scroll.setWidget(control_panel)
        control_scroll.setMaximumWidth(320)

        # Plot grid container
        self.plots_container = QWidget()
        self.plots_layout = QGridLayout(self.plots_container)
        self.plots_layout.setSpacing(6)
        self._rebuild_plot_grid()

        plots_scroll = QScrollArea()
        plots_scroll.setWidgetResizable(True)
        plots_scroll.setWidget(self.plots_container)

        main_layout.addWidget(control_scroll, 0)
        main_layout.addWidget(plots_scroll, 1)
        central.setLayout(main_layout)

        self.update_display()

    def _rebuild_plot_grid(self):
        # Clear existing
        while self.plots_layout.count():
            item = self.plots_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        self.plot_cards = {}
        plot_list = list(CORE_PLOT_TYPES)
        if self.show_extended:
            plot_list.extend(EXTENDED_PLOT_TYPES)

        for i, (ptype, title) in enumerate(plot_list):
            card = self._make_card(title)
            self.plot_cards[ptype] = card
            r, c = divmod(i, 3)
            self.plots_layout.addWidget(card["widget"], r, c)

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
        img_lbl.setMinimumSize(380, 280)
        img_lbl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(3, 3, 3, 3)
        layout.addWidget(title_lbl)
        layout.addWidget(meta_lbl)
        layout.addWidget(img_lbl)

        return {"widget": widget, "title": title_lbl, "meta": meta_lbl, "image": img_lbl}

    def _set_defaults(self):
        defaults = {"A": 1.0, "B": 1.0, "m": -1.0, "gamma": 0.0, "W": 0.0, "kappa": 0.0, "E": 0.0}
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

    def _toggle_extended(self):
        self.show_extended = not self.show_extended
        self._update_button_labels()
        self._rebuild_plot_grid()
        self.update_display()

    def _apply_ylims(self):
        tmin = self.ymin_input.text().strip()
        tmax = self.ymax_input.text().strip()

        if tmin == "" and tmax == "":
            self.spectrum_ymin = None
            self.spectrum_ymax = None
        else:
            try:
                ymin = float(tmin)
                ymax = float(tmax)
                if ymin >= ymax:
                    return
                self.spectrum_ymin = ymin
                self.spectrum_ymax = ymax
            except ValueError:
                return

        self._spectrum_ylim_version += 1
        self.rendered_cache.clear()
        self.update_display()

    def _update_button_labels(self):
        self.perturb_btn.setText(f"Distortion: {self.selected_perturbation}")
        self.extended_btn.setText(f"Extended plots: {'ON' if self.show_extended else 'OFF'}")

    def _relevant_params(self, plot_type):
        """Return which parameters are relevant for matching this plot type."""
        base_keys = ["A", "B", "m", "B_y", "winding", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y"]
        
        if plot_type in ("specloc_signature_vs_E", "specloc_loggap_vs_E"):
            return base_keys + ["gamma", "W", "kappa"]
        elif plot_type in ("specloc_gap_gamma_E", "specloc_sig_gamma_E"):
            return base_keys + ["W", "kappa"]
        elif plot_type == "specloc_spectrum_vs_gamma":
            return base_keys + ["W", "kappa", "E"]
        elif plot_type == "specloc_spectrum_vs_W":
            return base_keys + ["gamma", "kappa", "E"]
        elif plot_type == "specloc_spectrum_vs_kappa":
            return base_keys + ["gamma", "W", "E"]
        elif plot_type in ("specloc_gap_gamma_W", "specloc_sig_gamma_W"):
            return base_keys + ["kappa", "E"]
        elif plot_type in ("specloc_gap_gamma_kappa", "specloc_sig_gamma_kappa"):
            return base_keys + ["W", "E"]
        elif plot_type in ("specloc_gap_W_kappa", "specloc_sig_W_kappa"):
            return base_keys + ["gamma", "E"]
        else:
            return base_keys

    def _best_record(self, plot_type):
        target = dict(self.current_params)
        relevant_keys = self._relevant_params(plot_type)
        
        candidates = [r for r in self.records 
                     if r.get("plot_type") == plot_type
                     and str(r.get("perturbation_type", "")) == self.selected_perturbation]
        
        if not candidates:
            return None

        def distance(rec):
            d = 0.0
            count = 0
            for k in relevant_keys:
                if k not in target:
                    continue
                rv = rec.get(k)
                if rv is None or math.isnan(rv):
                    continue
                d += (rv - target[k]) ** 2
                count += 1
            return d if count > 0 else float("inf")

        return min(candidates, key=distance)

    def _render_variant(self, plot_type):
        if plot_type == "specloc_spectrum_vs_gamma":
            return f"ylim{self._spectrum_ylim_version}"
        elif plot_type in ("specloc_signature_vs_E", "specloc_loggap_vs_E"):
            E_ref = self.current_params.get("E", 0.0)
            return f"Eref{E_ref:.6f}"
        return "base"

    def _render_images(self, records_needed):
        pairs = [(r["record_id"], r["plot_type"], self._render_variant(r["plot_type"])) for r in records_needed]
        missing = [p for p in pairs if p not in self.rendered_cache]
        
        if not missing:
            return

        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            tf.write("record_id\tplot_type\n")
            for rid, ptype, _ in missing:
                tf.write(f"{rid}\t{ptype}\n")
            req_file = tf.name

        script = Path(__file__).resolve().parents[1] / "simulations" / "data_processing" / "lazy_render.jl"
        cmd = ["julia", str(script), str(self.records_file), str(self.groups_file), str(self.cache_dir), req_file]
        
        env = dict(os.environ)
        env.setdefault("GKSwstype", "png")
        env["SPECLOC_SUPPRESS_GR_WARNINGS"] = "true"
        
        # Pass E reference for signature/loggap plots
        E_ref = self.current_params.get("E")
        if E_ref is not None:
            env["SPECLOC_SIGGAP_E_REF"] = str(E_ref)
        
        # Pass spectrum y-limits
        if self.spectrum_ymin is not None and self.spectrum_ymax is not None:
            env["SPECLOC_SPECTRUM_YMIN"] = str(self.spectrum_ymin)
            env["SPECLOC_SPECTRUM_YMAX"] = str(self.spectrum_ymax)
            env["SPECLOC_LAZY_FORCE_REBUILD"] = "true"

        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
            if proc.returncode == 0:
                for line in proc.stdout.splitlines():
                    parts = line.split("\t")
                    if len(parts) == 3:
                        rid, ptype, path = parts
                        variant = self._render_variant(ptype)
                        self.rendered_cache[(rid, ptype, variant)] = path
        finally:
            Path(req_file).unlink(missing_ok=True)

    def update_display(self):
        target_str = ", ".join(f"{k}={v:.3f}" for k, v in sorted(self.current_params.items()))
        self.setWindowTitle(f"Specloc Viewer - {self.selected_perturbation} - {target_str}")

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
            variant = self._render_variant(ptype)
            path = self.rendered_cache.get((rid, ptype, variant), "")
            
            # Build metadata string with relevant params
            meta_parts = [f"ID={rid}"]
            relevant = self._relevant_params(ptype)
            for k in relevant:
                if k in ["A", "B", "m", "gamma", "W", "kappa", "E"]:
                    v = rec.get(k)
                    if v is not None and not math.isnan(v):
                        meta_parts.append(f"{k}={v:.3f}")
            card["meta"].setText(", ".join(meta_parts))

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
        print("Usage: python3 qwz_specloc_viewer.py <lazy_packet_dir>")
        sys.exit(1)

    packet = Path(sys.argv[1])
    if not packet.exists():
        print(f"Error: {packet} does not exist")
        sys.exit(1)

    app = QApplication(sys.argv[:1])
    viewer = SpeclocViewer(packet)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
