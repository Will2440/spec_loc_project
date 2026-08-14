#!/usr/bin/env python3
"""
Index-driven interactive viewer for pre-rendered plots.

Preferred workflow:
1) Run simulations/data_processing/main.jl to create plots and plot_index.tsv
2) Launch this viewer pointed at either the index file or its parent folder

Examples:
    python3 useful_scripts/qwz_plt_viewer.py simulations/data_processing/plots/20260814_120000/plot_index.tsv
    python3 useful_scripts/qwz_plt_viewer.py simulations/data_processing/plots/

Fallback:
    If no index file is found, the viewer falls back to legacy ribbon filename parsing.
"""

import sys
import csv
import math
import re
from pathlib import Path
from collections import defaultdict

from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QSlider,
    QGridLayout,
    QPushButton,
    QScrollArea,
    QSizePolicy,
)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QTimer
from PIL import Image
import io


class PlotViewer(QMainWindow):
    def __init__(self, source_path):
        super().__init__()
        self.source_path = Path(source_path)
        self.current_params = {}
        self.param_ranges = {}
        self.records = []
        self.mode = "index"
        self.plot_files = {}
        self.sliders = {}
        self.labels = {}
        self.plot_cards = {}
        self.available_plot_types = []
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_display)
        self.update_timer.setSingleShot(True)

        self.load_source()
        self.init_ui()

    @staticmethod
    def _try_float(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return float("nan")

    def _find_index_file(self, source):
        if source.is_file() and source.name == "plot_index.tsv":
            return source
        if source.is_dir():
            direct = source / "plot_index.tsv"
            if direct.exists():
                return direct
            matches = sorted(source.rglob("plot_index.tsv"), key=lambda p: p.stat().st_mtime)
            if matches:
                return matches[-1]
        return None

    def load_source(self):
        index_file = self._find_index_file(self.source_path)
        if index_file is not None:
            self.mode = "index"
            self.load_index(index_file)
            return

        if self.source_path.is_dir():
            self.mode = "legacy"
            self.scan_legacy_folder(self.source_path)
            return

        raise RuntimeError(f"Could not resolve source path: {self.source_path}")

    def load_index(self, index_file):
        numeric_keys = {
            "A",
            "B",
            "m",
            "B_y",
            "winding",
            "Lx_ribbon",
            "Lx_obc",
            "Ly_obc",
            "specloc_x",
            "specloc_y",
            "gamma",
            "W",
            "kappa",
            "E",
        }

        self.records = []
        with open(index_file, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rec = dict(row)
                for k in numeric_keys:
                    if k in rec:
                        rec[k] = self._try_float(rec[k])
                self.records.append(rec)

        # Build slider ranges from finite numeric values.
        for k in sorted(numeric_keys):
            vals = sorted({r[k] for r in self.records if k in r and not math.isnan(r[k])})
            if len(vals) > 1:
                self.param_ranges[k] = vals

        self.available_plot_types = sorted({r["plot_type"] for r in self.records})
        print(f"Loaded index: {index_file}")
        print(f"Records: {len(self.records)}")
        print(f"Plot types ({len(self.available_plot_types)}): {self.available_plot_types}")

    def scan_legacy_folder(self, plots_dir):
        pattern = r'ribbon_spectrum_A([\d.-]+)_B([\d.-]+)_m([\d.-]+)_gamma([\d.-]+)_perturbation_\w+_B_y([\d.-]+)_winding_(\d+)_(\w+)\.png'
        plots_dict = defaultdict(dict)

        for png_file in Path(plots_dir).glob("*.png"):
            match = re.match(pattern, png_file.name)
            if not match:
                continue
            A, B, m, gamma, B_y, winding, plot_type = match.groups()
            params = {
                "A": float(A),
                "B": float(B),
                "m": float(m),
                "gamma": float(gamma),
                "B_y": float(B_y),
                "winding": int(winding),
            }
            key = tuple((k, params[k]) for k in sorted(params.keys()))
            plots_dict[key][plot_type] = str(png_file)

        for key in plots_dict.keys():
            pd = dict(key)
            for pname, pval in pd.items():
                if pname not in self.param_ranges:
                    self.param_ranges[pname] = set()
                self.param_ranges[pname].add(pval)
        for pname in list(self.param_ranges.keys()):
            self.param_ranges[pname] = sorted(list(self.param_ranges[pname]))

        self.plot_files = plots_dict
        self.available_plot_types = ["ipr", "berry_curv", "chern_accumulation"]
        print(f"Legacy mode: {len(plots_dict)} parameter combinations in {plots_dir}")
    
    def init_ui(self):
        self.setWindowTitle("SpecLoc Plot Viewer")
        self.setGeometry(80, 80, 1700, 980)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout()

        # Controls column
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        control_layout.setSpacing(10)
        control_layout.setAlignment(Qt.AlignTop)
        control_panel.setMinimumWidth(360)

        for param_name in sorted(self.param_ranges.keys()):
            if not self.param_ranges[param_name]:
                continue
            param_layout = QVBoxLayout()

            label = QLabel(f"{param_name}: {self.param_ranges[param_name][0]:.3f}")
            label.setMinimumWidth(150)
            self.labels[param_name] = label
            param_layout.addWidget(label)

            slider = QSlider(Qt.Horizontal)
            slider.setMinimum(0)
            slider.setMaximum(len(self.param_ranges[param_name]) - 1)
            slider.setValue(0)
            slider.setTickPosition(QSlider.TicksBelow)
            slider.setTickInterval(max(1, len(self.param_ranges[param_name]) // 10))
            slider.valueChanged.connect(lambda value, p=param_name: self.on_slider_move(p, value))

            self.sliders[param_name] = slider
            param_layout.addWidget(slider)
            control_layout.addLayout(param_layout)

        reset_button = QPushButton("Reset")
        reset_button.clicked.connect(self.reset_params)
        control_layout.addWidget(reset_button)
        control_layout.addStretch()
        control_scroll = QScrollArea()
        control_scroll.setWidgetResizable(True)
        control_scroll.setMinimumWidth(340)
        control_scroll.setWidget(control_panel)

        # Plot grid (scrollable)
        plots_container = QWidget()
        self.image_layout = QGridLayout(plots_container)
        self.image_layout.setSpacing(12)

        for i, plot_type in enumerate(self.available_plot_types):
            title = QLabel(plot_type)
            title.setAlignment(Qt.AlignCenter)
            title.setStyleSheet("font-weight: bold;")

            image = QLabel("No image")
            image.setAlignment(Qt.AlignCenter)
            image.setStyleSheet("border: 1px solid #999; background: white;")
            image.setMinimumSize(420, 320)
            image.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

            card = QWidget()
            card_layout = QVBoxLayout(card)
            card_layout.addWidget(title)
            card_layout.addWidget(image)

            self.plot_cards[plot_type] = {"title": title, "image": image}
            r = i // 3
            c = i % 3
            self.image_layout.addWidget(card, r, c)

        plots_scroll = QScrollArea()
        plots_scroll.setWidgetResizable(True)
        plots_scroll.setWidget(plots_container)

        main_layout.addWidget(control_scroll, 0)
        main_layout.addWidget(plots_scroll, 1)

        central_widget.setLayout(main_layout)
        self.update_display()

    def on_slider_move(self, param_name, value):
        param_value = self.param_ranges[param_name][value]
        self.current_params[param_name] = param_value
        self.labels[param_name].setText(f"{param_name}: {param_value:.3f}")
        self.update_timer.start(80)

    def _distance(self, rec, target_params):
        d = 0.0
        used = 0
        for k, v in target_params.items():
            if k not in rec:
                continue
            rv = rec[k]
            if isinstance(rv, float) and math.isnan(rv):
                continue
            d += (float(rv) - float(v)) ** 2
            used += 1
        if used == 0:
            return float("inf")
        return d

    def _best_record_for_type(self, plot_type, target_params):
        candidates = [r for r in self.records if r.get("plot_type") == plot_type]
        if not candidates:
            return None
        return min(candidates, key=lambda r: self._distance(r, target_params))

    def _legacy_best_plot_map(self, target_params):
        best_key = None
        best_distance = float("inf")
        for key in self.plot_files.keys():
            pd = dict(key)
            d = 0.0
            for pname, pval in target_params.items():
                if pname in pd:
                    d += (float(pd[pname]) - float(pval)) ** 2
            if d < best_distance:
                best_distance = d
                best_key = key
        if best_key is None:
            return {}
        return self.plot_files.get(best_key, {})

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

        param_str = ", ".join(f"{k}={v:.3f}" for k, v in sorted(target_params.items()))
        self.setWindowTitle(f"SpecLoc Plot Viewer ({self.mode}) - {param_str}")

        if self.mode == "index":
            for plot_type, card in self.plot_cards.items():
                rec = self._best_record_for_type(plot_type, target_params)
                if rec is None:
                    card["image"].setText("No matching record")
                    continue
                card["title"].setText(plot_type)
                self.display_image(card["image"], rec.get("plot_path", ""))
        else:
            legacy_map = self._legacy_best_plot_map(target_params)
            for plot_type, card in self.plot_cards.items():
                path = legacy_map.get(plot_type, "")
                if path:
                    self.display_image(card["image"], path)
                else:
                    card["image"].setText("No matching image")

    def display_image(self, label, image_path):
        try:
            if not image_path:
                label.setText("No image path")
                return
            img = Image.open(image_path)
            max_width = 1000
            max_height = 760
            img.thumbnail((max_width, max_height), Image.Resampling.LANCZOS)

            img_byte_arr = io.BytesIO()
            img.save(img_byte_arr, format="PNG")
            img_byte_arr.seek(0)

            pixmap = QPixmap()
            pixmap.loadFromData(img_byte_arr.getvalue(), "PNG")
            label.setPixmap(pixmap)
            label.setText("")
        except Exception as e:
            label.setText(f"Error loading image:\n{str(e)}")

    def reset_params(self):
        for param_name in self.sliders:
            self.sliders[param_name].setValue(0)
            self.current_params[param_name] = self.param_ranges[param_name][0]
            self.labels[param_name].setText(f"{param_name}: {self.param_ranges[param_name][0]:.3f}")
        self.update_display()


def main():
    # Resolve source argument.
    if len(sys.argv) > 1:
        source = Path(sys.argv[1])
    else:
        possible_paths = [
            Path("simulations/data_processing/plots"),
            Path("notebook_calcs/plots/ribbon_spectrum/all_for_interactive_plt"),
        ]
        source = None
        for p in possible_paths:
            if p.exists():
                source = p
                break
        if source is None:
            print("Could not auto-detect plot source.")
            print("Usage: python3 useful_scripts/qwz_plt_viewer.py <plot_index.tsv or folder>")
            sys.exit(1)

    app = QApplication(sys.argv[:1])
    viewer = PlotViewer(source)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
