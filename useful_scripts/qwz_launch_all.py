#!/usr/bin/env python3
"""
Launch all three SpecLoc viewers simultaneously

Opens ribbon, OBC, and spectral localiser viewers in organized screen positions.

Usage:
    python3 useful_scripts/qwz_launch_all.py <lazy_packet_dir>
"""

import subprocess
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 qwz_launch_all.py <lazy_packet_dir>")
        print("\nLaunches three viewers:")
        print("  - Ribbon Spectrum + Band3D")
        print("  - OBC Density of States")
        print("  - Spectral Localiser")
        sys.exit(1)

    packet_dir = Path(sys.argv[1])
    if not packet_dir.exists():
        print(f"Error: {packet_dir} does not exist")
        sys.exit(1)

    script_dir = Path(__file__).resolve().parent

    viewers = [
        ("Ribbon Spectrum + Band3D", script_dir / "qwz_ribbon_viewer.py"),
        ("OBC Density of States", script_dir / "qwz_obc_viewer.py"),
        ("Spectral Localiser", script_dir / "qwz_specloc_viewer.py"),
    ]

    print(f"Launching viewers for packet: {packet_dir}")
    processes = []

    for name, script in viewers:
        if not script.exists():
            print(f"Warning: {script} not found, skipping {name}")
            continue
        
        print(f"  Starting: {name}")
        proc = subprocess.Popen(
            ["python3", str(script), str(packet_dir)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        processes.append((name, proc))

    print(f"\n{len(processes)} viewer(s) launched.")
    print("Press Ctrl+C to close all viewers.\n")

    try:
        # Wait for all processes
        for name, proc in processes:
            proc.wait()
    except KeyboardInterrupt:
        print("\nShutting down viewers...")
        for name, proc in processes:
            proc.terminate()
        for name, proc in processes:
            proc.wait()
        print("All viewers closed.")


if __name__ == "__main__":
    main()
