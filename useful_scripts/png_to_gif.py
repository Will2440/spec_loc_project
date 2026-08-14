import os
import re
from PIL import Image

def extract_energy(filename):
    """
    Extracts the numerical value of E from filenames like:
    'localiser_position_heatmaps_A1.0_B1.0_m-1.0_E-0.6122448979591837_kappa0.2_gamma0.0.png'
    """
    # Regex looks for 'E' followed by an optional minus sign and numbers/decimals
    match = re.search(r'E(-?\d+\.?\d*(?:[eE][-+]?\d+)?)', filename)
    if match:
        return float(match.group(1))
    return None

def create_gif(image_paths, output_path, duration=500, loop=0, deletion=False, bounce=False):
    """
    Creates a GIF from a list of PNG images.
    """
    if not image_paths:
        raise ValueError("The list of image paths is empty.")
    
    images = [Image.open(img) for img in image_paths]
    images = [img.convert("RGB") if img.mode != "RGB" else img for img in images]
    
    if bounce:
        images += images[::-1]

    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop
    )
    print(f"GIF saved successfully at: {output_path}")

    if deletion:
        for img_path in image_paths:
            try:
                os.remove(img_path)
                print(f"Deleted: {img_path}")
            except OSError as e:
                print(f"Error deleting {img_path}: {e}")

# Main execution
if __name__ == "__main__":
    folder_path = "/Users/Will/Documents/spec_loc_project/plots/perturbed_spatial_scan/perturbed_spatial_scan_m-1.0_gamma2.0/"
    
    # 1. Find all PNG files in the folder
    all_files = [
        f for f in os.listdir(folder_path) 
        if f.endswith(".png")
    ]

    # 2. Extract E and pair with full path
    files_with_energy = []
    for filename in all_files:
        e_val = extract_energy(filename)
        if e_val is not None:
            full_path = os.path.join(folder_path, filename)
            files_with_energy.append((e_val, full_path))

    # 3. Sort files numerically by E value (ascending order)
    files_with_energy.sort(key=lambda x: x[0])

    # Extract ordered paths
    png_files = [file_path for e_val, file_path in files_with_energy]

    print(f"Found {len(png_files)} frames sorted by energy E.")

    # 4. Generate the GIF
    if png_files:
        output_gif = os.path.join(folder_path, "localiser_heatmaps_energy_sorted.gif")
        create_gif(png_files, output_gif, duration=500, loop=0, deletion=False, bounce=False)
    else:
        print("No matching PNG files containing 'E' were found in the specified folder.")