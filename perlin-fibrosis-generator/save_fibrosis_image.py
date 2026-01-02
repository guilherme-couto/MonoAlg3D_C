import sys
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.animation import FuncAnimation, PillowWriter
import pyvista as pv

def save_image_from_mat(mat_file, output_filename, title=None, dpi=200):
    # --- 1. Load Data ---
    try:
        data = scipy.io.loadmat(mat_file)
    except Exception as e:
        print(f"Error Python: Could not open {mat_file}. Details: {e}")
        return

    if 'presence' in data:
        matrix = data['presence']
    else:
        keys = [k for k in data.keys() if not k.startswith('__')]
        if not keys:
            print("Error Python: No matrix found.")
            return
        matrix = data[keys[0]]

    matrix = np.squeeze(matrix)
    matrix = (matrix > 0).astype(int)

    # Setup Visualization for Matplotlib
    fibro_cmap = ListedColormap([[0.1, 0.5, 0.8], [0.9, 0.5, 0.1]]) # Blue/Orange

    # Setup Visualization for PyVista
    color_healthy = "#1A80CC"  # Blue
    color_fibrosis = "#E6801A" # Orange

    ndim = matrix.ndim

    # =========================================================================
    # 2D PLOT: MATPLOTLIB
    # =========================================================================
    if ndim == 2:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=dpi)
        ax.imshow(matrix, cmap=fibro_cmap, interpolation='nearest', origin='lower')
        ax.axis('off')
        if title: ax.set_title(title, fontsize=12)
        plt.savefig(f"{output_filename}.png", dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        print(f"Success Python: 2D Image saved as {output_filename}.png")

    # =========================================================================
    # 3D PROCESSING: HYBRID (Matplotlib for GIF + PyVista for Renders)
    # =========================================================================
    elif ndim == 3:
        Ny, Nx, Nz = matrix.shape

        # --- 1. GIF Slice-by-Slice ---
        print(f"Generating 3D slice animation...")
        fig_gif, ax_gif = plt.subplots(figsize=(6, 6), dpi=dpi)
        if title: fig_gif.suptitle(title, fontsize=14)
        ax_gif.axis('off')

        im_gif = ax_gif.imshow(matrix[:, :, 0], cmap=fibro_cmap, interpolation='nearest',
                               origin='lower', vmin=0, vmax=1)
        title_txt = ax_gif.set_title(f"Z-Slice: 0 / {Nz-1}")

        def update_frame(k):
            im_gif.set_array(matrix[:, :, k])
            title_txt.set_text(f"Z-Slice: {k} / {Nz-1}")
            return im_gif, title_txt

        anim = FuncAnimation(fig_gif, update_frame, frames=Nz, interval=100, blit=False)
        try:
            anim.save(f"{output_filename}_movie.gif", writer=PillowWriter(fps=10))
            print(f"Success Python: GIF saved.")
        except: pass
        plt.close(fig_gif)

        # --- SETUP PYVISTA ---
        pv.set_plot_theme("document") # White background

        # Create VTK Grid
        grid = pv.ImageData()
        grid.dimensions = np.array(matrix.shape) + 1
        grid.spacing = (1, 1, 1)
        grid.cell_data["values"] = matrix.flatten(order="F")

        # --- 2. ORTHOGONAL SLICES ---
        print("Generating PyVista Orthogonal Slices...")

        # Create cross-sectional slices in the middle of the volume
        slices = grid.slice_orthogonal()

        p = pv.Plotter(off_screen=True, window_size=[1024, 1024])
        p.add_mesh(grid.outline(), color="k")

        # Add the slices with the defined colors
        p.add_mesh(slices, cmap=[color_healthy, color_fibrosis],
                   show_scalar_bar=False, lighting=False)

        p.add_axes()
        p.view_isometric()
        if title: p.add_text(title, font_size=14)
        p.screenshot(f"{output_filename}_ortho_slices.png")
        p.close()
        print(f"Success Python: Ortho Slices saved.")

        # --- 3. VOLUMETRIC RENDER (Fibrosis Only) ---
        print("Generating PyVista Volumetric Render...")

        # Filter only fibrosis (> 0.5)
        # Occlude healthy tissue to see internal structure
        fibrosis_mesh = grid.threshold(0.5)

        p = pv.Plotter(off_screen=True, window_size=[1024, 1024])
        p.add_mesh(grid.outline(), color="k")

        # Render fibrosis mesh with smooth shading for better appearance
        p.add_mesh(fibrosis_mesh,
                   color=color_fibrosis,
                   smooth_shading=True,
                   specular=0.5,     # Brightness
                   diffuse=0.8,      # Base color
                   ambient=0.3)      # Ambient light

        p.add_axes()
        p.view_isometric()
        if title: p.add_text(title, font_size=14)
        p.screenshot(f"{output_filename}_fibrosis.png")
        p.close()
        print(f"Success Python: Volumetric Render saved.")

    else:
        print(f"Error: Unsupported dimensions {ndim}")

    # Cleanup
    try: os.remove(mat_file)
    except: pass

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python save_fibrosis_image.py <input.mat> <output_name> [title]")
        exit(1)

    input_mat = sys.argv[1]
    output = sys.argv[2]
    ttl = sys.argv[3] if len(sys.argv) > 3 else None

    save_image_from_mat(input_mat, output, ttl)
