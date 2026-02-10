#!/usr/bin/env python3
"""
Simple GUI for SASModels Shape Visualizer

Features:
- Lists available shape models (form-factor models with geometry)
- Lets you select a model and plot it with default parameters
- Optional toggle for cross-section subplots

This is a lightweight front-end around `sasmodels_shape_visualizer`.

Note on macOS:
- All matplotlib/Tkinter GUI work is kept on the main thread to avoid
  NSException / Cocoa threading issues.
"""

import os
import sys

# Ensure sasmodels package is importable when running from explore/ directory
_explore_dir = os.path.dirname(os.path.abspath(__file__))
_parent_dir = os.path.dirname(_explore_dir)
if _parent_dir not in sys.path:
    sys.path.insert(0, _parent_dir)

import tkinter as tk
from tkinter import messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  needed to enable 3D projection
from shape_visualizer import SASModelsLoader, SASModelsShapeDetector


class ShapeVisualizerGUI(tk.Tk):
    """Main GUI window for selecting and plotting models."""

    def __init__(self) -> None:
        super().__init__()
        self.title("SASModels Shape Visualizer")

        # Make window reasonably sized (wider to accommodate new layout)
        self.geometry("1200x700")

        # Track last plotted model so we don't replot on the same selection
        self._last_plotted_name: str | None = None
        self._canvas: FigureCanvasTkAgg | None = None
        # Map listbox row index -> index into self.models, or None for category headers
        self._list_index_map: list[int | None] = []

        # UI Elements
        self._create_widgets()

        # Load model list
        self.models = self._load_models()
        self._populate_model_list()

    def _create_widgets(self) -> None:
        """Create all widgets."""
        # Overall grid:
        # - Column 0: model list (full height)
        # - Column 1: description on top, plot below
        # - Row 2: controls
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=0)

        # Left: list of models (full height)
        list_frame = ttk.Frame(self)
        list_frame.grid(row=0, column=0, rowspan=2, sticky="nsew", padx=(10, 5), pady=10)
        list_frame.columnconfigure(0, weight=1)
        list_frame.rowconfigure(1, weight=1)

        ttk.Label(list_frame, text="Shape Models").grid(row=0, column=0, sticky="w")

        self.model_listbox = tk.Listbox(list_frame, exportselection=False)
        self.model_listbox.grid(row=1, column=0, sticky="nsew", pady=(5, 0))

        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.model_listbox.yview)
        scrollbar.grid(row=1, column=1, sticky="ns")
        self.model_listbox.configure(yscrollcommand=scrollbar.set)

        self.model_listbox.bind("<<ListboxSelect>>", self._on_model_select)
        self.model_listbox.bind("<Double-Button-1>", self._on_model_double_click)

        # Right-top: details (model description at the top)
        info_frame = ttk.Frame(self)
        info_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 10), pady=(10, 5))
        info_frame.columnconfigure(0, weight=1)
        info_frame.rowconfigure(0, weight=1)

        self.info_text = tk.Text(info_frame, wrap="word", height=6)
        self.info_text.grid(row=0, column=0, sticky="nsew")
        self.info_text.configure(state="disabled")

        # Right-bottom: plot area (embedded matplotlib figure)
        self.plot_frame = ttk.Frame(self)
        self.plot_frame.grid(row=1, column=1, sticky="nsew", padx=(5, 10), pady=(5, 10))
        self.plot_frame.columnconfigure(0, weight=1)
        self.plot_frame.rowconfigure(0, weight=1)

        # Bottom: buttons and options
        bottom_frame = ttk.Frame(self)
        bottom_frame.grid(row=2, column=0, columnspan=2, sticky="ew", padx=10, pady=(0, 10))
        bottom_frame.columnconfigure(0, weight=1)
        bottom_frame.columnconfigure(1, weight=0)
        bottom_frame.columnconfigure(2, weight=0)

        self.cross_sections_var = tk.BooleanVar(value=True)
        self.wireframe_var = tk.BooleanVar(value=False)

        options_frame = ttk.Frame(bottom_frame)
        options_frame.grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(
            options_frame,
            text="Show cross-sections",
            variable=self.cross_sections_var,
        ).grid(row=0, column=0, padx=(0, 15))
        ttk.Checkbutton(
            options_frame,
            text="Wireframe",
            variable=self.wireframe_var,
        ).grid(row=0, column=1)

        ttk.Button(
            bottom_frame,
            text="Plot Selected",
            command=self._on_plot_clicked,
        ).grid(row=0, column=1, padx=(10, 5))

        ttk.Button(
            bottom_frame,
            text="Close",
            command=self.destroy,
        ).grid(row=0, column=2, padx=(5, 0))

    def _load_models(self):
        """Load available shape models."""
        models = SASModelsLoader.list_available_models()
        shape_models = []
        for name in models:
            info = SASModelsLoader.load_model_info(name)
            if not info:
                continue
            category = info.get("category", "")
            if isinstance(category, str) and category.startswith("shape:"):
                shape_models.append((name, info))
        # Sort by category, then model name
        shape_models.sort(key=lambda x: (x[1].get("category", ""), x[0]))
        return shape_models

    def _populate_model_list(self) -> None:
        """Fill listbox with model names."""
        self.model_listbox.delete(0, tk.END)
        self._list_index_map = []

        current_category: str | None = None
        for idx, (name, info) in enumerate(self.models):
            category = info.get("category", "shape:other")
            if category != current_category:
                current_category = category
                # Insert category header row
                header = f"[{category}]"
                self.model_listbox.insert(tk.END, header)
                self._list_index_map.append(None)
            # Insert model row, indented under category
            # Check if visualization is available
            has_viz = info.get("has_shape_visualization", False)
            if has_viz:
                display_name = f"  {name}"
            else:
                # Add visual indicator for models without visualization
                display_name = f"  ○ {name}"
            self.model_listbox.insert(tk.END, display_name)
            self._list_index_map.append(idx)

    def _get_selected_model(self):
        """Return (name, info) for current selection, or (None, None)."""
        selection = self.model_listbox.curselection()
        if not selection:
            return None, None
        row = selection[0]
        if row >= len(self._list_index_map):
            return None, None
        model_idx = self._list_index_map[row]
        if model_idx is None:
            # Category header selected; treat as no model selected
            return None, None
        return self.models[model_idx]

    def _on_model_select(self, event=None) -> None:  # noqa: ARG002
        """Update info pane when selection changes."""
        name, info = self._get_selected_model()
        self.info_text.configure(state="normal")
        self.info_text.delete("1.0", tk.END)
        if not name:
            self.info_text.insert("1.0", "Select a model to see details.")
        else:
            self._fill_info_text(name, info)
        self.info_text.configure(state="disabled")

        # Automatically plot when a new model is selected
        # Use a short defer to allow listbox state to settle
        self.after(100, self._auto_plot_if_changed)

    def _auto_plot_if_changed(self) -> None:
        """Plot automatically if the selection changed since last plot."""
        name, info = self._get_selected_model()
        if not name:
            return
        if name == self._last_plotted_name:
            return
        # Plot synchronously on main thread (safe for Tk/macos)
        self._plot_model(name, info)

    def _fill_info_text(self, name: str, info: dict) -> None:
        """Fill the info text widget with model metadata."""
        title = info.get("title", "")
        category = info.get("category", "")
        description = info.get("description", "").strip()

        text_lines = []
        text_lines.append(f"Model: {name}")
        if title:
            text_lines.append(f"Title: {title}")
        if category:
            text_lines.append(f"Category: {category}")
        text_lines.append("")

        # Show volume/orientation parameters
        text_lines.append("Key parameters:")
        params = info.get("parameters", [])
        for p in params:
            if len(p) < 5:
                continue
            pname, units, default, _bounds, ptype, *pdesc = p
            if ptype not in ("volume", "orientation"):
                continue
            desc = pdesc[0] if pdesc else ""
            unit_str = f" [{units}]" if units else ""
            text_lines.append(f"  - {pname}{unit_str} ({ptype}), default={default}  {desc}")

        if description:
            text_lines.append("")
            text_lines.append("Description:")
            text_lines.append(description.splitlines()[0])

        self.info_text.insert("1.0", "\n".join(text_lines))

    def _on_model_double_click(self, event=None) -> None:  # noqa: ARG002
        """Double-click to plot."""
        self._on_plot_clicked()

    def _on_plot_clicked(self) -> None:
        """Plot selected model with default parameters."""
        name, info = self._get_selected_model()
        if not name:
            messagebox.showinfo("No selection", "Please select a model first.")
            return
        self._plot_model(name, info)

    def _plot_model(self, name: str, info: dict) -> None:
        """Create visualizer and show plot embedded in the GUI."""
        try:
            visualizer = SASModelsShapeDetector.create_visualizer(info)
            if visualizer is None:
                messagebox.showerror("Visualizer error", f"No visualizer available for model '{name}'.")
                return
            params = visualizer.get_default_params()
            show_cross = self.cross_sections_var.get()
            show_wire = self.wireframe_var.get()

            # Create a new matplotlib Figure
            # Use larger figure size to accommodate the new layout
            fig = Figure(figsize=(10, 7), dpi=100)

            if show_cross:
                # Use GridSpec for flexible layout:
                # - Left column (wider): 3D plot
                # - Right column (narrower): 3 cross-sections stacked vertically
                # Increased hspace to prevent legend overlap between cross-sections
                gs = GridSpec(3, 2, figure=fig, width_ratios=[3, 1], hspace=0.6, wspace=0.3)

                # 3D plot takes entire left column (all 3 rows)
                ax_3d = fig.add_subplot(gs[:, 0], projection="3d")

                # Cross-sections stacked in right column
                ax_xy = fig.add_subplot(gs[0, 1])
                ax_xz = fig.add_subplot(gs[1, 1])
                ax_yz = fig.add_subplot(gs[2, 1])
            else:
                ax_3d = fig.add_subplot(111, projection="3d")
                ax_xy = ax_xz = ax_yz = None  # type: ignore[assignment]

            # Use same mesh logic as CLI visualizer
            mesh_data = visualizer.create_mesh(params)
            visualizer._plot_mesh_components(ax_3d, mesh_data, show_wire)  # type: ignore[attr-defined]
            visualizer._setup_3d_axis(ax_3d, mesh_data, params)  # type: ignore[attr-defined]

            # Remove axes completely, title, and add parameter text
            ax_3d.set_axis_off()  # Hide all axes
            ax_3d.set_title('')  # Remove title (redundant with parameter box)

            # Get volume parameters to display
            volume_params = visualizer.get_volume_params()
            param_text_lines = []

            # Add model name
            model_name = name.replace('_', ' ').title()
            param_text_lines.append(f"{model_name}")
            param_text_lines.append("")  # Empty line

            # Add key parameters
            for param in volume_params[:5]:  # Show up to 5 key parameters
                if param in params:
                    value = params[param]
                    # Format nicely: remove underscores, capitalize, show value
                    param_display = param.replace('_', ' ').title()
                    param_text_lines.append(f"{param_display} = {value:.1f} Å")

            if param_text_lines:
                param_text = '\n'.join(param_text_lines)
                # Add text annotation in upper left corner of the plot
                ax_3d.text2D(0.02, 0.98, param_text, transform=ax_3d.transAxes,
                           fontsize=10, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            if show_cross and hasattr(visualizer, "_plot_cross_sections"):
                visualizer._plot_cross_sections(ax_xy, ax_xz, ax_yz, params)  # type: ignore[attr-defined]

            fig.tight_layout()

            # Replace any existing canvas
            if self._canvas is not None:
                self._canvas.get_tk_widget().destroy()

            self._canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
            self._canvas.draw()
            self._canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

            self._last_plotted_name = name
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("Plot error", f"Error plotting model '{name}':\n{exc}")


def main() -> None:
    """Run the GUI."""
    app = ShapeVisualizerGUI()
    app.mainloop()


if __name__ == "__main__":
    main()


