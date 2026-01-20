#!/usr/bin/env python3
"""
Generalized SASModels Shape Visualizer

This module provides a comprehensive framework for visualizing 3D shapes
from any sasmodels model that has geometric parameters. It automatically
detects shape types and generates appropriate visualizations.

Features:
- Automatic shape detection and classification
- Support for all major sasmodels shape categories
- Extensible architecture for adding new shapes
- Interactive 3D visualization with parameter controls
- Cross-section views and parameter comparison
- Export capabilities for high-quality figures
"""

import argparse
import importlib
import os
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


class ShapeVisualizer(ABC):
    """Abstract base class for shape visualizers."""

    def __init__(self, model_info: Dict[str, Any]):
        self.model_info = model_info
        self.name = model_info.get('name', 'unknown')
        self.parameters = model_info.get('parameters', [])
        self.category = model_info.get('category', 'unknown')
        self._model_module = None  # Cache for model module

    def _get_model_module(self):
        """Get the model module, importing it if necessary."""
        if self._model_module is None:
            try:
                module_path = f'sasmodels.models.{self.name}'
                self._model_module = importlib.import_module(module_path)
            except ImportError:
                self._model_module = False  # Mark as failed
        return self._model_module if self._model_module else None

    @abstractmethod
    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create 3D mesh data for the shape."""
        pass

    @abstractmethod
    def get_default_params(self) -> Dict[str, float]:
        """Get default parameter values for the shape."""
        pass

    def get_volume_params(self) -> List[str]:
        """Get list of volume-related parameters."""
        volume_params = []
        for param in self.parameters:
            if len(param) >= 5 and param[4] == "volume":
                volume_params.append(param[0])
        return volume_params

    def get_orientation_params(self) -> List[str]:
        """Get list of orientation parameters."""
        orientation_params = []
        for param in self.parameters:
            if len(param) >= 5 and param[4] == "orientation":
                orientation_params.append(param[0])
        return orientation_params

    def plot_3d(self, params: Dict[str, float] = None,
                save_file: str = None, show_wireframe: bool = False,
                show_cross_sections: bool = True,
                figsize: Tuple[int, int] = (16, 10)) -> None:
        """Create 3D visualization of the shape with optional cross-sections."""
        if params is None:
            params = self.get_default_params()

        try:
            mesh_data = self.create_mesh(params)
        except Exception as e:
            print(f"Error creating mesh for {self.name}: {e}")
            return

        if show_cross_sections:
            # Create subplot layout: 3D view + cross-sections
            fig = plt.figure(figsize=figsize)

            # Main 3D plot (left side, larger)
            ax_3d = fig.add_subplot(2, 3, (1, 4), projection='3d')

            # Cross-section plots (right side)
            ax_xy = fig.add_subplot(2, 3, 2)  # XY plane (top view)
            ax_xz = fig.add_subplot(2, 3, 3)  # XZ plane (side view)
            ax_yz = fig.add_subplot(2, 3, 5)  # YZ plane (front view)

            # Plot 3D shape
            self._plot_mesh_components(ax_3d, mesh_data, show_wireframe)
            self._setup_3d_axis(ax_3d, mesh_data, params)

            # Plot cross-sections
            self._plot_cross_sections(ax_xy, ax_xz, ax_yz, params)

        else:
            # Original single 3D plot
            fig = plt.figure(figsize=figsize)
            ax_3d = fig.add_subplot(111, projection='3d')

            # Plot the shape components
            self._plot_mesh_components(ax_3d, mesh_data, show_wireframe)
            self._setup_3d_axis(ax_3d, mesh_data, params)

        # Add parameter info box
        self._add_parameter_info(fig, params)

        plt.tight_layout()

        if save_file:
            plt.savefig(save_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved as {save_file}")

        plt.show()

    def _setup_3d_axis(self, ax, mesh_data: Dict[str, Any], params: Dict[str, float]):
        """Setup 3D axis labels, title, and limits."""
        # Set labels and title
        ax.set_xlabel('X (Å)', fontsize=12)
        ax.set_ylabel('Y (Å)', fontsize=12)
        ax.set_zlabel('Z (Å)', fontsize=12)

        # Create title with parameters
        title = f'{self.name.replace("_", " ").title()}\n'
        volume_params = self.get_volume_params()
        for param in volume_params[:3]:  # Show first 3 volume parameters
            if param in params:
                title += f'{param} = {params[param]:.1f} Å, '
        title = title.rstrip(', ')
        ax.set_title(title, fontsize=14, pad=20)

        # Set aspect ratio and limits
        self._set_plot_limits(ax, mesh_data, params)

    @abstractmethod
    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot 2D cross-sections of the shape."""
        pass

    def _plot_mesh_components(self, ax, mesh_data: Dict[str, Any], show_wireframe: bool):
        """Plot mesh components on the 3D axis."""
        colors = ['lightblue', 'lightcoral', 'lightgreen', 'lightyellow', 'lightpink']
        color_idx = 0

        for component_name, component_data in mesh_data.items():
            if component_name.startswith('_'):  # Skip metadata
                continue

            if isinstance(component_data, tuple) and len(component_data) == 3:
                x, y, z = component_data
                color = colors[color_idx % len(colors)]

                if show_wireframe:
                    ax.plot_wireframe(x, y, z, alpha=0.6, color=color, linewidth=0.5)
                else:
                    ax.plot_surface(x, y, z, alpha=0.7, color=color,
                                  edgecolor='none', shade=True)
                color_idx += 1

    def _set_plot_limits(self, ax, mesh_data: Dict[str, Any], params: Dict[str, float]):
        """Set appropriate plot limits based on shape dimensions."""
        # Default limits - subclasses can override
        max_dim = 100
        volume_params = self.get_volume_params()
        if volume_params and volume_params[0] in params:
            max_dim = max(max_dim, params[volume_params[0]] * 1.2)

        ax.set_xlim([-max_dim, max_dim])
        ax.set_ylim([-max_dim, max_dim])
        ax.set_zlim([-max_dim, max_dim])

    def _add_parameter_info(self, fig, params: Dict[str, float]):
        """Add parameter information box to the figure."""
        param_text = f'Parameters ({self.category}):\n'
        volume_params = self.get_volume_params()
        for param in volume_params:
            if param in params:
                param_text += f'• {param}: {params[param]:.1f} Å\n'

        plt.figtext(0.02, 0.98, param_text, fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))


class SphereVisualizer(ShapeVisualizer):
    """Visualizer for spherical shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 50)
        phi = np.linspace(0, np.pi, resolution//2)
        theta = np.linspace(0, 2*np.pi, resolution)
        phi_mesh, theta_mesh = np.meshgrid(phi, theta)
        x = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
        y = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
        z = radius * np.cos(phi_mesh)
        return {'sphere': (x, y, z)}

    def get_default_params(self) -> Dict[str, float]:
        # Try to get default from model parameters
        for param in self.parameters:
            if param[0] == 'radius':
                return {'radius': param[2]}
        return {'radius': 50}

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 50)
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x = radius * np.cos(theta)
        circle_y = radius * np.sin(theta)
        ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2)
        ax_xy.set_xlim(-radius*1.2, radius*1.2)
        ax_xy.set_ylim(-radius*1.2, radius*1.2)
        ax_xy.set_xlabel('X (Å)')
        ax_xy.set_ylabel('Y (Å)')
        ax_xy.set_title('XY Cross-section (Top View)')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)
        ax_xz.plot(circle_x, circle_y, 'r-', linewidth=2)
        ax_xz.set_xlim(-radius*1.2, radius*1.2)
        ax_xz.set_ylim(-radius*1.2, radius*1.2)
        ax_xz.set_xlabel('X (Å)')
        ax_xz.set_ylabel('Z (Å)')
        ax_xz.set_title('XZ Cross-section (Side View)')
        ax_xz.set_aspect('equal')
        ax_xz.grid(True, alpha=0.3)
        ax_yz.plot(circle_x, circle_y, 'g-', linewidth=2)
        ax_yz.set_xlim(-radius*1.2, radius*1.2)
        ax_yz.set_ylim(-radius*1.2, radius*1.2)
        ax_yz.set_xlabel('Y (Å)')
        ax_yz.set_ylabel('Z (Å)')
        ax_yz.set_title('YZ Cross-section (Front View)')
        ax_yz.set_aspect('equal')
        ax_yz.grid(True, alpha=0.3)


class CylinderVisualizer(ShapeVisualizer):
    """Visualizer for cylindrical shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)
        length = params.get('length', 400)
        theta = np.linspace(0, 2*np.pi, resolution)
        z = np.linspace(-length/2, length/2, resolution//2)
        theta_mesh, z_mesh = np.meshgrid(theta, z)
        x = radius * np.cos(theta_mesh)
        y = radius * np.sin(theta_mesh)
        r_cap = np.linspace(0, radius, resolution//4)
        theta_cap = np.linspace(0, 2*np.pi, resolution)
        r_cap_mesh, theta_cap_mesh = np.meshgrid(r_cap, theta_cap)
        x_cap = r_cap_mesh * np.cos(theta_cap_mesh)
        y_cap = r_cap_mesh * np.sin(theta_cap_mesh)
        z_cap_top = np.full_like(x_cap, length/2)
        z_cap_bottom = np.full_like(x_cap, -length/2)
        return {
            'cylinder': (x, y, z_mesh),
            'cap_top': (x_cap, y_cap, z_cap_top),
            'cap_bottom': (x_cap, y_cap, z_cap_bottom)
        }

    def get_default_params(self) -> Dict[str, float]:
        # Extract defaults from model parameters
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'length':
                defaults['length'] = param[2]

        # Fallback defaults
        if 'radius' not in defaults:
            defaults['radius'] = 20
        if 'length' not in defaults:
            defaults['length'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)
        length = params.get('length', 400)
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x = radius * np.cos(theta)
        circle_y = radius * np.sin(theta)
        ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2, label='Cylinder')
        ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)
        ax_xy.set_xlim(-radius*1.5, radius*1.5)
        ax_xy.set_ylim(-radius*1.5, radius*1.5)
        ax_xy.set_xlabel('X (Å)')
        ax_xy.set_ylabel('Y (Å)')
        ax_xy.set_title('XY Cross-section (Top View)')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)
        rect_x = [-length/2, -length/2, length/2, length/2, -length/2]
        rect_z = [-radius, radius, radius, -radius, -radius]
        ax_xz.plot(rect_x, rect_z, 'r-', linewidth=2, label='Cylinder')
        ax_xz.fill(rect_x, rect_z, 'lightcoral', alpha=0.3)
        ax_xz.set_xlim(-length/2*1.2, length/2*1.2)
        ax_xz.set_ylim(-radius*1.5, radius*1.5)
        ax_xz.set_xlabel('Z (Å)')
        ax_xz.set_ylabel('X (Å)')
        ax_xz.set_title('XZ Cross-section (Side View)')
        ax_xz.grid(True, alpha=0.3)
        ax_yz.plot(rect_x, rect_z, 'g-', linewidth=2, label='Cylinder')
        ax_yz.fill(rect_x, rect_z, 'lightgreen', alpha=0.3)
        ax_yz.set_xlim(-length/2*1.2, length/2*1.2)
        ax_yz.set_ylim(-radius*1.5, radius*1.5)
        ax_yz.set_xlabel('Z (Å)')
        ax_yz.set_ylabel('Y (Å)')
        ax_yz.set_title('YZ Cross-section (Front View)')
        ax_yz.grid(True, alpha=0.3)


class EllipsoidVisualizer(ShapeVisualizer):
    """Visualizer for ellipsoidal shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius_polar = params.get('radius_polar', 20)
        radius_equatorial = params.get('radius_equatorial', 400)
        phi = np.linspace(0, np.pi, resolution//2)
        theta = np.linspace(0, 2*np.pi, resolution)
        phi_mesh, theta_mesh = np.meshgrid(phi, theta)
        x = radius_equatorial * np.sin(phi_mesh) * np.cos(theta_mesh)
        y = radius_equatorial * np.sin(phi_mesh) * np.sin(theta_mesh)
        z = radius_polar * np.cos(phi_mesh)
        return {'ellipsoid': (x, y, z)}

    def get_default_params(self) -> Dict[str, float]:
        # Extract defaults from model parameters
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius_polar':
                defaults['radius_polar'] = param[2]
            elif param[0] == 'radius_equatorial':
                defaults['radius_equatorial'] = param[2]

        # Fallback defaults
        if 'radius_polar' not in defaults:
            defaults['radius_polar'] = 20
        if 'radius_equatorial' not in defaults:
            defaults['radius_equatorial'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        radius_polar = params.get('radius_polar', 20)
        radius_equatorial = params.get('radius_equatorial', 400)
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x = radius_equatorial * np.cos(theta)
        circle_y = radius_equatorial * np.sin(theta)
        ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2, label='Equatorial')
        ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)
        ax_xy.set_xlim(-radius_equatorial*1.2, radius_equatorial*1.2)
        ax_xy.set_ylim(-radius_equatorial*1.2, radius_equatorial*1.2)
        ax_xy.set_xlabel('X (Å)')
        ax_xy.set_ylabel('Y (Å)')
        ax_xy.set_title('XY Cross-section (Equatorial)')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)
        ellipse_x = radius_equatorial * np.cos(theta)
        ellipse_z = radius_polar * np.sin(theta)
        ax_xz.plot(ellipse_x, ellipse_z, 'r-', linewidth=2, label='Meridional')
        ax_xz.fill(ellipse_x, ellipse_z, 'lightcoral', alpha=0.3)
        ax_xz.set_xlim(-radius_equatorial*1.2, radius_equatorial*1.2)
        ax_xz.set_ylim(-radius_polar*1.2, radius_polar*1.2)
        ax_xz.set_xlabel('X (Å)')
        ax_xz.set_ylabel('Z (Å)')
        ax_xz.set_title('XZ Cross-section (Meridional)')
        ax_xz.set_aspect('equal')
        ax_xz.grid(True, alpha=0.3)
        ax_yz.plot(ellipse_x, ellipse_z, 'g-', linewidth=2, label='Meridional')
        ax_yz.fill(ellipse_x, ellipse_z, 'lightgreen', alpha=0.3)
        ax_yz.set_xlim(-radius_equatorial*1.2, radius_equatorial*1.2)
        ax_yz.set_ylim(-radius_polar*1.2, radius_polar*1.2)
        ax_yz.set_xlabel('Y (Å)')
        ax_yz.set_ylabel('Z (Å)')
        ax_yz.set_title('YZ Cross-section (Meridional)')
        ax_yz.set_aspect('equal')
        ax_yz.grid(True, alpha=0.3)
        ax_xz.annotate('', xy=(-radius_equatorial, 0), xytext=(radius_equatorial, 0),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xz.text(0, -radius_polar*0.3, f'R_eq = {radius_equatorial:.0f} Å', ha='center', fontsize=10)
        ax_xz.annotate('', xy=(0, -radius_polar), xytext=(0, radius_polar),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xz.text(radius_equatorial*0.3, 0, f'R_pol = {radius_polar:.0f} Å', ha='center', fontsize=10, rotation=90)


class ParallelepipedVisualizer(ShapeVisualizer):
    """Visualizer for parallelepiped shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        length_a = params.get('length_a', 35)
        length_b = params.get('length_b', 75)
        length_c = params.get('length_c', 400)
        faces = {}
        y = np.linspace(-length_b/2, length_b/2, resolution//4)
        z = np.linspace(-length_c/2, length_c/2, resolution//2)
        y_mesh, z_mesh = np.meshgrid(y, z)
        faces['front'] = (np.full_like(y_mesh, length_a/2), y_mesh, z_mesh)
        faces['back'] = (np.full_like(y_mesh, -length_a/2), y_mesh, z_mesh)
        x = np.linspace(-length_a/2, length_a/2, resolution//4)
        z = np.linspace(-length_c/2, length_c/2, resolution//2)
        x_mesh, z_mesh = np.meshgrid(x, z)
        faces['right'] = (x_mesh, np.full_like(x_mesh, length_b/2), z_mesh)
        faces['left'] = (x_mesh, np.full_like(x_mesh, -length_b/2), z_mesh)
        x = np.linspace(-length_a/2, length_a/2, resolution//4)
        y = np.linspace(-length_b/2, length_b/2, resolution//4)
        x_mesh, y_mesh = np.meshgrid(x, y)
        faces['top'] = (x_mesh, y_mesh, np.full_like(x_mesh, length_c/2))
        faces['bottom'] = (x_mesh, y_mesh, np.full_like(x_mesh, -length_c/2))
        return faces

    def get_default_params(self) -> Dict[str, float]:
        # Extract defaults from model parameters
        defaults = {}
        for param in self.parameters:
            if param[0] == 'length_a':
                defaults['length_a'] = param[2]
            elif param[0] == 'length_b':
                defaults['length_b'] = param[2]
            elif param[0] == 'length_c':
                defaults['length_c'] = param[2]

        # Fallback defaults
        if 'length_a' not in defaults:
            defaults['length_a'] = 35
        if 'length_b' not in defaults:
            defaults['length_b'] = 75
        if 'length_c' not in defaults:
            defaults['length_c'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        length_a = params.get('length_a', 35)
        length_b = params.get('length_b', 75)
        length_c = params.get('length_c', 400)
        rect_xy_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
        rect_xy_y = [-length_b/2, length_b/2, length_b/2, -length_b/2, -length_b/2]
        ax_xy.plot(rect_xy_x, rect_xy_y, 'b-', linewidth=2, label='A × B face')
        ax_xy.fill(rect_xy_x, rect_xy_y, 'lightblue', alpha=0.3)
        ax_xy.set_xlim(-length_a/2*1.3, length_a/2*1.3)
        ax_xy.set_ylim(-length_b/2*1.3, length_b/2*1.3)
        ax_xy.set_xlabel('X (Å) - Length A')
        ax_xy.set_ylabel('Y (Å) - Length B')
        ax_xy.set_title('XY Cross-section (Top View)')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)
        rect_xz_x = [-length_a/2, -length_a/2, length_a/2, length_a/2, -length_a/2]
        rect_xz_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]
        ax_xz.plot(rect_xz_x, rect_xz_z, 'r-', linewidth=2, label='A × C face')
        ax_xz.fill(rect_xz_x, rect_xz_z, 'lightcoral', alpha=0.3)
        ax_xz.set_xlim(-length_a/2*1.3, length_a/2*1.3)
        ax_xz.set_ylim(-length_c/2*1.3, length_c/2*1.3)
        ax_xz.set_xlabel('X (Å) - Length A')
        ax_xz.set_ylabel('Z (Å) - Length C')
        ax_xz.set_title('XZ Cross-section (Side View)')
        ax_xz.grid(True, alpha=0.3)
        rect_yz_y = [-length_b/2, -length_b/2, length_b/2, length_b/2, -length_b/2]
        rect_yz_z = [-length_c/2, length_c/2, length_c/2, -length_c/2, -length_c/2]
        ax_yz.plot(rect_yz_y, rect_yz_z, 'g-', linewidth=2, label='B × C face')
        ax_yz.fill(rect_yz_y, rect_yz_z, 'lightgreen', alpha=0.3)
        ax_yz.set_xlim(-length_b/2*1.3, length_b/2*1.3)
        ax_yz.set_ylim(-length_c/2*1.3, length_c/2*1.3)
        ax_yz.set_xlabel('Y (Å) - Length B')
        ax_yz.set_ylabel('Z (Å) - Length C')
        ax_yz.set_title('YZ Cross-section (Front View)')
        ax_yz.grid(True, alpha=0.3)
        ax_xy.annotate('', xy=(-length_a/2, -length_b/2*1.4), xytext=(length_a/2, -length_b/2*1.4),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xy.text(0, -length_b/2*1.5, f'A = {length_a:.0f} Å', ha='center', fontsize=10)
        ax_xy.annotate('', xy=(-length_a/2*1.4, -length_b/2), xytext=(-length_a/2*1.4, length_b/2),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xy.text(-length_a/2*1.5, 0, f'B = {length_b:.0f} Å', ha='center', fontsize=10, rotation=90)
        ax_xz.annotate('', xy=(-length_a/2, -length_c/2*1.2), xytext=(length_a/2, -length_c/2*1.2),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xz.text(0, -length_c/2*1.25, f'A = {length_a:.0f} Å', ha='center', fontsize=10)
        ax_xz.annotate('', xy=(-length_a/2*1.2, -length_c/2), xytext=(-length_a/2*1.2, length_c/2),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xz.text(-length_a/2*1.25, 0, f'C = {length_c:.0f} Å', ha='center', fontsize=10, rotation=90)


class CappedCylinderVisualizer(ShapeVisualizer):
    """Visualizer for capped cylinder shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)
        radius_cap = params.get('radius_cap', 25)
        length = params.get('length', 400)
        if radius_cap < radius:
            raise ValueError(f"Cap radius ({radius_cap}) must be >= cylinder radius ({radius})")
        h = np.sqrt(radius_cap**2 - radius**2)
        theta = np.linspace(0, 2*np.pi, resolution)
        z_cyl = np.linspace(-length/2, length/2, resolution//2)
        theta_cyl, z_cyl_mesh = np.meshgrid(theta, z_cyl)
        x_cyl = radius * np.cos(theta_cyl)
        y_cyl = radius * np.sin(theta_cyl)
        phi_max = np.arccos(h / radius_cap)
        phi = np.linspace(0, phi_max, resolution//4)
        phi_mesh, theta_mesh = np.meshgrid(phi, theta)
        x_cap_top = radius_cap * np.sin(phi_mesh) * np.cos(theta_mesh)
        y_cap_top = radius_cap * np.sin(phi_mesh) * np.sin(theta_mesh)
        z_cap_top = length/2 - h + radius_cap * np.cos(phi_mesh)
        x_cap_bottom = radius_cap * np.sin(phi_mesh) * np.cos(theta_mesh)
        y_cap_bottom = radius_cap * np.sin(phi_mesh) * np.sin(theta_mesh)
        z_cap_bottom = -length/2 + h - radius_cap * np.cos(phi_mesh)
        return {
            'cylinder': (x_cyl, y_cyl, z_cyl_mesh),
            'cap_top': (x_cap_top, y_cap_top, z_cap_top),
            'cap_bottom': (x_cap_bottom, y_cap_bottom, z_cap_bottom)
        }

    def get_default_params(self) -> Dict[str, float]:
        # Extract defaults from model parameters
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'radius_cap':
                defaults['radius_cap'] = param[2]
            elif param[0] == 'length':
                defaults['length'] = param[2]

        # Fallback defaults
        if 'radius' not in defaults:
            defaults['radius'] = 20
        if 'radius_cap' not in defaults:
            defaults['radius_cap'] = 25
        if 'length' not in defaults:
            defaults['length'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)
        radius_cap = params.get('radius_cap', 25)
        length = params.get('length', 400)
        if radius_cap < radius:
            return  # Skip if invalid parameters

        h = np.sqrt(radius_cap**2 - radius**2)

        # XY plane (top view) - circle (same as cylinder)
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x = radius * np.cos(theta)
        circle_y = radius * np.sin(theta)

        ax_xy.plot(circle_x, circle_y, 'b-', linewidth=2, label='Cylinder')
        ax_xy.fill(circle_x, circle_y, 'lightblue', alpha=0.3)

        # Show cap outline if significantly larger
        if radius_cap > radius * 1.1:
            cap_circle_x = radius_cap * np.cos(theta)
            cap_circle_y = radius_cap * np.sin(theta)
            ax_xy.plot(cap_circle_x, cap_circle_y, 'r--', linewidth=1, alpha=0.7, label='Cap outline')

        ax_xy.set_xlim(-radius_cap*1.2, radius_cap*1.2)
        ax_xy.set_ylim(-radius_cap*1.2, radius_cap*1.2)
        ax_xy.set_xlabel('X (Å)')
        ax_xy.set_ylabel('Y (Å)')
        ax_xy.set_title('XY Cross-section (Top View)')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)
        ax_xy.legend()

        # XZ plane (side view) - cylinder + caps
        # Cylinder body
        cyl_x = [-length/2, -length/2, length/2, length/2, -length/2]
        cyl_z = [-radius, radius, radius, -radius, -radius]
        ax_xz.plot(cyl_x, cyl_z, 'b-', linewidth=2, label='Cylinder')
        ax_xz.fill(cyl_x, cyl_z, 'lightblue', alpha=0.3)

        # Spherical caps
        cap_angles = np.linspace(0, 2*np.pi, 100)

        # Top cap
        cap_center_top = length/2 - h
        cap_x_top = cap_center_top + radius_cap * np.cos(cap_angles)
        cap_z_top = radius_cap * np.sin(cap_angles)

        # Only show the part that extends beyond cylinder
        mask_top = cap_x_top >= length/2
        ax_xz.plot(cap_x_top[mask_top], cap_z_top[mask_top], 'r-', linewidth=2, label='Caps')
        ax_xz.fill_between(cap_x_top[mask_top], cap_z_top[mask_top], 0, alpha=0.3, color='lightcoral')

        # Bottom cap
        cap_center_bottom = -length/2 + h
        cap_x_bottom = cap_center_bottom + radius_cap * np.cos(cap_angles)
        cap_z_bottom = radius_cap * np.sin(cap_angles)

        mask_bottom = cap_x_bottom <= -length/2
        ax_xz.plot(cap_x_bottom[mask_bottom], cap_z_bottom[mask_bottom], 'r-', linewidth=2)
        ax_xz.fill_between(cap_x_bottom[mask_bottom], cap_z_bottom[mask_bottom], 0, alpha=0.3, color='lightcoral')

        # Mark cap centers
        ax_xz.plot(cap_center_top, 0, 'ro', markersize=6, label='Cap centers')
        ax_xz.plot(cap_center_bottom, 0, 'ro', markersize=6)

        ax_xz.set_xlim((-length/2 - radius_cap*0.5), (length/2 + radius_cap*0.5))
        ax_xz.set_ylim(-radius_cap*1.2, radius_cap*1.2)
        ax_xz.set_xlabel('Z (Å)')
        ax_xz.set_ylabel('X (Å)')
        ax_xz.set_title('XZ Cross-section (Side View)')
        ax_xz.grid(True, alpha=0.3)
        ax_xz.legend()

        # YZ plane (front view) - same as XZ
        ax_yz.plot(cyl_x, cyl_z, 'g-', linewidth=2, label='Cylinder')
        ax_yz.fill(cyl_x, cyl_z, 'lightgreen', alpha=0.3)

        ax_yz.plot(cap_x_top[mask_top], cap_z_top[mask_top], 'orange', linewidth=2, label='Caps')
        ax_yz.fill_between(cap_x_top[mask_top], cap_z_top[mask_top], 0, alpha=0.3, color='moccasin')
        ax_yz.plot(cap_x_bottom[mask_bottom], cap_z_bottom[mask_bottom], 'orange', linewidth=2)
        ax_yz.fill_between(cap_x_bottom[mask_bottom], cap_z_bottom[mask_bottom], 0, alpha=0.3, color='moccasin')

        ax_yz.plot(cap_center_top, 0, 'o', color='orange', markersize=6, label='Cap centers')
        ax_yz.plot(cap_center_bottom, 0, 'o', color='orange', markersize=6)

        ax_yz.set_xlim((-length/2 - radius_cap*0.5), (length/2 + radius_cap*0.5))
        ax_yz.set_ylim(-radius_cap*1.2, radius_cap*1.2)
        ax_yz.set_xlabel('Z (Å)')
        ax_yz.set_ylabel('Y (Å)')
        ax_yz.set_title('YZ Cross-section (Front View)')
        ax_yz.grid(True, alpha=0.3)
        ax_yz.legend()

        # Add dimension annotations
        ax_xz.annotate('', xy=(-length/2, -radius*1.4), xytext=(length/2, -radius*1.4),
                      arrowprops=dict(arrowstyle='<->', color='black'))
        ax_xz.text(0, -radius*1.6, f'L = {length:.0f} Å', ha='center', fontsize=10)

        ax_xz.text(cap_center_top + radius_cap*0.3, radius_cap*0.7, f'R = {radius_cap:.0f} Å',
                  fontsize=10, rotation=45)
        ax_xz.text(-length/4, radius*0.7, f'r = {radius:.0f} Å', fontsize=10)
        ax_xz.text(cap_center_top, -radius*0.3, f'h = {h:.1f} Å', fontsize=10, ha='center')


class CoreShellCylinderVisualizer(ShapeVisualizer):
    """Visualizer for core-shell cylinder shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)
        thickness = params.get('thickness', 20)
        length = params.get('length', 400)

        # Outer dimensions
        outer_radius = radius + thickness
        outer_length = length + 2 * thickness

        # Create core cylinder
        theta = np.linspace(0, 2*np.pi, resolution)
        z_core = np.linspace(-length/2, length/2, resolution//2)
        theta_core, z_core_mesh = np.meshgrid(theta, z_core)
        x_core = radius * np.cos(theta_core)
        y_core = radius * np.sin(theta_core)

        # Create shell cylinder (outer surface)
        z_shell = np.linspace(-outer_length/2, outer_length/2, resolution//2)
        theta_shell, z_shell_mesh = np.meshgrid(theta, z_shell)
        x_shell = outer_radius * np.cos(theta_shell)
        y_shell = outer_radius * np.sin(theta_shell)

        # Create end caps (shell only, as annular disks)
        r_cap_inner = np.linspace(0, radius, resolution//8)
        r_cap_outer = np.linspace(radius, outer_radius, resolution//8)
        theta_cap = np.linspace(0, 2*np.pi, resolution)

        # Inner caps (on core)
        r_inner_mesh, theta_inner_mesh = np.meshgrid(r_cap_inner, theta_cap)
        x_cap_core = r_inner_mesh * np.cos(theta_inner_mesh)
        y_cap_core = r_inner_mesh * np.sin(theta_inner_mesh)
        z_cap_core_top = np.full_like(x_cap_core, length/2)
        z_cap_core_bottom = np.full_like(x_cap_core, -length/2)

        # Outer shell caps (annular rings on ends)
        r_outer_mesh, theta_outer_mesh = np.meshgrid(r_cap_outer, theta_cap)
        x_cap_shell = r_outer_mesh * np.cos(theta_outer_mesh)
        y_cap_shell = r_outer_mesh * np.sin(theta_outer_mesh)
        z_cap_shell_top = np.full_like(x_cap_shell, outer_length/2)
        z_cap_shell_bottom = np.full_like(x_cap_shell, -outer_length/2)

        # Middle shell caps (between core and outer shell)
        r_full = np.linspace(0, outer_radius, resolution//4)
        r_full_mesh, theta_full_mesh = np.meshgrid(r_full, theta_cap)
        x_cap_middle = r_full_mesh * np.cos(theta_full_mesh)
        y_cap_middle = r_full_mesh * np.sin(theta_full_mesh)
        z_cap_middle_top = np.full_like(x_cap_middle, length/2)
        z_cap_middle_bottom = np.full_like(x_cap_middle, -length/2)

        return {
            'core_cylinder': (x_core, y_core, z_core_mesh),
            'shell_cylinder': (x_shell, y_shell, z_shell_mesh),
            'shell_cap_top': (x_cap_middle, y_cap_middle, z_cap_middle_top),
            'shell_cap_bottom': (x_cap_middle, y_cap_middle, z_cap_middle_bottom),
            'end_cap_top': (x_cap_shell, y_cap_shell, z_cap_shell_top),
            'end_cap_bottom': (x_cap_shell, y_cap_shell, z_cap_shell_bottom),
        }

    def get_default_params(self) -> Dict[str, float]:
        # Extract defaults from model parameters
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'thickness':
                defaults['thickness'] = param[2]
            elif param[0] == 'length':
                defaults['length'] = param[2]

        # Fallback defaults
        if 'radius' not in defaults:
            defaults['radius'] = 20
        if 'thickness' not in defaults:
            defaults['thickness'] = 20
        if 'length' not in defaults:
            defaults['length'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)
        thickness = params.get('thickness', 20)
        length = params.get('length', 400)
        outer_radius = radius + thickness
        outer_length = length + 2 * thickness
        theta = np.linspace(0, 2*np.pi, 100)
        core_x = radius * np.cos(theta)
        core_y = radius * np.sin(theta)
        shell_x = outer_radius * np.cos(theta)
        shell_y = outer_radius * np.sin(theta)
        ax_xy.plot(shell_x, shell_y, 'r-', linewidth=2, label='Shell')
        ax_xy.fill(shell_x, shell_y, 'lightcoral', alpha=0.3)
        ax_xy.plot(core_x, core_y, 'b-', linewidth=2, label='Core')
        ax_xy.fill(core_x, core_y, 'lightblue', alpha=0.5)
        ax_xy.set_xlim(-outer_radius*1.3, outer_radius*1.3)
        ax_xy.set_ylim(-outer_radius*1.3, outer_radius*1.3)
        ax_xy.set_xlabel('X (Å)')
        ax_xy.set_ylabel('Y (Å)')
        ax_xy.set_title('XY Cross-section (Top View)')
        ax_xy.set_aspect('equal')
        ax_xy.grid(True, alpha=0.3)
        ax_xy.legend()
        core_rect_x = [-length/2, -length/2, length/2, length/2, -length/2]
        core_rect_z = [-radius, radius, radius, -radius, -radius]
        shell_rect_x = [-outer_length/2, -outer_length/2, outer_length/2, outer_length/2, -outer_length/2]
        shell_rect_z = [-outer_radius, outer_radius, outer_radius, -outer_radius, -outer_radius]
        ax_xz.plot(shell_rect_x, shell_rect_z, 'r-', linewidth=2, label='Shell')
        ax_xz.fill(shell_rect_x, shell_rect_z, 'lightcoral', alpha=0.3)
        ax_xz.plot(core_rect_x, core_rect_z, 'b-', linewidth=2, label='Core')
        ax_xz.fill(core_rect_x, core_rect_z, 'lightblue', alpha=0.5)
        ax_xz.set_xlim(-outer_length/2*1.2, outer_length/2*1.2)
        ax_xz.set_ylim(-outer_radius*1.3, outer_radius*1.3)
        ax_xz.set_xlabel('Z (Å)')
        ax_xz.set_ylabel('X (Å)')
        ax_xz.set_title('XZ Cross-section (Side View)')
        ax_xz.grid(True, alpha=0.3)
        ax_xz.legend()
        ax_xz.annotate('', xy=(-length/2, -outer_radius*1.5), xytext=(length/2, -outer_radius*1.5),
                      arrowprops=dict(arrowstyle='<->', color='blue'))
        ax_xz.text(0, -outer_radius*1.6, f'L = {length:.0f} Å (core)', ha='center', fontsize=9)
        ax_xz.annotate('', xy=(-outer_length/2, -outer_radius*1.8), xytext=(outer_length/2, -outer_radius*1.8),
                      arrowprops=dict(arrowstyle='<->', color='red'))
        ax_xz.text(0, -outer_radius*1.9, f'L+2T = {outer_length:.0f} Å (total)', ha='center', fontsize=9, color='red')
        ax_xz.annotate('', xy=(outer_length/2*0.7, 0), xytext=(outer_length/2*0.7, radius),
                      arrowprops=dict(arrowstyle='<->', color='blue'))
        ax_xz.text(outer_length/2*0.75, radius/2, f'r={radius:.0f}', fontsize=9, color='blue')
        ax_xz.annotate('', xy=(outer_length/2*0.85, 0), xytext=(outer_length/2*0.85, outer_radius),
                      arrowprops=dict(arrowstyle='<->', color='red'))
        ax_xz.text(outer_length/2*0.9, outer_radius/2, f'R+T={outer_radius:.0f}', fontsize=9, color='red')
        ax_yz.plot(shell_rect_x, shell_rect_z, 'g-', linewidth=2, label='Shell')
        ax_yz.fill(shell_rect_x, shell_rect_z, 'lightgreen', alpha=0.3)
        ax_yz.plot(core_rect_x, core_rect_z, 'orange', linewidth=2, label='Core')
        ax_yz.fill(core_rect_x, core_rect_z, 'moccasin', alpha=0.5)
        ax_yz.set_xlim(-outer_length/2*1.2, outer_length/2*1.2)
        ax_yz.set_ylim(-outer_radius*1.3, outer_radius*1.3)
        ax_yz.set_xlabel('Z (Å)')
        ax_yz.set_ylabel('Y (Å)')
        ax_yz.set_title('YZ Cross-section (Front View)')
        ax_yz.grid(True, alpha=0.3)
        ax_yz.legend()


class CoreShellSphereVisualizer(ShapeVisualizer):
    """Visualizer for core-shell sphere shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 60)
        thickness = params.get('thickness', 10)
        phi = np.linspace(0, np.pi, resolution//2)
        theta = np.linspace(0, 2*np.pi, resolution)
        phi_mesh, theta_mesh = np.meshgrid(phi, theta)
        x_core = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
        y_core = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
        z_core = radius * np.cos(phi_mesh)
        shell_radius = radius + thickness
        x_shell = shell_radius * np.sin(phi_mesh) * np.cos(theta_mesh)
        y_shell = shell_radius * np.sin(phi_mesh) * np.sin(theta_mesh)
        z_shell = shell_radius * np.cos(phi_mesh)
        return {
            'core': (x_core, y_core, z_core),
            'shell': (x_shell, y_shell, z_shell)
        }

    def get_default_params(self) -> Dict[str, float]:
        # Try to get defaults from model_info if available
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]  # Default value is at index 2
            elif param[0] == 'thickness':
                defaults['thickness'] = param[2]

        # Fallback to hardcoded defaults if not found
        if 'radius' not in defaults:
            defaults['radius'] = 60
        if 'thickness' not in defaults:
            defaults['thickness'] = 10

        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class CoreShellBicelleEllipticalVisualizer(ShapeVisualizer):
    """Visualizer for core-shell bicelle elliptical shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 30)  # r_minor
        x_core = params.get('x_core', 3)  # r_major/r_minor ratio
        thick_rim = params.get('thick_rim', 8)
        thick_face = params.get('thick_face', 14)
        length = params.get('length', 50)

        r_minor = radius
        r_major = radius * x_core

        # Outer dimensions
        outer_r_minor = r_minor + thick_rim
        outer_r_major = r_major + thick_rim
        outer_length = length + 2 * thick_face

        # Create core elliptical cylinder
        theta = np.linspace(0, 2*np.pi, resolution)
        z_core = np.linspace(-length/2, length/2, resolution//2)
        theta_core, z_core_mesh = np.meshgrid(theta, z_core)
        x_core = r_major * np.cos(theta_core)
        y_core = r_minor * np.sin(theta_core)

        # Create shell elliptical cylinder (outer surface)
        z_shell = np.linspace(-outer_length/2, outer_length/2, resolution//2)
        theta_shell, z_shell_mesh = np.meshgrid(theta, z_shell)
        x_shell = outer_r_major * np.cos(theta_shell)
        y_shell = outer_r_minor * np.sin(theta_shell)

        # Create end caps
        # Core end caps (elliptical)
        u = np.linspace(0, 1, resolution//4)
        theta_cap = np.linspace(0, 2*np.pi, resolution)
        u_mesh, theta_cap_mesh = np.meshgrid(u, theta_cap)

        x_cap_core = u_mesh * r_major * np.cos(theta_cap_mesh)
        y_cap_core = u_mesh * r_minor * np.sin(theta_cap_mesh)
        z_cap_core_top = np.full_like(x_cap_core, length/2)
        z_cap_core_bottom = np.full_like(x_cap_core, -length/2)

        # Shell end caps (elliptical)
        x_cap_shell = u_mesh * outer_r_major * np.cos(theta_cap_mesh)
        y_cap_shell = u_mesh * outer_r_minor * np.sin(theta_cap_mesh)
        z_cap_shell_top = np.full_like(x_cap_shell, outer_length/2)
        z_cap_shell_bottom = np.full_like(x_cap_shell, -outer_length/2)

        return {
            'core_cylinder': (x_core, y_core, z_core_mesh),
            'shell_cylinder': (x_shell, y_shell, z_shell_mesh),
            'core_cap_top': (x_cap_core, y_cap_core, z_cap_core_top),
            'core_cap_bottom': (x_cap_core, y_cap_core, z_cap_core_bottom),
            'shell_cap_top': (x_cap_shell, y_cap_shell, z_cap_shell_top),
            'shell_cap_bottom': (x_cap_shell, y_cap_shell, z_cap_shell_bottom),
        }

    def get_default_params(self) -> Dict[str, float]:
        # Extract defaults from model parameters
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'x_core':
                defaults['x_core'] = param[2]
            elif param[0] == 'thick_rim':
                defaults['thick_rim'] = param[2]
            elif param[0] == 'thick_face':
                defaults['thick_face'] = param[2]
            elif param[0] == 'length':
                defaults['length'] = param[2]

        # Fallback defaults
        if 'radius' not in defaults:
            defaults['radius'] = 30
        if 'x_core' not in defaults:
            defaults['x_core'] = 3
        if 'thick_rim' not in defaults:
            defaults['thick_rim'] = 8
        if 'thick_face' not in defaults:
            defaults['thick_face'] = 14
        if 'length' not in defaults:
            defaults['length'] = 50
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class EllipticalCylinderVisualizer(ShapeVisualizer):
    """Visualizer for elliptical cylinder shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius_minor = params.get('radius_minor', 20)
        axis_ratio = params.get('axis_ratio', 1.5)
        length = params.get('length', 400)

        radius_major = radius_minor * axis_ratio

        # Create elliptical cylinder
        theta = np.linspace(0, 2*np.pi, resolution)
        z = np.linspace(-length/2, length/2, resolution//2)
        theta_mesh, z_mesh = np.meshgrid(theta, z)

        x = radius_major * np.cos(theta_mesh)
        y = radius_minor * np.sin(theta_mesh)

        # Create end caps (elliptical)
        u = np.linspace(0, 1, resolution//4)
        theta_cap = np.linspace(0, 2*np.pi, resolution)
        u_mesh, theta_cap_mesh = np.meshgrid(u, theta_cap)

        x_cap = u_mesh * radius_major * np.cos(theta_cap_mesh)
        y_cap = u_mesh * radius_minor * np.sin(theta_cap_mesh)
        z_cap_top = np.full_like(x_cap, length/2)
        z_cap_bottom = np.full_like(x_cap, -length/2)

        return {
            'cylinder': (x, y, z_mesh),
            'cap_top': (x_cap, y_cap, z_cap_top),
            'cap_bottom': (x_cap, y_cap, z_cap_bottom)
        }

    def get_default_params(self) -> Dict[str, float]:
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius_minor':
                defaults['radius_minor'] = param[2]
            elif param[0] == 'axis_ratio':
                defaults['axis_ratio'] = param[2]
            elif param[0] == 'length':
                defaults['length'] = param[2]

        if 'radius_minor' not in defaults:
            defaults['radius_minor'] = 20
        if 'axis_ratio' not in defaults:
            defaults['axis_ratio'] = 1.5
        if 'length' not in defaults:
            defaults['length'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class HollowCylinderVisualizer(ShapeVisualizer):
    """Visualizer for hollow cylinder shapes."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 20)  # Inner/core radius
        thickness = params.get('thickness', 10)
        length = params.get('length', 400)

        outer_radius = radius + thickness

        # Create inner cylinder surface
        theta = np.linspace(0, 2*np.pi, resolution)
        z = np.linspace(-length/2, length/2, resolution//2)
        theta_mesh, z_mesh = np.meshgrid(theta, z)

        x_inner = radius * np.cos(theta_mesh)
        y_inner = radius * np.sin(theta_mesh)

        # Create outer cylinder surface
        x_outer = outer_radius * np.cos(theta_mesh)
        y_outer = outer_radius * np.sin(theta_mesh)

        # Create end caps (annular disks)
        r_cap = np.linspace(radius, outer_radius, resolution//4)
        theta_cap = np.linspace(0, 2*np.pi, resolution)
        r_cap_mesh, theta_cap_mesh = np.meshgrid(r_cap, theta_cap)

        x_cap = r_cap_mesh * np.cos(theta_cap_mesh)
        y_cap = r_cap_mesh * np.sin(theta_cap_mesh)
        z_cap_top = np.full_like(x_cap, length/2)
        z_cap_bottom = np.full_like(x_cap, -length/2)

        return {
            'inner_cylinder': (x_inner, y_inner, z_mesh),
            'outer_cylinder': (x_outer, y_outer, z_mesh),
            'cap_top': (x_cap, y_cap, z_cap_top),
            'cap_bottom': (x_cap, y_cap, z_cap_bottom)
        }

    def get_default_params(self) -> Dict[str, float]:
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'thickness':
                defaults['thickness'] = param[2]
            elif param[0] == 'length':
                defaults['length'] = param[2]

        if 'radius' not in defaults:
            defaults['radius'] = 20
        if 'thickness' not in defaults:
            defaults['thickness'] = 10
        if 'length' not in defaults:
            defaults['length'] = 400
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class PearlNecklaceVisualizer(ShapeVisualizer):
    """Visualizer for pearl_necklace: chain of spheres connected by a thin string."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 40) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 80.0)
        edge_sep = params.get('edge_sep', 350.0)  # surface-to-surface distance
        thick_string = params.get('thick_string', 2.5)
        num_pearls = int(round(params.get('num_pearls', 3)))
        num_pearls = max(num_pearls, 1)

        # Spacing between pearl centers
        center_step = 2 * radius + edge_sep
        z_positions = [
            (i - (num_pearls - 1) / 2.0) * center_step
            for i in range(num_pearls)
        ]

        # Sphere (pearl) mesh
        phi = np.linspace(0, np.pi, resolution // 2)
        theta = np.linspace(0, 2 * np.pi, resolution)
        phi_mesh, theta_mesh = np.meshgrid(phi, theta)

        pearls = {}
        for i, z0 in enumerate(z_positions):
            x = radius * np.sin(phi_mesh) * np.cos(theta_mesh)
            y = radius * np.sin(phi_mesh) * np.sin(theta_mesh)
            z = radius * np.cos(phi_mesh) + z0
            pearls[f'pearl_{i}'] = (x, y, z)

        # String segments as thin cylinders between neighboring pearls
        string_radius = thick_string / 2.0
        strings = {}
        if num_pearls > 1 and string_radius > 0:
            theta_c = np.linspace(0, 2 * np.pi, resolution)
            for i in range(num_pearls - 1):
                z_start = z_positions[i] + radius
                z_end = z_positions[i + 1] - radius
                z_seg = np.linspace(z_start, z_end, resolution // 2)
                theta_mesh_c, z_seg_mesh = np.meshgrid(theta_c, z_seg)
                x_c = string_radius * np.cos(theta_mesh_c)
                y_c = string_radius * np.sin(theta_mesh_c)
                strings[f'string_{i}'] = (x_c, y_c, z_seg_mesh)

        mesh = {}
        mesh.update(pearls)
        mesh.update(strings)
        return mesh

    def get_default_params(self) -> Dict[str, float]:
        defaults = {}
        for p in self.parameters:
            if p[0] == 'radius':
                defaults['radius'] = p[2]
            elif p[0] == 'edge_sep':
                defaults['edge_sep'] = p[2]
            elif p[0] == 'thick_string':
                defaults['thick_string'] = p[2]
            elif p[0] == 'num_pearls':
                defaults['num_pearls'] = p[2]
        defaults.setdefault('radius', 80.0)
        defaults.setdefault('edge_sep', 350.0)
        defaults.setdefault('thick_string', 2.5)
        defaults.setdefault('num_pearls', 3)
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class StackedDisksVisualizer(ShapeVisualizer):
    """Visualizer for stacked_disks: stack of core+layer disks."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 40) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        thick_core = params.get('thick_core', 10.0)
        thick_layer = params.get('thick_layer', 10.0)
        radius = params.get('radius', 15.0)
        n_stacking = max(int(round(params.get('n_stacking', 1.0))), 1)

        # Period between disk centers (approximate): core + two layers
        d = thick_core + 2 * thick_layer
        total_height = n_stacking * d
        z_centers = [
            (i - (n_stacking - 1) / 2.0) * d
            for i in range(n_stacking)
        ]

        theta = np.linspace(0, 2 * np.pi, resolution)
        z_disk = np.linspace(-d / 2, d / 2, resolution // 3)
        theta_mesh, z_local = np.meshgrid(theta, z_disk)
        x = radius * np.cos(theta_mesh)
        y = radius * np.sin(theta_mesh)

        mesh = {}
        for i, z0 in enumerate(z_centers):
            z = z_local + z0
            mesh[f'disk_{i}'] = (x, y, z)

        return mesh

    def get_default_params(self) -> Dict[str, float]:
        defaults = {}
        for p in self.parameters:
            if p[0] == 'thick_core':
                defaults['thick_core'] = p[2]
            elif p[0] == 'thick_layer':
                defaults['thick_layer'] = p[2]
            elif p[0] == 'radius':
                defaults['radius'] = p[2]
            elif p[0] == 'n_stacking':
                defaults['n_stacking'] = p[2]
        defaults.setdefault('thick_core', 10.0)
        defaults.setdefault('thick_layer', 10.0)
        defaults.setdefault('radius', 15.0)
        defaults.setdefault('n_stacking', 1.0)
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class PringleVisualizer(ShapeVisualizer):
    """Visualizer for pringle (saddle-shaped) disc."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        radius = params.get('radius', 60)
        thickness = params.get('thickness', 10)
        alpha = params.get('alpha', 0.001)
        beta = params.get('beta', 0.02)

        # Create saddle surface (hyperbolic paraboloid approximation)
        r = np.linspace(0, radius, resolution//2)
        theta = np.linspace(0, 2*np.pi, resolution)
        r_mesh, theta_mesh = np.meshgrid(r, theta)

        x = r_mesh * np.cos(theta_mesh)
        y = r_mesh * np.sin(theta_mesh)

        # Saddle shape: z = alpha * x^2 - beta * y^2
        z = alpha * x**2 - beta * y**2

        # Top and bottom surfaces (offset by thickness)
        z_top = z + thickness/2
        z_bottom = z - thickness/2

        return {
            'top_surface': (x, y, z_top),
            'bottom_surface': (x, y, z_bottom)
        }

    def get_default_params(self) -> Dict[str, float]:
        defaults = {}
        for param in self.parameters:
            if param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'thickness':
                defaults['thickness'] = param[2]
            elif param[0] == 'alpha':
                defaults['alpha'] = param[2]
            elif param[0] == 'beta':
                defaults['beta'] = param[2]

        if 'radius' not in defaults:
            defaults['radius'] = 60
        if 'thickness' not in defaults:
            defaults['thickness'] = 10
        if 'alpha' not in defaults:
            defaults['alpha'] = 0.001
        if 'beta' not in defaults:
            defaults['beta'] = 0.02
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class FlexibleCylinderEllipticalVisualizer(ShapeVisualizer):
    """Visualizer for flexible cylinder with elliptical cross-section (simplified as elliptical cylinder)."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        # Fallback to old implementation if model function not available
        length = params.get('length', 1000)
        radius = params.get('radius', 20)
        axis_ratio = params.get('axis_ratio', 1.5)
        kuhn_length = params.get('kuhn_length', 100)

        radius_minor = radius
        radius_major = radius * axis_ratio

        # Show as straight elliptical cylinder with annotation about flexibility
        theta = np.linspace(0, 2*np.pi, resolution)
        z = np.linspace(-length/2, length/2, resolution//2)
        theta_mesh, z_mesh = np.meshgrid(theta, z)

        x = radius_major * np.cos(theta_mesh)
        y = radius_minor * np.sin(theta_mesh)

        return {
            'cylinder': (x, y, z_mesh),
            '_note': f'Simplified view - actual model is flexible with Kuhn length {kuhn_length} Å'
        }

    def get_default_params(self) -> Dict[str, float]:
        defaults = {}
        for param in self.parameters:
            if param[0] == 'length':
                defaults['length'] = param[2]
            elif param[0] == 'radius':
                defaults['radius'] = param[2]
            elif param[0] == 'axis_ratio':
                defaults['axis_ratio'] = param[2]
            elif param[0] == 'kuhn_length':
                defaults['kuhn_length'] = param[2]

        if 'length' not in defaults:
            defaults['length'] = 1000
        if 'radius' not in defaults:
            defaults['radius'] = 20
        if 'axis_ratio' not in defaults:
            defaults['axis_ratio'] = 1.5
        if 'kuhn_length' not in defaults:
            defaults['kuhn_length'] = 100
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        # Fallback to old implementation if model function not available
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class GenericModelVisualizer(ShapeVisualizer):
    """Generic visualizer that delegates to model functions."""

    def create_mesh(self, params: Dict[str, float], resolution: int = 50) -> Dict[str, Any]:
        """Create mesh by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'create_shape_mesh'):
            return model_module.create_shape_mesh(params, resolution)
        raise NotImplementedError(f"Model '{self.name}' does not have create_shape_mesh function")

    def get_default_params(self) -> Dict[str, float]:
        """Extract defaults from model parameters."""
        defaults = {}
        for param in self.parameters:
            if len(param) >= 5 and param[4] == 'volume':
                defaults[param[0]] = param[2]
        return defaults

    def _plot_cross_sections(self, ax_xy, ax_xz, ax_yz, params: Dict[str, float]):
        """Plot cross-sections by calling model function."""
        model_module = self._get_model_module()
        if model_module and hasattr(model_module, 'plot_shape_cross_sections'):
            model_module.plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params)
            return
        raise NotImplementedError(f"Model '{self.name}' does not have plot_shape_cross_sections function")


class BarbellVisualizer(GenericModelVisualizer):
    """Visualizer for barbell shapes."""
    pass


class VesicleVisualizer(GenericModelVisualizer):
    """Visualizer for vesicle (hollow sphere) shapes."""
    pass


class TriaxialEllipsoidVisualizer(GenericModelVisualizer):
    """Visualizer for triaxial ellipsoid shapes."""
    pass


class RectangularPrismVisualizer(GenericModelVisualizer):
    """Visualizer for rectangular prism shapes."""
    pass


class HollowRectangularPrismVisualizer(GenericModelVisualizer):
    """Visualizer for hollow rectangular prism shapes."""
    pass


class LinearPearlsVisualizer(GenericModelVisualizer):
    """Visualizer for linear pearls shapes."""
    pass


class CoreShellEllipsoidVisualizer(GenericModelVisualizer):
    """Visualizer for core-shell ellipsoid shapes."""
    pass


class FuzzySphereVisualizer(GenericModelVisualizer):
    """Visualizer for fuzzy sphere shapes."""
    pass


class FlexibleCylinderVisualizer(GenericModelVisualizer):
    """Visualizer for flexible cylinder shapes."""
    pass


class MultilayerVesicleVisualizer(GenericModelVisualizer):
    """Visualizer for multilayer vesicle shapes."""
    pass


class SuperballVisualizer(GenericModelVisualizer):
    """Visualizer for superball shapes."""
    pass


class CoreShellParallelepipedVisualizer(GenericModelVisualizer):
    """Visualizer for core-shell parallelepiped shapes."""
    pass


class HollowRectangularPrismThinWallsVisualizer(GenericModelVisualizer):
    """Visualizer for hollow rectangular prism with thin walls."""
    pass


class CoreShellBicelleVisualizer(GenericModelVisualizer):
    """Visualizer for core-shell bicelle shapes."""
    pass


class SASModelsShapeDetector:
    """Automatically detect and classify sasmodels shapes."""

    # Shape type mappings (order matters - more specific shapes first!)
    SHAPE_MAPPINGS = {
        'core_shell_bicelle_elliptical': CoreShellBicelleEllipticalVisualizer,  # Very specific
        'flexible_cylinder_elliptical': FlexibleCylinderEllipticalVisualizer,  # Very specific
        'hollow_rectangular_prism_thin_walls': HollowRectangularPrismThinWallsVisualizer,  # Most specific
        'hollow_rectangular_prism': HollowRectangularPrismVisualizer,  # Must come before rectangular_prism
        'core_shell_parallelepiped': CoreShellParallelepipedVisualizer,  # Must come before parallelepiped
        'core_shell_bicelle': CoreShellBicelleVisualizer,  # Must come before cylinder
        'core_shell_ellipsoid': CoreShellEllipsoidVisualizer,  # Must come before ellipsoid
        'triaxial_ellipsoid': TriaxialEllipsoidVisualizer,  # Must come before ellipsoid
        'multilayer_vesicle': MultilayerVesicleVisualizer,  # Must come before vesicle
        'pearl_necklace': PearlNecklaceVisualizer,
        'linear_pearls': LinearPearlsVisualizer,
        'stacked_disks': StackedDisksVisualizer,
        'fuzzy_sphere': FuzzySphereVisualizer,  # Must come before sphere
        'core_shell_sphere': CoreShellSphereVisualizer,  # Must come before 'sphere'
        'vesicle': VesicleVisualizer,  # Must come before sphere
        'superball': SuperballVisualizer,  # Must come before sphere
        'core_shell_cylinder': CoreShellCylinderVisualizer,  # Must come before 'cylinder'
        'flexible_cylinder': FlexibleCylinderVisualizer,  # Must come before 'cylinder'
        'elliptical_cylinder': EllipticalCylinderVisualizer,  # Must come before 'cylinder'
        'hollow_cylinder': HollowCylinderVisualizer,  # Must come before 'cylinder'
        'capped_cylinder': CappedCylinderVisualizer,  # Must come before 'cylinder'
        'barbell': BarbellVisualizer,  # Similar to capped_cylinder
        'rectangular_prism': RectangularPrismVisualizer,
        'pringle': PringleVisualizer,
        'sphere': SphereVisualizer,
        'cylinder': CylinderVisualizer,
        'ellipsoid': EllipsoidVisualizer,
        'parallelepiped': ParallelepipedVisualizer,
    }

    @classmethod
    def detect_shape_type(cls, model_info: Dict[str, Any]) -> str:
        """Detect the shape type from model information."""
        name = model_info.get('name', '').lower()
        category = model_info.get('category', '').lower()
        param_names = [p[0] for p in model_info.get('parameters', [])]

        # Direct name matching (most specific first)
        for shape_name in cls.SHAPE_MAPPINGS:
            if shape_name in name:
                return shape_name

        # Parameter-based detection for specific shapes (before category matching)
        # Check for core-shell cylinder (has radius, thickness, and length)
        if 'radius' in param_names and 'thickness' in param_names and 'length' in param_names:
            return 'core_shell_cylinder'

        # Check for core-shell sphere (has radius and thickness, but no length)
        if 'radius' in param_names and 'thickness' in param_names and 'length' not in param_names:
            return 'core_shell_sphere'

        # Check for capped cylinder
        if 'radius_cap' in param_names:
            return 'capped_cylinder'

        # Category-based detection
        if 'sphere' in category:
            return 'sphere'
        elif 'cylinder' in category:
            return 'cylinder'
        elif 'ellipsoid' in category:
            return 'ellipsoid'
        elif 'parallelepiped' in category:
            return 'parallelepiped'

        # Fallback parameter-based detection
        if 'radius' in param_names and 'length' not in param_names:
            return 'sphere'
        elif 'radius' in param_names and 'length' in param_names:
            return 'cylinder'
        elif 'radius_polar' in param_names and 'radius_equatorial' in param_names:
            return 'ellipsoid'
        elif 'length_a' in param_names and 'length_b' in param_names:
            return 'parallelepiped'

        return 'unknown'

    @classmethod
    def create_visualizer(cls, model_info: Dict[str, Any]) -> Optional[ShapeVisualizer]:
        """Create appropriate visualizer for the model."""
        shape_type = cls.detect_shape_type(model_info)

        if shape_type in cls.SHAPE_MAPPINGS:
            visualizer_class = cls.SHAPE_MAPPINGS[shape_type]
            return visualizer_class(model_info)

        print(f"Warning: No visualizer available for shape type '{shape_type}'")
        return None


class SASModelsLoader:
    """Load sasmodels model information."""

    @staticmethod
    def load_model_info(model_name: str) -> Optional[Dict[str, Any]]:
        """Load model information from sasmodels."""
        try:
            # Import the model (we're already in sasmodels package)
            module_path = f'sasmodels.models.{model_name}'
            model_module = importlib.import_module(module_path)

            # Extract model information
            model_info = {
                'name': getattr(model_module, 'name', model_name),
                'title': getattr(model_module, 'title', ''),
                'description': getattr(model_module, 'description', ''),
                'category': getattr(model_module, 'category', ''),
                'parameters': getattr(model_module, 'parameters', []),
                'has_shape_visualization': getattr(model_module, 'has_shape_visualization', False),
            }

            return model_info

        except ImportError as e:
            print(f"Error loading model '{model_name}': {e}")
            return None

    @staticmethod
    def list_available_models() -> List[str]:
        """List all available sasmodels."""
        try:
            # We're in explore/shape_visualizer.py, so models are in ../sasmodels/models/
            current_dir = os.path.dirname(__file__)
            models_dir = os.path.join(os.path.dirname(current_dir), 'sasmodels', 'models')

            if not os.path.exists(models_dir):
                print("Warning: Could not find sasmodels/models directory")
                return []

            models = []
            for file in os.listdir(models_dir):
                if file.endswith('.py') and not file.startswith('_'):
                    model_name = file[:-3]  # Remove .py extension
                    models.append(model_name)

            return sorted(models)

        except Exception as e:
            print(f"Error listing models: {e}")
            return []


def create_comparison_plot(model_names: List[str], figsize: Tuple[int, int] = (15, 10)):
    """Create a comparison plot of multiple shapes."""
    n_models = len(model_names)
    if n_models == 0:
        return

    # Calculate subplot layout
    cols = min(3, n_models)
    rows = (n_models + cols - 1) // cols

    fig = plt.figure(figsize=figsize)

    for i, model_name in enumerate(model_names):
        model_info = SASModelsLoader.load_model_info(model_name)
        if model_info is None:
            continue

        visualizer = SASModelsShapeDetector.create_visualizer(model_info)
        if visualizer is None:
            continue

        ax = fig.add_subplot(rows, cols, i+1, projection='3d')

        try:
            params = visualizer.get_default_params()
            mesh_data = visualizer.create_mesh(params)
            visualizer._plot_mesh_components(ax, mesh_data, False)

            ax.set_title(f'{model_name.replace("_", " ").title()}', fontsize=12)
            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')

            # Set limits
            visualizer._set_plot_limits(ax, mesh_data, params)

        except Exception as e:
            ax.text(0.5, 0.5, 0.5, f'Error: {str(e)[:50]}...',
                   transform=ax.transAxes, ha='center', va='center')

    plt.tight_layout()
    plt.show()


def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(description='Generalized SASModels Shape Visualizer')
    parser.add_argument('model', nargs='?', help='Model name to visualize')
    parser.add_argument('--list', action='store_true', help='List available models')
    parser.add_argument('--compare', nargs='+', help='Compare multiple models')
    parser.add_argument('--save', type=str, help='Save plot to file')
    parser.add_argument('--wireframe', action='store_true', help='Show wireframe')
    parser.add_argument('--no-cross-sections', action='store_true', help='Hide cross-section plots')
    parser.add_argument('--resolution', type=int, default=50, help='Mesh resolution')

    # Parameter overrides
    parser.add_argument('--params', type=str, help='Parameter overrides as JSON string')

    args = parser.parse_args()

    if args.list:
        print("Available SASModels:")
        print("=" * 50)
        models = SASModelsLoader.list_available_models()

        # Group by category
        categorized = {}
        for model_name in models:
            model_info = SASModelsLoader.load_model_info(model_name)
            if model_info:
                category = model_info.get('category', 'unknown')
                if category not in categorized:
                    categorized[category] = []
                categorized[category].append(model_name)

        for category, model_list in sorted(categorized.items()):
            print(f"\n{category}:")
            for model in sorted(model_list):
                print(f"  {model}")

        return

    if args.compare:
        print(f"Comparing models: {', '.join(args.compare)}")
        create_comparison_plot(args.compare)
        return

    if not args.model:
        print("Please specify a model name or use --list to see available models")
        return

    # Load and visualize single model
    print(f"Loading model: {args.model}")
    model_info = SASModelsLoader.load_model_info(args.model)

    if model_info is None:
        print(f"Could not load model '{args.model}'")
        return

    print(f"Model: {model_info['name']}")
    print(f"Category: {model_info['category']}")
    print(f"Title: {model_info['title']}")

    visualizer = SASModelsShapeDetector.create_visualizer(model_info)

    if visualizer is None:
        print(f"No visualizer available for model '{args.model}'")
        return

    # Get parameters
    params = visualizer.get_default_params()

    # Override parameters if provided
    if args.params:
        try:
            import json
            param_overrides = json.loads(args.params)
            params.update(param_overrides)
            print(f"Parameter overrides: {param_overrides}")
        except json.JSONDecodeError as e:
            print(f"Error parsing parameters: {e}")

    print(f"Using parameters: {params}")

    # Create visualization
    show_cross_sections = not args.no_cross_sections
    visualizer.plot_3d(params, args.save, args.wireframe, show_cross_sections)


def generate_shape_image(model_name, params=None, output_file=None,
                         show_cross_sections=True, show_wireframe=False):
    """
    Generate a PNG image of the model shape.
    
    This is the public API function for sasview to use.
    
    Args:
        model_name: Name of the sasmodels model
        params: Dictionary of parameter values (uses defaults if None)
        output_file: Path to save PNG (returns BytesIO if None)
        show_cross_sections: Whether to include cross-section views
        show_wireframe: Whether to show wireframe instead of solid
        
    Returns:
        BytesIO buffer with PNG data if output_file is None, or None if saved to file
        
    Raises:
        ValueError: If model doesn't support visualization or parameters are invalid
    """
    from io import BytesIO

    # Load model info
    model_info = SASModelsLoader.load_model_info(model_name)
    if model_info is None:
        raise ValueError(f"Could not load model '{model_name}'")

    # Check if visualization is supported
    # Note: This will be checked via model_info.has_shape_visualization once flags are added
    visualizer = SASModelsShapeDetector.create_visualizer(model_info)
    if visualizer is None:
        raise ValueError(f"Model '{model_name}' does not support shape visualization")

    # Get parameters
    if params is None:
        params = visualizer.get_default_params()

    # Create figure
    if show_cross_sections:
        fig = plt.figure(figsize=(16, 10))
        ax_3d = fig.add_subplot(2, 3, (1, 4), projection='3d')
        ax_xy = fig.add_subplot(2, 3, 2)
        ax_xz = fig.add_subplot(2, 3, 3)
        ax_yz = fig.add_subplot(2, 3, 5)
    else:
        fig = plt.figure(figsize=(10, 8))
        ax_3d = fig.add_subplot(111, projection='3d')
        ax_xy = ax_xz = ax_yz = None

    # Create mesh and plot
    try:
        mesh_data = visualizer.create_mesh(params)
        visualizer._plot_mesh_components(ax_3d, mesh_data, show_wireframe)
        visualizer._setup_3d_axis(ax_3d, mesh_data, params)

        if show_cross_sections:
            visualizer._plot_cross_sections(ax_xy, ax_xz, ax_yz, params)

        visualizer._add_parameter_info(fig, params)
        plt.tight_layout()

        # Save or return
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            return None
        else:
            buf = BytesIO()
            plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            return buf

    except Exception as e:
        plt.close(fig)
        raise ValueError(f"Error generating visualization: {e}")


if __name__ == "__main__":
    main()
