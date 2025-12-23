#!/usr/bin/env python3
"""
Script to extract visualization code from shape_visualizer.py and add to model files.
"""
import os
import re

# Read shape_visualizer.py
with open('shape_visualizer.py') as f:
    viz_content = f.read()

# Model to visualizer class mapping
MODELS = {
    'cylinder': ('CylinderVisualizer', 250),
    'core_shell_cylinder': ('CoreShellCylinderVisualizer', 721),
    'capped_cylinder': ('CappedCylinderVisualizer', 553),
    'ellipsoid': ('EllipsoidVisualizer', 341),
    'parallelepiped': ('ParallelepipedVisualizer', 430),
    'core_shell_bicelle_elliptical': ('CoreShellBicelleEllipticalVisualizer', 1012),
    'elliptical_cylinder': ('EllipticalCylinderVisualizer', 1189),
    'hollow_cylinder': ('HollowCylinderVisualizer', 1291),
    'pringle': ('PringleVisualizer', 1674),
    'pearl_necklace': ('PearlNecklaceVisualizer', 1427),
    'stacked_disks': ('StackedDisksVisualizer', 1555),
    'flexible_cylinder_elliptical': ('FlexibleCylinderEllipticalVisualizer', 1784),
}

def extract_class_methods(class_name, start_line, content_lines):
    """Extract create_mesh and _plot_cross_sections methods from a class."""
    # Find class start
    class_start = None
    for i in range(start_line - 1, min(start_line + 50, len(content_lines))):
        if f'class {class_name}' in content_lines[i]:
            class_start = i
            break

    if class_start is None:
        return None, None

    # Find class end (next class or end of file)
    class_end = len(content_lines)
    for i in range(class_start + 1, len(content_lines)):
        if content_lines[i].startswith('class ') and i > class_start:
            class_end = i
            break

    class_lines = content_lines[class_start:class_end]
    class_text = '\n'.join(class_lines)

    # Extract create_mesh
    create_mesh_match = re.search(
        r'def create_mesh\(self[^)]*\):.*?(?=\n    def |\nclass |\Z)',
        class_text, re.DOTALL
    )

    # Extract _plot_cross_sections
    plot_match = re.search(
        r'def _plot_cross_sections\(self[^)]*\):.*?(?=\n    def |\nclass |\Z)',
        class_text, re.DOTALL
    )

    if not create_mesh_match or not plot_match:
        return None, None

    create_mesh_body = create_mesh_match.group(0)
    plot_body = plot_match.group(0)

    # Remove 'self' parameter and fix indentation
    create_mesh_body = re.sub(r'def create_mesh\(self[^)]*\):', 'def create_shape_mesh(params, resolution=50):', create_mesh_body)
    create_mesh_body = re.sub(r'\bself\.', '', create_mesh_body)
    # Fix indentation - remove 4 spaces from each line
    create_lines = create_mesh_body.split('\n')
    create_fixed = []
    for line in create_lines:
        if line.startswith('    '):
            create_fixed.append(line[4:])
        else:
            create_fixed.append(line)
    create_mesh_body = '\n'.join(create_fixed)

    plot_body = re.sub(r'def _plot_cross_sections\(self[^)]*\):', 'def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):', plot_body)
    plot_body = re.sub(r'\bself\.', '', plot_body)
    plot_lines = plot_body.split('\n')
    plot_fixed = []
    for line in plot_lines:
        if line.startswith('    '):
            plot_fixed.append(line[4:])
        else:
            plot_fixed.append(line)
    plot_body = '\n'.join(plot_fixed)

    return create_mesh_body, plot_body

# Process each model
content_lines = viz_content.split('\n')
for model_name, (class_name, approx_line) in MODELS.items():
    model_file = f'../sasmodels/models/{model_name}.py'

    if not os.path.exists(model_file):
        print(f"Warning: {model_file} not found")
        continue

    # Check if already has functions
    with open(model_file) as f:
        model_content = f.read()

    if 'def create_shape_mesh' in model_content:
        print(f"Skipping {model_name} - already has visualization functions")
        continue

    # Extract methods
    create_mesh_code, plot_code = extract_class_methods(class_name, approx_line, content_lines)

    if not create_mesh_code or not plot_code:
        print(f"Warning: Could not extract code for {model_name}")
        continue

    # Add import numpy statements
    create_mesh_code = '    import numpy as np\n' + create_mesh_code.split('\n', 1)[1] if '\n' in create_mesh_code else create_mesh_code
    plot_code = '    import numpy as np\n' + plot_code.split('\n', 1)[1] if '\n' in plot_code else plot_code

    # Add docstrings
    create_mesh_code = f'    """Create 3D mesh for {model_name} visualization."""\n' + create_mesh_code
    plot_code = f'    """Plot 2D cross-sections of the {model_name}."""\n' + plot_code

    # Wrap in function definitions
    create_mesh_func = f'def create_shape_mesh(params, resolution=50):\n{create_mesh_code}'
    plot_func = f'def plot_shape_cross_sections(ax_xy, ax_xz, ax_yz, params):\n{plot_code}'

    # Find insertion point
    lines = model_content.split('\n')
    insert_idx = None
    for i, line in enumerate(lines):
        if 'has_shape_visualization = True' in line:
            insert_idx = i + 1
            break

    if insert_idx is None:
        print(f"Warning: Could not find insertion point for {model_name}")
        continue

    # Insert functions
    new_lines = lines[:insert_idx] + [''] + [create_mesh_func] + [''] + [plot_func] + [''] + lines[insert_idx:]
    new_content = '\n'.join(new_lines)

    # Write back
    with open(model_file, 'w') as f:
        f.write(new_content)

    print(f"Added visualization functions to {model_name}")

print("\nDone!")

