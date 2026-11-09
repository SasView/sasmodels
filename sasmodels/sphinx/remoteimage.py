"""
==================
Intersphinx images
==================

.. role:: code-py(code)
  :language: Python

Insert remote images into documents.

Overrides the ``figure`` and ``image`` directives so that it
retrieves the images from a remote URI if they are not available locally. This
effectively extends intersphinx to support images as well.

Configuration
-------------

To use remote image linking, add ``'sasmodels.sphinx.remoteimage'`` to your ``extensions``
config value, and use these config values to activate linking:

**remote_image_url**

.. container:: confval

    **type**: :code-py:`str`

    **default**: :code-py:`""`

    URI for the remote image directory created for a sphinx/docutils html document.

    This assumes that images are stored in a single directory. Only one remote image
    directory is supported.
"""

from pathlib import Path
from typing import Any

# TODO: Make sphinx a sasmodels dependency
# TODO: Split sphinx extentions into their own pypi packages
from docutils.parsers.rst.directives.images import Image
from sphinx.application import Sphinx
from sphinx.directives.patches import Figure


class RemoteImageMixin:
    """Overridden Image directive to support remote fallback."""
    def run(self):
        env = self.state.document.settings.env
        # The first argument to the image directive is the URI
        uri = self.arguments[0]

        # 1. Determine the absolute path of the image relative to the current doc
        # env.docname is the path to the current .rst file
        doc_dir = Path(env.doc2path(env.docname)).parent
        local_path = Path(doc_dir, uri).resolve()

        # 2. Check if it exists locally
        if not local_path.exists():
            # Use the remote base URL from conf.py
            remote_base = env.config.remote_image_url
            filename = local_path.name

            # Construct remote URL: {remote}/_image/myimg.png
            # Adjust this path logic if your remote structure differs
            self.arguments[0] = f"{remote_base}/{filename}"

        # print(f"Calling {super()}.run() with {self.arguments[0]}")
        return super().run()

class RemoteFigure(RemoteImageMixin, Figure):
    pass

class RemoteImage(RemoteImageMixin, Image):
    pass


def setup(app: Sphinx) -> dict[str, Any]:
    # Add a config value so you can define the remote URL in conf.py
    app.add_config_value(
        name='remote_image_url',
        default='_images', # TODO: use intersphinx inventory uri.parent / "_images"?
        rebuild='html', # all html files change if you modify the remote
        types=[str],
        description='Path to the image directory on an intersphinx document server',
        )

    # Override the standard directives
    app.add_directive('image', RemoteImage, override=True)
    app.add_directive('figure', RemoteFigure, override=True)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
